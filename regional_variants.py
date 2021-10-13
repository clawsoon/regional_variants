#!/usr/bin/env python3

# Treeseq files from:
# https://zenodo.org/record/5495535

import sys
import os
import subprocess
import shutil
import collections
import argparse
import json
import sqlite3

import numpy
import tskit
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
plt.rcParams['svg.fonttype'] = 'none'

import populations


def download(url_template, download_folder, tsz):
    sys.stderr.write(f'Downloading {tsz}...\n')
    url = url_template.format(treeseq_file=treeseq_file)
    dst = os.path.join(download_folder, tsz)
    subprocess.run(['wget', '-O', dst, url])


def extract(download_folder, extraction_folder, tsz):
    sys.stderr.write(f'Extracting {tsz}...\n')
    download = os.path.join(download_folder, tsz)
    if download_folder != extraction_folder:
        # tsunzip only extracts to the folder that the treeseq file is in,
        # so if we want it somewhere else we have to copy.
        tsz_path = os.path.join(extraction_folder, tsz)
        shutil.copy(download, tsz_path)
        # No -k means the .tsz file goes away and we don't have to delete
        # it later.
        subprocess.run(['tsunzip', '-v', tsz_path])
    else:
        subprocess.run(['tsunzip', '-v', '-k', download])


def delete(download_folder, extraction_folder, tsz, treeseq_file):
    treeseq_path = os.path.join(extraction_folder, treeseq_file)
    if os.path.exists(treeseq_path):
        sys.stderr.write(f'Deleting {treeseq_path}...\n')
        os.remove(treeseq_path)


def init(database_path):
    sys.stderr.write(f'Initializing {database_path}...\n')

    if os.path.exists(database_path):
        os.remove(database_path)

    database = sqlite3.connect(database_path)
    cursor = database.cursor()

    cursor.execute('''
        CREATE TABLE IF NOT EXISTS individual (
            id INT PRIMARY KEY,
            name TEXT,
            region TEXT,
            all_regions INTEGER,
            some_regions INTEGER,
            one_region INTEGER,
            one_person INTEGER
        );
    ''')

    cursor.execute('''
        CREATE TABLE IF NOT EXISTS num_regions (
            name TEXT PRIMARY KEY,
            number INTEGER
        );
    ''')

    for field in ['all_regions', 'some_regions', 'one_region', 'one_person']:
        cursor.execute('INSERT INTO num_regions (name, number) VALUES (?, 0)', (field,))

    database.commit()
    database.close()


def fill(database_path, treeseq_file):

    sys.stderr.write(f'Filling {database_path}...\n')

    database = sqlite3.connect(database_path)
    database.row_factory = sqlite3.Row
    cursor = database.cursor()

    num_regions = len(set(populations.regions.values()))
    individual_regions = {}

    sys.stderr.write(f'Loading {treeseq_file}\n')
    sys.stderr.flush()

    ts = tskit.load(treeseq_file)

    individuals = {}

    # We have names in ts which are guaranteed to be consistent,
    # and IDs which aren't.  However, the ts IDs should be faster
    # to work with.
    # So we'll use ts IDs from the first ts file we encounter for
    # our SQL IDs.  After that, the mapping may be arbitrary.
    # We need to use names from each ts file to get a mapping
    # from ts IDs to SQL IDs.
    # So the mapping is ts ID -> SQL ID, and we get there via
    # name.

    # Individual names are all different, but population names
    # duplicated between sources, so we use source+name to identify
    # populations.
    for individual in ts.individuals():
        individual_metadata = json.loads(individual.metadata)
        if 'sgdp_id' in individual_metadata:
            source = 'SGDP'
            individual_name = individual_metadata['sgdp_id']
        elif 'sample' in individual_metadata:
            source = 'HGDP'
            individual_name = individual_metadata['sample']
        elif 'individual_id' in individual_metadata:
            source = 'TGP'
            individual_name = individual_metadata['individual_id']

        population = ts.population(ts.node(individual.nodes[0]).population)
        population_metadata = json.loads(population.metadata)
        population_name = population_metadata['name']

        region = populations.regions[(source, population_name)]
        individual_regions[individual.id] = region

        individuals[individual_name] = {'id': individual.id, 'name': individual_name, 'region': region}

        #cursor.execute('INSERT OR IGNORE INTO individual (id, name, region, all_regions, some_regions, one_region, one_person) VALUES (?, ?, ?, 0, 0, 0, 0)', (individual.id, individual_name, region))

    cursor.execute('SELECT * FROM individual')
    rows = cursor.fetchall()
    if rows:
        sql_ids = {individuals[row['name']]['id']: row['id'] for row in rows}
    else:
        sql_ids = {i['id']: i['id'] for i in individuals.values()}
        entries = [(i['id'], i['name'], i['region']) for i in individuals.values()]
        cursor.executemany('INSERT INTO individual (id, name, region, all_regions, some_regions, one_region, one_person) VALUES (?, ?, ?, 0, 0, 0, 0)', entries)


    for vnum, variant in enumerate(ts.variants()):
        if vnum % 1000 == 0:
            sys.stderr.write(f'\rLoaded {vnum}            ')
            sys.stderr.flush()

        for allele_index, allele in enumerate(variant.alleles):

            # Skip ancestral states and missing data.
            if allele_index == 0 or not allele:
                continue

            node_ids = numpy.where(variant.genotypes == allele_index)[0]
            full_nodes = [ts.node(node) for node in node_ids]
            individual_ids = set(node.individual for node in full_nodes)

            if len(individual_ids) == 1:
                column = 'one_person'
            else:
                region_count = len(set(individual_regions[id] for id in individual_ids))
                if region_count == 1:
                    column = 'one_region'
                elif region_count < num_regions:
                    column = 'some_regions'
                else:
                    column = 'all_regions'

            sql = f'UPDATE individual SET {column}={column}+1 WHERE id=?'
            cursor.executemany(sql, [(sql_ids[id],) for id in individual_ids])
            sql = f'UPDATE num_regions SET number=number+1 WHERE name=?'
            cursor.execute(sql, (column,))

    sys.stderr.write(f'\rLoaded {vnum}            \n')
    sys.stderr.flush()

    sys.stderr.write(f'Done\n')
    database.commit()
    database.close()

    fetch(database_path)


def fetch(database_path):

    database = sqlite3.connect(database_path)
    database.row_factory = sqlite3.Row
    cursor = database.cursor()

    cursor.execute('''
        SELECT
            region,
            AVG(all_regions) AS "All regions",
            AVG(some_regions) AS "Some regions",
            AVG(one_region) AS "One region",
            AVG(one_person) AS "One person"
        FROM
            individual
        GROUP BY
            region
        ORDER BY
            region
    ''')

    fields = [d[0] for d in cursor.description[1:]]
    rows = [row for row in cursor]
    for row in rows:
        print(row['region'])
        for field in fields:
            print(field, row[field])

    return fields, rows


def draw(database_path, image_path):

    sys.stderr.write(f'Drawing {database_path}...\n')

    fields, rows = fetch(database_path)

    regions = [row['region'].replace(' ', '\n') for row in rows]
    bottoms = [0 for r in regions]

    fig, ax = plt.subplots()

    for field in fields:
        values = [row[field] for row in rows]
        ax.bar(regions, values, bottom=bottoms, label=field)
        bottoms = [bottoms[i]+values[i] for i in range(len(values))]

    database = sqlite3.connect(database_path)
    cursor = database.cursor()
    individual_count = list(cursor.execute('SELECT COUNT(*) FROM individual'))[0][0]
    variant_count = list(cursor.execute('SELECT SUM(number) FROM num_regions'))[0][0]
    plt.plot([], [], ' ', label=f'Based on {variant_count:,} genetic variants')
    plt.plot([], [], ' ', label=f'from {individual_count:,} fully sequenced individuals')

    ax.set_ylabel('Variant count')
    #ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.0e'))
    ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda value, position: f'{value:.0e}'.replace('+0', '')))
    ax.set_title('''Variant counts in average individual from a given region
by worldwide distribution of genetic variant''')

    ax.legend()
    #'''based on 3,754 fully-sequenced individuals'''

    plt.savefig(image_path)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('commands', nargs='+', choices=['all', 'init', 'download', 'extract', 'fill', 'delete', 'draw'], default='all')
    parser.add_argument('--database-path', default='regional_variants.sqlite')
    parser.add_argument('--image-path', default='regional_variants.svg')
    parser.add_argument('--url-template', default='https://zenodo.org/record/5495535/files/{treeseq_file}.tsz?download=1')
    parser.add_argument('--download-folder', default=os.getcwd())
    parser.add_argument('--extraction-folder', default=os.getcwd())
    parser.add_argument('--treeseq-template', default='hgdp_tgp_sgdp_chr{chromosome}.dated.trees')
    parser.add_argument('--chromosomes', default='1_p,1_q,2_p,2_q,3_p,3_q,4_p,4_q,5_p,5_q,6_p,6_q,7_p,7_q,8_p,8_q,9_p,9_q,10_p,10_q,11_p,11_q,12_p,12_q,13_q,14_q,15_q,16_p,16_q,17_p,17_q,18_p,18_q,19_p,19_q,20_p,20_q,21_q,22_q')

    args = parser.parse_args()

    chromosomes = args.chromosomes.split(',')

    if 'all' in args.commands or 'init' in args.commands:
        init(args.database_path)

    for chromosome in chromosomes:

        treeseq_file = args.treeseq_template.format(chromosome=chromosome)
        tsz = treeseq_file + '.tsz'
        treeseq_path = os.path.join(args.extraction_folder, treeseq_file)

        if 'all' in args.commands or 'download' in args.commands:
            download(args.url_template, args.download_folder, tsz)

        if 'all' in args.commands or 'extract' in args.commands:
            extract(args.download_folder, args.extraction_folder, tsz)

        if 'all' in args.commands or 'fill' in args.commands:
            fill(args.database_path, treeseq_path)

        if 'all' in args.commands or 'delete' in args.commands:
            delete(args.download_folder, args.extraction_folder, tsz, treeseq_file)

    if 'all' in args.commands or 'draw' in args.commands:
        draw(args.database_path, args.image_path)


