#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-




if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-t", '--linear_tax_table', type=str, help="linear taxon table generated by ncbitax2lin.py")

    args = parser.parse_args()

    import MySQLdb
    import os
    sqlpsw = os.environ['SQLPSW']
    conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                user="root", # your username
                                passwd=sqlpsw, # your password
                                db="blastnr") # name of the data base
    cursor = conn.cursor()

    headers = ["tax_id","superkingdom","phylum","class","order","family","genus",
               "species","cohort","forma","infraclass","infraorder","kingdom","no rank",
               "no rank1","no rank10","no rank11","no rank12","no rank13","no rank14",
               "no rank15","no rank16","no rank17","no rank18","no rank19","no rank2",
               "no rank20","no rank21","no rank22","no rank3","no rank4","no rank5","no rank6",
               "no rank7","no rank8","no rank9","parvorder","species group","species subgroup",
               "species1","subclass","subfamily","subgenus","subkingdom","suborder","subphylum",
               "subspecies","subtribe","superclass","superfamily","superorder","superphylum","tribe","varietas"]

    '''
    0 "tax_id"
    1 "superkingdom"
    2 "phylum",
    3 "class",
    4 "order",
    5 "family",
    6 "genus",
    7 "species",
    8 "cohort",
    9 "forma",
    10 "infraclass",
    11 "infraorder",
    12 "kingdom",
    13 "no rank",
    14 "no rank1",
    15 "no rank10",
    16 "no rank11",
    17 "no rank12",
    18 "no rank13",
    19 "no rank14",
    20 "no rank15",
    21 "no rank16",
    22 "no rank17",
    23 "no rank18",
    24 "no rank19",
    25 "no rank2",
    26 "no rank20",
    27 "no rank21",
    28 "no rank22",
    29 "no rank3",
    30 "no rank4",
    31 "no rank5",
    32 "no rank6",
    33 "no rank7",
    34 "no rank8",
    35 "no rank9",
    36 "parvorder",
    37 "species group",
    38 "species subgroup",
    39 "species1",
    40 "subclass",
    41 "subfamily",
    41 "subgenus",
    42 "subkingdom",
    43 "suborder",
    44 "subphylum",
    45 "subspecies",
    46 "subtribe",
    47 "superclass",
    48 "superfamily",
    49 "superorder",
    50 "superphylum",
    51 "tribe",
    52 "varietas"]


    0 "tax_id"
    1 "superkingdom"
    12 "kingdom",
    50 "superphylum",
    2 "phylum",
    3 "class",
    4 "order",
    5 "family",
    6 "genus",
    37 "species group",
    38 "species subgroup",
    7 "species",
    45 "subspecies",
    ]

    '''

    sql = 'create table blastnr_blastnr_taxonomy(taxon_id INT,' \
          ' superkingdom varchar(400),' \
          ' kingdom varchar(400),' \
          ' superphylum varchar(400),' \
          ' phylum varchar(400),' \
          ' class varchar(400),' \
          ' `order` varchar(400),' \
          ' family varchar(400),' \
          ' genus varchar(400),' \
          ' species_group varchar(400),' \
          ' species_subgroup varchar(400),' \
          ' species varchar(400),' \
          ' subspecies varchar(400),' \
          ' index taxon_id(taxon_id))'

    print sql
    cursor.execute(sql,)
    conn.commit()

    import csv
    import re

    with open(args.linear_tax_table, 'r') as f:
        reader = csv.reader(f, delimiter=',')
        for n, row in enumerate(reader):
            if n == 0:
                continue
            if n % 1000 == 0:
                print '%s' % (n)
            species_edit = re.sub('"', '', row[7])
            species = re.sub("'", "", species_edit)
            sql = 'insert into blastnr_blastnr_taxonomy values (%s, "%s", "%s", ' \
                  ' "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s")' % (row[0],
                                                                                    row[1],
                                                                                    row[12],
                                                                                    row[50],
                                                                                    row[2],
                                                                                    row[3],
                                                                                    row[4],
                                                                                    row[5],
                                                                                    row[6],
                                                                                    row[37],
                                                                                    row[38],
                                                                                    species,
                                                                                    row[45],
                                                                                    )

            cursor.execute(sql,)
        conn.commit()