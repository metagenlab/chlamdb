#! /usr/bin/env python


import sys
from Bio import GenBank
from Bio import Entrez
from Bio import SeqIO
from BioSQL import BioSeqDatabase
from chlamdb.biosqldb.manipulate_biosqldb import load_db
from chlamdb.biosqldb.manipulate_biosqldb import query_yes_no

def import_phylo(phylo_list, biodb):
    from chlamdb.biosqldb import manipulate_biosqldb
    from chlamdb.biosqldb import biosql_own_sql_tables
    import ete3
    import re
    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'create database if not exists biosqldb_phylogenies'
    server.adaptor.execute(sql,)

    sql = 'create table IF NOT EXISTS biosqldb_phylogenies.%s (orthogroup varchar(100), phylogeny longtext);' % biodb
    server.adaptor.execute(sql,)

    locuslag2orthogroup = biosql_own_sql_tables.locus_tag2orthogroup(biodb)
    l = len(phylo_list)
    for n, phylo in enumerate(phylo_list):
        print ("%s/%s" % (n, l))
        t = ete3.Tree(phylo)
        leaves = [i for i in t.iter_leaves()]
        orthogroup = locuslag2orthogroup[leaves[0].name]
        #print t.write()
        sql = 'insert into biosqldb_phylogenies.%s values ("%s", "%s");' % (biodb, orthogroup, t.write())
        #print sql
        server.adaptor.execute(sql,)
    server.commit()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", '--trees', type=str, help="tree files", nargs='+')
    parser.add_argument("-d", '--db_name', type=str, help="db name", required=True)

    args = parser.parse_args()

    import_phylo(args.trees, args.db_name)
