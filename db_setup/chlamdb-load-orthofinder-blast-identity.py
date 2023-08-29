#!/usr/bin/env python

import pandas
import os 
from sqlalchemy import create_engine

def parse_sequences_id(sequenceid):
    lst = []
    f = open(sequenceid, "r")
    for row in f:
        data = row.split(": ")
        id = data[0]
        locus_tag = data[1].split(" ")[0]
        lst.append([id, locus_tag])
        
    return pandas.DataFrame(lst, columns = ['id', 'locus_tag']).set_index("id")

def parse_blast_files(blast_files, keep=["qseqid", "sseqid", "pident", "length"]):
    columns = [ "qseqid",
                "sseqid",
                "pident",
                "length",
                "mismatch",
                "gapopen",
                "qstart",
                "qend",
                "sstart",
                "send",
                "evalue",
                "bitscore"]
    df_list = []
    for blast_file in blast_files:
        df = pandas.read_csv(blast_file, sep='\t', header=None, names=columns)[keep]
        df_list.append(df)
    df_merged = pandas.concat(df_list)
    return df_merged

def parse_orthogroups(orthofinder_orthogroups):
    # Orthogroups.txt
    # OG0000000: AOM43_RS00760 AOM43_RS02410 AOM43_RS0
    f = open(orthofinder_orthogroups, "r")
    og2locus_list = {}
    for row in f:
        data = row.strip().split(": ")
        og_id = data[0]
        locus_list = data[1].split(" ")
        og2locus_list[og_id] = locus_list
    return og2locus_list    

def create_table(biodb):
    sqlpsw = os.environ['SQLPSW']
    engine = create_engine(f"mysql://root:{sqlpsw}@127.0.0.1/{biodb}")
    conn_mysql = engine.connect()
    sql = 'CREATE TABLE if not exists orthology_identity_v2 (' \
                ' seqfeature_id_a INTEGER,' \
                ' seqfeature_id_b INTEGER,' \
                ' pident FLOAT(5) ,' \
                ' length INTEGER , ' \
                ' n_hsps INTEGER )'
    conn_mysql.execute(sql)



def get_parwise_id(orthofinder_orthogroups, 
                   blast_files,
                   sequenceid, 
                   biodb):
    import itertools
    
    sqlpsw = os.environ['SQLPSW']
    engine = create_engine(f"mysql://root:{sqlpsw}@127.0.0.1/{biodb}")
    conn_mysql = engine.connect()
    print("parse_blast_files")
    blast_df = parse_blast_files(blast_files)
    print("parse_orthogroups")
    og2locus_list = parse_orthogroups(orthofinder_orthogroups)
    print("parse_sequences_id")
    id_locus_tag_df = parse_sequences_id(sequenceid)
    print("merge with orthofinder ids with locus_tag df")
    #print(id_locus_tag_df.head())
    blast_df = blast_df.set_index("qseqid").join(id_locus_tag_df.rename(columns={"locus_tag":"locus_a"})).reset_index(drop=True).set_index("sseqid").join(id_locus_tag_df.rename(columns={"locus_tag":"locus_b"})).reset_index(drop=True)
    #print(blast_df.head())
    sql = 'select locus_tag,seqfeature_id from custom_tables_locus2seqfeature_id'
    locus2_seqfeature_id_df = pandas.read_sql(sql, conn_mysql).set_index("locus_tag")
    #print("Merge with seafeature ids df")
    blast_df = blast_df.set_index("locus_a").join(locus2_seqfeature_id_df.rename(columns={"seqfeature_id": "seqfeature_id_a"})) #.reset_index(col_fill="locus_a")
    blast_df.index.name = "locus_a"
    #print(blast_df.head())
    blast_df = blast_df.reset_index().set_index("locus_b").join(locus2_seqfeature_id_df.rename(columns={"seqfeature_id": "seqfeature_id_b"})) #.reset_index(col_fill="locus_b")
    blast_df.index.name = "locus_b"
    #print(blast_df.head())
    blast_df = blast_df.reset_index()[["locus_a", "locus_b","seqfeature_id_a", "seqfeature_id_b", "pident", "length"]].set_index(["locus_a", "locus_b"])
    #print("blast_df", blast_df.head())
    for n,og in enumerate(og2locus_list):
        if n % 100 == 0:
            print(f"{n} / {len(og2locus_list)}")
        locus_list = og2locus_list[og]
        if len(locus_list) == 1:
            continue
        print(len(og), og)
        pairs = [i for i in itertools.permutations(locus_list, 2)]
        df_subset = blast_df.loc[pairs]
        n_hsps = df_subset.groupby(["seqfeature_id_a", "seqfeature_id_b"]).count()[["pident"]]
        n_hsps.columns = ['n_hsps']
        top_hsp = df_subset.groupby(["seqfeature_id_a", "seqfeature_id_b"]).head(1).reset_index(drop=True).set_index(["seqfeature_id_a", "seqfeature_id_b"])
        top_hsp = top_hsp.join(n_hsps)
        top_hsp.to_sql("orthology_identity_v2", conn_mysql, if_exists="append")


if __name__ == '__main__':
    import argparse
    import logging
    logging.basicConfig(filename="ortho_load.log", level=logging.DEBUG)
    from chlamdb.biosqldb import biosql_own_sql_tables
    from chlamdb.biosqldb import get_locus2seqfeature_table
    from chlamdb.biosqldb import manipulate_biosqldb

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", '--db_name', type=str, help="db name")
    parser.add_argument("-b", '--blast_files', nargs="+", help="blast_files")
    parser.add_argument("-s", '--sequences', type=str, help="SequenceIDs.txt")
    parser.add_argument("-o", '--ortho', type=str, help="Orthogroups.txt")
    
    args = parser.parse_args()
    create_table(args.db_name)
    get_parwise_id(args.ortho,
                   args.blast_files,
                   args.sequences,
                   args.db_name)