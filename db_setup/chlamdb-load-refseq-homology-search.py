#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


# Load blast/plast/diamond results into biosql databas
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: april 2019
# ---------------------------------------------------------------------------
import logging 
logging.basicConfig(filename="load_data.log", level=logging.DEBUG)

def load_blastnr_file_into_db(locus_tag2taxon_id,
                              locus_tag2seqfeature_id,
                              locus_tag2bioentry_id,
                              biodb,
                              hash2locus_list,
                              linear_taxonomy,
                              diamond_refseq, 
                              uniref_db):
    import pandas

    '''
    Load tabulated blast results into sql table blastnr_`db_name`
    Ab unique identifier (primary sql key) is attributed to each blast hsp

    :param seqfeature_id2locus_tag: dictionnary with whole data for `biodb` biodatabase
    :param locus_tag2seqfeature_id: dictionnary with whole data for `biodb` biodatabase
    :param db_name: name of the biodatabase
    :param input_blast_files: all input tabulated blast files
    :return: None
    '''
    from ete3 import NCBITaxa
    ncbi = NCBITaxa()

    import time
    import sqlite3
    from chlamdb.biosqldb import manipulate_biosqldb
    sqlpsw = os.environ['SQLPSW']
    from sqlalchemy import create_engine
    
    engine = create_engine(f"mysql://root:{sqlpsw}@127.0.0.1/{biodb}")
    engine_conn = engine.connect()
    conn_mysql = engine.raw_connection()


    sqlite3_conn = sqlite3.connect(diamond_refseq)
    sqlite3_cursor = sqlite3_conn.cursor()
    # engine = create_engine("mysql+pymysql://…")
    
    sql1 = 'attach "%s" as linear_taxonomy' % linear_taxonomy
    sql2 = 'attach "%s" as uniref' % uniref_db
    
    sqlite3_cursor.execute(sql1)
    sqlite3_cursor.execute(sql2)

    sql3 = 'select t1.*,t2.superkingdom, t3.sequence_length, t3.description from diamond_uniref t1 ' \
           ' inner join linear_taxonomy.ncbi_taxonomy t2 on t1.taxon_id=t2.tax_id' \
            ' inner join uniref.uniref t3 on t1.sseqid=t3.accession'
    print(logging.info("Retrieve table"))
    
    # use chunksize iterator
    hash2locus_list= hash2locus_list[["locus_tag", "hash"]].set_index(["hash"])
    for n,df in enumerate(pandas.read_sql(sql3, sqlite3_conn, chunksize=50000)):
        memory_usage = df.memory_usage().sum()
        logging.info(f"Memory df {n} {memory_usage}")
        
        columns = ["hit_number", 
                "query_hash", 
                "subject_accession", 
                "percent_identity", 
                "length", 
                "mismatch", 
                "gaps", 
                "query_start", 
                "query_end", 
                "subject_start", 
                "subject_end", 
                "evalue", 
                "bit_score", 
                "subject_taxid", 
                "subject_kingdom", 
                "subject_length", 
                "subject_title"]
        
        df.columns = columns
        
        logging.info("Join table with locus2hash")
        df_merged = df.set_index(["query_hash"]).join(hash2locus_list)
        memory_usage_join = df_merged.memory_usage().sum()
        logging.info(f"Memory df join {n} {memory_usage_join}")
        
        df_merged["query_taxon_id"] = [locus_tag2taxon_id[locus_tag] for locus_tag in df_merged.locus_tag]
        df_merged["query_bioentry_id"] = [locus_tag2bioentry_id[locus_tag] for locus_tag in df_merged.locus_tag] 
        df_merged["seqfeature_id"]  = [locus_tag2seqfeature_id[locus_tag] for locus_tag in df_merged.locus_tag]
        df_merged["subject_taxid"] = df_merged["subject_taxid"].astype(int)
        nr_hit_taxid = list(set(df_merged["subject_taxid"].to_list()))
        rank = ncbi.get_rank(nr_hit_taxid)
        subject_scientific_names = ncbi.get_taxid_translator(nr_hit_taxid)
        
        df_merged["subject_scientific_name"]  = [f"{subject_scientific_names[subject_taxid]} ({rank[subject_taxid]})" if subject_taxid in rank else '-' for subject_taxid in df_merged.subject_taxid]

        df_merged[["query_taxon_id","query_bioentry_id", 
                "seqfeature_id", "hit_number", "subject_accession", 
                "subject_kingdom", "subject_scientific_name", 
                "subject_taxid", "subject_title", "evalue", 
                "bit_score", "percent_identity", "gaps", "length", 
                "length", "query_start", "query_end", "subject_start", "subject_end", "subject_length"]].to_sql("blastnr_blastnr", engine, if_exists="append", index=False)



def create_sql_plastnr_tables(db_name):
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)
    conn = server.adaptor.conn
    cursor = server.adaptor.cursor
    
    sql_plast = 'CREATE TABLE IF NOT EXISTS blastnr_blastnr (query_taxon_id INT,' \
                ' query_bioentry_id INT,' \
                ' seqfeature_id INT,' \
                ' hit_number int,' \
                ' subject_accession varchar(200),' \
                ' subject_kingdom varchar(200),' \
                ' subject_scientific_name TEXT(2000000), ' \
                ' subject_taxid INT,' \
                ' subject_title VARCHAR(2000),' \
                ' evalue varchar(200),' \
                ' bit_score float,' \
                ' percent_identity float,' \
                ' gaps int,' \
                ' length int,' \
                ' query_start int,' \
                ' query_end int,' \
                ' subject_start int,' \
                ' subject_end int,' \
                ' subject_length int,' \
                ' INDEX query_taxon_id (query_taxon_id),' \
                ' INDEX query_bioentry_id (query_bioentry_id),' \
                ' INDEX hit_number (hit_number),' \
                ' INDEX seqfeature_id (seqfeature_id),' \
                ' INDEX subject_taxid(subject_taxid))'

    try:

        cursor.execute(sql_plast)
        conn.commit()
    except:
        logging.info(sql_plast)
        logging.info('not created')


if __name__ == '__main__':
    import argparse
    from chlamdb.biosqldb import manipulate_biosqldb
    import sys
    import os
    import json
    import numpy
    import re
    import time
    from datetime import datetime
    from multiprocessing import Process
    import  chlamdb_setup_utils

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", '--mysql_database', type=str, help="Biosql biodatabase name")
    parser.add_argument("-t", '--create_tables', action='store_true', help="Create SQL tables")
    parser.add_argument("-f", '--filter_n_hits', action='store_true', help="filter_n_hits (max 100 hits/locus)")
    parser.add_argument("-u", '--hash2locus_tag', type=str, help="Tab separated file with correspondance between sequence hashes and locus tags")
    parser.add_argument("-lt", '--linear_taxonomy', type=str, help="linear_taxonomy.db")
    parser.add_argument("-rd", '--diamond_refseq', type=str, help="diamond_refseq.db")
    parser.add_argument("-ud", '--uniref_db', type=str, help="uniref100.db")
    

    args = parser.parse_args()

    mysql_host = 'localhost'
    mysql_user = 'root'

    mysql_pwd = os.environ['SQLPSW']

    db_name = args.mysql_database

    biodb = args.mysql_database

    server, db = manipulate_biosqldb.load_db(biodb)

    hash2locus_list = chlamdb_setup_utils.get_hash2locus_list(args.hash2locus_tag, as_df=True)

    sys.stdout.write("creating locus_tag2seqfeature_id")
    locus_tag2seqfeature_id = manipulate_biosqldb.locus_tag2seqfeature_id_dict(server, biodb)

    create_sql_plastnr_tables(db_name)

    logging.info('get locus2taxon_id')
    sql = 'select locus_tag, taxon_id from orthology_detail'
    locus_tag2taxon_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    logging.info('get locus2bioentry')
    sql2 = 'select locus_tag,bioentry_id from biodatabase t1 ' \
           ' inner join bioentry as t2 on t1.biodatabase_id=t2.biodatabase_id' \
           ' inner join orthology_detail t3 on t2.accession=t3.accession where t1.name="%s"' % (db_name)

    locus_tag2bioentry_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2,))

    load_blastnr_file_into_db(locus_tag2taxon_id,
                              locus_tag2seqfeature_id,
                              locus_tag2bioentry_id,
                              db_name,
                              hash2locus_list,
                              args.linear_taxonomy,
                              args.diamond_refseq,
                              args.uniref_db)
    
    manipulate_biosqldb.update_config_table(biodb, "BLAST_refseq")