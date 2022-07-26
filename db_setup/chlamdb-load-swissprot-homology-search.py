#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


# Load blast results into biosql databas
# Tablulated blast files: qgi qacc sgi sacc sscinames sskingdoms staxids evalue nident pident positive gaps length qstart qend qcovs sstart send sstrand stitle
# Genbank files of query proteins must be leaded in the biosql database already (used to fetch locus tag, taxon_id, organism name, ...)
# Create 3 new tables (if they do not exit): - blastnr_taxonomy with txaon_id and corresponding full taxonomic path
#                                            todo: use taxonomy implemented in biosqldb shema instead: currently fetch taxonomic path from ncbi and add it to blast_nr_taxonomy
#                                           - blast_nr_`biodb_name`: contain data foreach hit
#                                           - blast_nr_taxonomy_`biodb_name`: contain taxonomical data oreach hit (as one hit can now include multiple taxons/MULTISPECIES hits)
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: april 2015
# ---------------------------------------------------------------------------

'''

    #### general data ####

    hit_number
    taxon_id
    locus_tag
    organism
    orthogroup

    #### NCBI taxonomical ranks ####

    query_gi
    subject_gi
    subject_scientific_name
    subject_taxid
    evalue
    n_identical
    percent_identity
    positive
    gaps
    length
    query_start
    query_end
    query_cov
    subject_start
    subject_end]
    subject_strand
    subject_title

    no_rank
    superkingdom
    kingdom
    subkingdom
    superphylum
    phylum
    subphylum
    superclass
    class
    subclass
    superorder
    order
    suborder
    superfamily
    family
    subfamily
    genus
    subgenus
    species
    species_subgroup
    species_group
    subspecies

    tribe
    infraorder
    subtribe
    forma
    infraclass
    varietas
    parvorder
'''

def _chunks(l, n):
    return [l[i:i+n] for i in range(0, len(l), n)]


def get_swissprot_annotation(accession_list):

    import urllib
    from urllib.error import URLError

    link = f'https://rest.uniprot.org/uniprotkb/search?query=accession:%s&fields=accession,organism_id,annotation_score,protein_name,gene_names,organism_name&format=tsv'  % ('+OR+accession:'.join(accession_list))
    link = f'https://rest.uniprot.org/uniprotkb/search?query=accession:%s&fields=accession,annotation_score&format=tsv&size={len(accession_list)}'  % ('+OR+accession:'.join(accession_list))
    #link = link.replace(' ', '%20')

    req = urllib.request.Request(link)

    page = urllib.request.urlopen(req)
    data = page.read().decode('utf-8').split('\n')
    
    rows = [i.rstrip().split('\t') for i in data]
    accession2score = {}
    for row in rows:
        if row != [''] and row[0] != 'Entry':
            accession2score[row[0]] = row[1]
    return accession2score
    #except:
    #    return (False)

def deleted_uniprot2new_identical_sequence(accession):
    import requests, sys
    import json
    requestURL = f"https://www.ebi.ac.uk/proteins/api/uniparc/accession/{accession}"
    r = requests.get(requestURL, headers={ "Accept" : "application/json"})
    print(r)
    if not r.ok:
      r.raise_for_status()
      sys.exit()
    responseBody = r.text
    data = json.loads(responseBody)
    print(data)
    active_entries = [i for i in data["dbReference"] if (i['active'] == 'Y')]
    inactive_entries = [i  for i in data["dbReference"] if (i['active'] != 'Y')]
    if len(active_entries) > 0:
        uniprot_trembl = [i for i in active_entries if (i['type'] == 'UniProtKB/TrEMBL')]
        uniprot_swissprot = [i for i in active_entries if (i['type'] == 'UniProtKB/Swiss-Prot')]
        if len(uniprot_swissprot) > 0:
            return uniprot_swissprot[0]['id']
        elif len(uniprot_trembl) > 0:
            return uniprot_trembl[0]['id']
        else:
            return False
    return False




def load_blastswissprot_file_into_db(locus_tag2taxon_id,
                                locus_tag2seqfeature_id,
                                locus_tag2bioentry_id,
                                input_blast_files,
                                biodb,
                                swissprot_sqlite,
                                hash2locus_list):

    from chlamdb.biosqldb import plastnr2sqltable
    import MySQLdb
    from chlamdb.biosqldb import manipulate_biosqldb
    import time
    import sqlite3
    import pandas
    server, db = manipulate_biosqldb.load_db(biodb)
    conn = server.adaptor.conn
    cursor = server.adaptor.cursor
    print("sqlite DB", swissprot_sqlite)
    sqlite3_conn = sqlite3.connect(swissprot_sqlite)
    sqlite3_cursor = sqlite3_conn.cursor()

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

    n_file = 0
    sp_accessions = []
    df_list = []
    for one_blast_file in input_blast_files:
        n_file +=1
        with open(one_blast_file, 'r') as f:
            df = pandas.read_csv(one_blast_file, sep='\t', header=None, names=columns)
            # keep only top 100 hits
            df = df.groupby(["qseqid"]).head(100)
            df_list.append(df)
            df["sseqid"] = [i.split("|")[1] for i in df["sseqid"]]
            sp_accessions += list(set(df["sseqid"].to_list()))

    nr_accessions = list(set(sp_accessions))
    print('Total number of accessions:', len(nr_accessions))
    hit_lists = _chunks(nr_accessions, 300)


    sql = 'select uniprot_accession, gene,recommendedName_fullName,annotation_score,reviewed,t2.accession as taxon_id from uniprot_annotation t1' \
          ' inner join cross_references t2 on t1.uniprot_id=t2.uniprot_id ' \
          ' inner join crossref_databases t3 on t2.db_id=t3.db_id ' \
          ' where t3.db_name="NCBI Taxonomy" and t1.uniprot_accession in ("%s");'
    accession2annotation={}
    acc2score = {}
    for n,chunk in enumerate(hit_lists):
        print(f"{n} / {len(hit_lists)}")
        acc_filter = '","'.join(chunk)
        acc2score_tmp = get_swissprot_annotation(chunk)
        acc2score.update(acc2score_tmp)
        data = pandas.read_sql(sql % acc_filter, sqlite3_conn).set_index(["uniprot_accession"]).to_dict(orient="index")
        accession2annotation.update(data)
    

    print ("getting locus2protein_length")
    sql = 'select locus_tag,char_length(translation) from orthology_detail'
    cursor.execute(sql,)
    locus_tag2protein_length = manipulate_biosqldb.to_dict(cursor.fetchall())

    print ('loading blast results into database...')

    sql_template = 'insert into blastnr_blast_swissprot '
    sql_template += 'values (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s);'

    deleted_accession2annotation = {}

    for n_df, df in enumerate(df_list):
        print(f"Table {n_df}")
        n_file +=1
        row_index = 0
        for n, row in df.iterrows():
            sp_accession = row.sseqid
            # taxon
            # annotation score
            # protein names
            # genes
            # organism

            query_accession = row.qseqid
            try:
                annot = accession2annotation[sp_accession]
                subject_taxon_id = annot["taxon_id"]
                subject_description = annot["recommendedName_fullName"]
                subject_genes = annot["gene"]
                
            except KeyError:
                print("Deleted entry!----- %s" % sp_accession)
                if sp_accession not in deleted_accession2annotation:
                    # get a new accession if possible, download annotation and store it
                    new_accession = deleted_uniprot2new_identical_sequence(sp_accession)
                    time.sleep(1)
                    if new_accession:
                        temp_accession2annotation = get_swissprot_annotation([new_accession])
                        # add new annotation to the dictionnary
                        deleted_accession2annotation[sp_accession] = temp_accession2annotation[new_accession]
                    else:
                        # if no match, false
                        deleted_accession2annotation[sp_accession] = False
                annot = deleted_accession2annotation[sp_accession]
                if annot:
                    subject_taxon_id = annot[0]
                    subject_annot_score = annot[1].split(' ')[0]
                    subject_description = annot[2]
                    subject_genes = annot[3]
                    subject_organism = annot[4]
                else:
                    # set to unknown
                    subject_taxon_id = '32644'
                    subject_annot_score = '0'
                    subject_description = '-'
                    subject_genes = '-'
            try:
                subject_annot_score = acc2score[sp_accession]
            except:
                subject_annot_score = 0
            
            evalue = row.evalue
            percent_identity = float(row.pident)
            gaps = int(row.gapopen)
            length = int(row.length)
            query_start = int(row.qstart)
            query_end = int(row.qend)

            subject_start = int(row.sstart)
            subject_end = int(row.send)
            bit_score = float(row.bitscore)

            # first hit of the file
            if n == 0:
                hit_n = 1
            # check if new query accession)
            
            elif query_accession != df.iloc[row_index-1, :].qseqid:
                # new query => insert its first hit into blastnrdb
                hit_n = 1
            # same query, new hit
            else:
                hit_n+=1
            row_index+=1
            '''
            1 query_taxon_id int
            2 query_bioentry_id INT
            3 seqfeature_id INT
            4 hit_number int
            5 subject_accession varchar(200)
            6 subject_kingdom varchar(200)
            7 subject_scientific_name TEXT(2000000)
            8 subject_taxid INT
            9 subject_title VARCHAR(2000)
            10 evalue varchar(200)
            11 bit_score float
            12 percent_identity float
            13 gaps int
            14 length int
            15 query_start int
            16 query_end int
            17 query_cov float
            18 subject_start int
            19 subject_end
            20 genes
            21 annot score
            '''
            for locus_tag in hash2locus_list[query_accession]:
                seqfeature_id = locus_tag2seqfeature_id[locus_tag]
                query_cov = round(((query_end-query_start)/float(locus_tag2protein_length[locus_tag]))*100,2)

                values = [locus_tag2taxon_id[locus_tag],
                            locus_tag2bioentry_id[locus_tag],
                            seqfeature_id,
                            hit_n,
                            sp_accession,
                            subject_taxon_id,
                            subject_description,
                            str(evalue),
                            bit_score,
                            percent_identity,
                            gaps,
                            length,
                            query_start,
                            query_end,
                            query_cov,
                            subject_start,
                            subject_end,
                            subject_genes,
                            subject_annot_score]
                #try:

                cursor.execute(sql_template, values)
                #except:
                #    print ('problem with sql query')
                #    print (sql_template)
                #    print(values)
        # commit entire file
        conn.commit()


def create_sql_blast_swissprot_tables(db_name):
    
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)
    conn = server.adaptor.conn
    cursor = server.adaptor.cursor


    '''
    1 query_taxon_id int
    2 query_bioentry_id INT
    3 seqfeature_id INT
    4 hit_number int
    5 subject_accession varchar(200)
    6 subject_kingdom varchar(200)
    7 subject_scientific_name TEXT(2000000)
    8 subject_taxid INT
    9 subject_title VARCHAR(2000)
    10 evalue varchar(200)
    11 bit_score float
    12 percent_identity float
    13 gaps int
    14 length int
    15 query_start int
    16 query_end int
    17 query_cov float
    18 subject_start int
    19 subject_end
    20 genes
    21 annot score
    '''

    sql_plast = 'CREATE TABLE IF NOT EXISTS blastnr_blast_swissprot (query_taxon_id INT,' \
                ' query_bioentry_id INT,' \
                ' seqfeature_id INT,' \
                ' hit_number int,' \
                ' subject_accession varchar(200),' \
                ' subject_taxid INT,' \
                ' subject_title VARCHAR(2000),' \
                ' evalue varchar(200),' \
                ' bit_score float,' \
                ' percent_identity float,' \
                ' gaps int,' \
                ' length int,' \
                ' query_start int,' \
                ' query_end int,' \
                ' query_cov float,' \
                ' subject_start int,' \
                ' subject_end int,' \
                ' genes TEXT,' \
                ' annot_score int,' \
                ' INDEX query_taxon_id (query_taxon_id),' \
                ' INDEX query_bioentry_id (query_bioentry_id),' \
                ' INDEX hit_number (hit_number),' \
                ' INDEX seqfeature_id (seqfeature_id),' \
                ' INDEX subject_taxid(subject_taxid))' 

    try:

        cursor.execute(sql_plast)
        print ('sql hits ok')
        conn.commit()
    except:
        print (sql_plast)
        print ('not created')


def blastswiss2biosql( locus_tag2seqfeature_id,
                    db_name,
                    swissprot_sqlite,
                    n_procs,
                    hash2locus_list,
                    *input_blast_files):

    import numpy
    from multiprocessing import Process
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(db_name)
    conn = server.adaptor.conn
    cursor = server.adaptor.cursor


    create_sql_blast_swissprot_tables(db_name)


    print ('get locus2taxon_id')
    sql = 'select locus_tag, taxon_id from orthology_detail'
    locus_tag2taxon_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
    print ('get locus2bioentry')
    sql2 = 'select locus_tag,bioentry_id from biodatabase t1 ' \
           ' inner join bioentry as t2 on t1.biodatabase_id=t2.biodatabase_id' \
           ' inner join orthology_detail t3 on t2.accession=t3.accession where t1.name="%s"' % (db_name)

    locus_tag2bioentry_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2,))

    load_blastswissprot_file_into_db(locus_tag2taxon_id,
                                     locus_tag2seqfeature_id,
                                     locus_tag2bioentry_id,
                                     input_blast_files,
                                     biodb,
                                     swissprot_sqlite,
                                     hash2locus_list)

    sys.stdout.write("done!")





if __name__ == '__main__':
    import argparse
    from chlamdb.biosqldb import manipulate_biosqldb
    import sys
    import os
    import chlamdb_setup_utils

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_blast', type=str, help="input blast tab files", nargs='+')
    parser.add_argument("-d", '--mysql_database', type=str, help="Biosql biodatabase name")
    parser.add_argument("-s", '--swissprot_sqlite', type=str, help="Swissprot sqlite database")
    parser.add_argument("-t", '--create_tables', action='store_true', help="Create SQL tables")
    parser.add_argument("-p", '--n_procs', type=int, help="Number of threads to use (default=8)", default=8)
    parser.add_argument("-c", '--clean_tables', action='store_true', help="delete all sql tables")
    parser.add_argument("-l", '--load_tables', action='store_true', help="load tab files into biodatabase")
    parser.add_argument("-f", '--filter_n_hits', action='store_true', help="filter_n_hits (max 100 hits/locus)")
    parser.add_argument("-u", '--hash2locus_tag', type=str, help="Tab separated file with correspondance between sequence hashes and locus tags")

    args = parser.parse_args()


    mysql_host = 'localhost'
    mysql_user = 'root'

    mysql_pwd = os.environ['SQLPSW']
    mysql_db = args.mysql_database

    biodb = args.mysql_database

    hash2locus_list = chlamdb_setup_utils.get_hash2locus_list(args.hash2locus_tag)

    if args.load_tables:

        server, db = manipulate_biosqldb.load_db(biodb)

        sys.stdout.write("creating locus_tag2seqfeature_id")
        locus_tag2seqfeature_id = manipulate_biosqldb.locus_tag2seqfeature_id_dict(server, biodb)

        sys.stdout.write("creating protein_id2seqfeature_id")
        protein_id2seqfeature_id = manipulate_biosqldb.protein_id2seqfeature_id_dict(server, biodb)

        blastswiss2biosql(locus_tag2seqfeature_id,
                        biodb,
                        args.swissprot_sqlite,
                        args.n_procs,
                        hash2locus_list,
                        *args.input_blast)

    manipulate_biosqldb.update_config_table(biodb, "BLAST_swissprot")