#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

def load_reference_genome_table_into_database(genome_refseq_file, biodb):

    import os
    import pandas
    sqlpsw = os.environ['SQLPSW']
    
    from sqlalchemy import create_engine
    
    # engine = create_engine("postgresql+psycopg2://scott:tiger@localhost:5432/mydatabase")
    engine = create_engine(f"mysql://root:{sqlpsw}@127.0.0.1/{biodb}")
    engine_conn = engine.connect()
    conn_mysql = engine.raw_connection()
    
    cursor = conn_mysql.cursor()

    # ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt

    sql = 'create table refseq_ref_repres_genomes (assembly_id INT AUTO_INCREMENT PRIMARY KEY,' \
          ' assembly_accession varchar(200),' \
          ' bioproject varchar(200),' \
          ' biosample varchar(200),' \
          ' wgs_master varchar(200),' \
          ' refseq_category varchar(200),' \
          ' taxid INT,' \
          ' species_taxid INT,' \
          ' organism_name varchar(200),' \
          ' infraspecific_name varchar(200),' \
          ' version_status varchar(200),' \
          ' assembly_level varchar(200),' \
          ' release_type varchar(200),' \
          ' genome_rep varchar(200),' \
          ' seq_rel_date varchar(200),' \
          ' ftp_path TEXT,' \
          ' n_mapped_uniparc INT,' \
          ' n_non_mapped_uniparc INT,' \
          ' keep_interpro BOOLEAN,' \
          ' index species_taxid(species_taxid))'

    print("create table")
    cursor.execute(sql,)
    conn_mysql.commit()
    print("ok")

    '''
    assembly_accession
    bioproject
    biosample
    wgs_master	
    refseq_category	
    taxid	
    species_taxid	
    organism_name	
    infraspecific_name	
    isolate	
    version_status	
    assembly_level	
    release_type
    genome_rep	
    seq_rel_date	
    asm_name	
    submitter	
    gbrs_paired_asm	
    paired_asm_comp	
    ftp_path	
    excluded_from_refseq	
    relation_to_type_material	
    asm_not_live_date
    
    '''

    cols_keep = [
        "assembly_accession",
        "bioproject",
        "biosample",
        "wgs_master",
        "refseq_category",
        "taxid",
        "species_taxid",
        "organism_name",
        "infraspecific_name",
        "version_status",
        "assembly_level",
        "release_type",
        "genome_rep",
        "seq_rel_date",
        "ftp_path",
    ]

    genome_refseq_table = pandas.read_csv(genome_refseq_file, sep="\t")[cols_keep]
    genome_refseq_table_filt = genome_refseq_table[genome_refseq_table["refseq_category"].isin(["reference genome", "representative genome"])]
    
    print("insert data")
    genome_refseq_table_filt.to_sql("refseq_ref_repres_genomes", 
                                    engine, 
                                    index=False, 
                                    if_exists="append")
    


def add_uniparc_mapping_stats(biodb, uniparc_match_stat_table):
    import os
    import pandas
    sqlpsw = os.environ['SQLPSW']
    
    from sqlalchemy import create_engine
    
    # engine = create_engine("postgresql+psycopg2://scott:tiger@localhost:5432/mydatabase")
    engine = create_engine(f"mysql://root:{sqlpsw}@127.0.0.1/{biodb}")
    engine_conn = engine.connect()
    conn_mysql = engine.raw_connection()   
    cursor = conn_mysql.cursor()
    
    sql = 'select assembly_accession, assembly_id from refseq_ref_repres_genomes'
    df_assembly = pandas.read_sql_query(sql, engine, index_col="assembly_accession")
    
    df = pandas.read_csv(uniparc_match_stat_table, sep="\t")
    
    for n,row in df.iterrows():
        if n % 100 == 0:
            print(n)
        asseembly_id = df_assembly.loc[row.genome, "assembly_id"]
        sql = f'update refseq_ref_repres_genomes set n_mapped_uniparc={row.n_mapped}, n_non_mapped_uniparc={row.n_non_mapped}, keep_interpro={row.keep} where assembly_id={asseembly_id}'
        cursor.execute(sql)
        #print(sql)
    conn_mysql.commit()

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", '--mysql_database', type=str, help="Biosql biodatabase name")
    parser.add_argument("-t", '--table',type=str,help="RefSeq table")
    parser.add_argument("-u", '--uniparc_mapping_stats',type=str,help="Unparc mapping table (interpro_stats.tsv)")
    
    args = parser.parse_args()
    
    add_uniparc_mapping_stats(args.mysql_database, args.uniparc_mapping_stats)
    
    #load_reference_genome_table_into_database(args.table, args.mysql_database)
    