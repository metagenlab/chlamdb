#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

def load_reference_genome_table_interpro_precomp(biodb, table_list, assembly_accession=False, genome_filter=False):

    import os
    import pandas
    sqlpsw = os.environ['SQLPSW']
    
    from sqlalchemy import create_engine
    
    # engine = create_engine("postgresql+psycopg2://scott:tiger@localhost:5432/mydatabase")
    engine = create_engine(f"mysql://root:{sqlpsw}@127.0.0.1/{biodb}")
    engine_conn = engine.connect()
    conn_mysql = engine.raw_connection()   
    cursor = conn_mysql.cursor()
    
    if not assembly_accession:
        acc = False

    sql = 'create table refseq_ref_repres_genomes_interpro_annots ( ' \
          ' assembly_id INT,'  \
          ' signature_id INT,'  \
          ' accession varchar(200),' \
          ' start INT,'  \
          ' end INT,'  \
          ' score float)' 
    
    try:
        cursor.execute(sql,)
        conn_mysql.commit()
    except:
        pass
    
    if genome_filter:
        filtering = pandas.read_csv(genome_filter, sep="\t", index_col="genome")
    
    df_signatures = pandas.read_sql_query('select accession as entry_accession,signature_id from refseq_ref_repres_genomes_interpro_entries', engine, index_col="entry_accession")

    df_assembly = pandas.read_sql_query('select assembly_accession,assembly_id from refseq_ref_repres_genomes', engine, index_col="assembly_accession")

    df_list = []
    for n,table in enumerate(table_list):
        #print(table)
        if n % 50 == 0 and n != 0:
            
            print(f"{n}\t{len(table_list)} INSERT")
            #df_merge = pandas.concat(df_list)
            #df_merge.to_sql("refseq_ref_repres_genomes_interpro_annots", 
            #                    engine, 
            #                    index=False, 
            #                    if_exists="append", chunksize=50000)
            df_list = []
        
        if not acc:
            assembly_accession = os.path.basename(table).split("_uniparc_mapping")[0]
        if genome_filter:
            if filtering.loc[assembly_accession, "keep"] == False:
                print(f"{assembly_accession}\tskip")
                continue
        
        
        assembly_id = df_assembly.loc[assembly_accession, "assembly_id"]
        
        df_annots = pandas.read_csv(table, sep="\t")
        
        df_annots["assembly_id"] = assembly_id
        
        df_annots["entry_accession"] = [signature_accession if not signature_accession.startswith("PTHR") else signature_accession.split(":")[0] for signature_accession in df_annots["entry_accession"]]
        
        signature_list = set(df_annots["entry_accession"].to_list())

        #db_list = set(df_annots["db_name"].to_list())
        
        #print("count missing")
        #missing = [i for i in signature_list if i not in df_signatures.index]
        #print(f"{assembly_accession}\t{len(signature_list)}\t{len(missing)}\t{missing[0:10]}", )
        print(f"{assembly_accession}\t{len(signature_list)}")
        # remove eventual missing entries
        #df_annots = df_annots[~df_annots["entry_accession"].isin(missing)]
        
        #print("DB list", len(db_list), db_list)
        
        df_annots = df_annots.set_index("entry_accession").join(df_signatures).reset_index(drop=True)
                
        df_annots = df_annots.drop("interpro_accession", axis=1)
        df_annots = df_annots.drop("evidence_name", axis=1)
        df_annots = df_annots.drop("db_name", axis=1)
        df_annots = df_annots.drop("name", axis=1)
        
        df_annots.to_sql("refseq_ref_repres_genomes_interpro_annots", 
                                        engine, 
                                        index=False, 
                                        if_exists="append")
                

    
    

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", '--mysql_database', type=str, help="Biosql biodatabase name")
    parser.add_argument("-t", '--table',type=str,help="interpro precomp table", nargs="+")
    parser.add_argument("-a", '--assembly_accession',type=str,help="RefSeq table")
    parser.add_argument("-s", '--genome_filter',type=str,help="interpro_stats.tsv")
    
    args = parser.parse_args()
    
    load_reference_genome_table_interpro_precomp(args.mysql_database, args.table, args.assembly_accession, args.genome_filter)