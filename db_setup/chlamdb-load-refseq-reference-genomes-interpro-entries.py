#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-



def get_all_entries():
    import requests
    import pandas
    import time
    db_list = {'PANTHER', 'PFAM', 'CDD', 'PROSITE', 'SMART', 'HAMAP', 'SSF', 'SFLD', 'PRINTS', 'PROFILE', 'PIRSF', 'NCBIfam', 'CATHGENE3D'}
    #db_list = ["HAMAP"]
    
    for db in db_list:
        print(db)
        link = f"https://www.ebi.ac.uk/interpro/api/entry/{db}/?page_size=30000"
        r = requests.get(url = link)
        j = r.json()
        data = []
        for entry in j["results"]:
            acc = entry["metadata"]["accession"]
            name = entry["metadata"]["name"]
            type = entry["metadata"]["type"]
            integrated = entry["metadata"]["integrated"]
            data.append([db, acc, name, type, integrated])
        df = pandas.DataFrame(data, columns = ["db","accession", "name", "type", "integrated"])
        df.to_csv(f"{db}_entries.tsv", sep="\t")
        time.sleep(3)
        

def insert_interpro_entries(biodb, tables):
    import pandas 
    import os 
    sqlpsw = os.environ['SQLPSW']
    
    from sqlalchemy import create_engine
    
    # engine = create_engine("postgresql+psycopg2://scott:tiger@localhost:5432/mydatabase")
    engine = create_engine(f"mysql://root:{sqlpsw}@127.0.0.1/{biodb}")
    engine_conn = engine.connect()
    conn_mysql = engine.raw_connection()
    cursor = conn_mysql.cursor()
    
    sql = 'create table refseq_ref_repres_genomes_interpro_entries ( ' \
          ' signature_id INT AUTO_INCREMENT PRIMARY KEY, '  \
          ' analysis_id INT, '  \
          ' accession varchar(200),' \
          ' name TEXT,' \
          ' type varchar(200),' \
          ' interpro_id INT,' \
          ' index signature_id(signature_id)) '
          
    cursor.execute(sql,)
    conn_mysql.commit()
    
    df_interpro_entries = pandas.read_sql_query('select name,interpro_id from interpro_entry', engine, index_col="name")
    df_db_entries = pandas.read_sql_query('select analysis_name,analysis_id from interpro_analysis', engine, index_col="analysis_name")
    
    for table in tables:
        df = pandas.read_csv(table, sep="\t", index_col=0)
        df["interpro_id"] = [df_interpro_entries.loc[interpro_accession, "interpro_id"] if not pandas.isna(interpro_accession) else None for interpro_accession in df["integrated"] ]
        df["analysis_id"] = [df_db_entries.loc[db, "analysis_id"] for db in df["db"] ]
        print(df.head())
        df = df.drop("integrated", axis=1)
        df = df.drop("db", axis=1)
        df.to_sql("refseq_ref_repres_genomes_interpro_entries", engine, if_exists="append", index=False)

def insert_databases(biodb, db_list):
    import pandas 
    import os 
    sqlpsw = os.environ['SQLPSW']
    
    from sqlalchemy import create_engine
    
    # engine = create_engine("postgresql+psycopg2://scott:tiger@localhost:5432/mydatabase")
    engine = create_engine(f"mysql://root:{sqlpsw}@127.0.0.1/{biodb}")
    engine_conn = engine.connect()
    conn_mysql = engine.raw_connection()
    cursor = conn_mysql.cursor()
    
    cursor.execute('select * from interpro_analysis;')
    db_in_table = [i[1] for i in cursor.fetchall()]
    print(db_in_table)

    for db in db_list:
        if db not in  db_in_table:
            print("db not in table", db)
            sql = 'insert into interpro_analysis(analysis_name) values (%s)'
            cursor.execute(sql, [db])
    conn_mysql.commit()


def update_interpro_entry_list(biodb):
    import pandas 
    import os 
    sqlpsw = os.environ['SQLPSW']
    
    from sqlalchemy import create_engine
    
    # engine = create_engine("postgresql+psycopg2://scott:tiger@localhost:5432/mydatabase")
    engine = create_engine(f"mysql://root:{sqlpsw}@127.0.0.1/{biodb}")
    engine_conn = engine.connect()
    conn_mysql = engine.raw_connection()
    cursor = conn_mysql.cursor()
    import urllib.request
    link = 'https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/entry.list'
    entry_list = urllib.request.urlopen(link).read().decode('utf-8').split("\n")
    cursor.execute('select name from interpro_entry')
    entry_list_in_db = set([i[0] for i in cursor.fetchall()])
    print(entry_list_in_db)
    count = 0
    for i, line in enumerate(entry_list):
        
        if not 'IPR' in line:
            continue
        
        data = line.rstrip().split("\t")
        accession = data[0]
        description = data[2]
        if accession not in entry_list_in_db:
            sql = 'insert into interpro_entry (name, description) values (%s, %s)'
            cursor.execute(sql,[accession, description])
            count+=1
        conn_mysql.commit()
    print(count)



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", '--mysql_database', type=str, help="Biosql biodatabase name")
    parser.add_argument("-t", '--table',type=str,help="RefSeq table")
    parser.add_argument("-e", '--interpro_entries',type=str,help="DB entries", nargs="+")
    
    
    args = parser.parse_args()
    
    #load_reference_genome_table_interpro_precomp(args.assembly_accession, args.mysql_database, args.table)
    #get_all_entries()
    #insert_databases( args.mysql_database, {'PANTHER', 'PFAM', 'CDD', 'PROSITE', 'SMART', 'HAMAP', 'SSF', 'SFLD', 'PRINTS', 'PROFILE', 'PIRSF', 'NCBIfam', 'CATHGENE3D'})
    insert_interpro_entries(args.mysql_database, args.interpro_entries)
    #update_interpro_entry_list(args.mysql_database)