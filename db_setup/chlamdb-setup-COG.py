#!/usr/bin/env python


def create_COG_tables(biodb):
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)
    conn = server.adaptor.conn
    cursor = server.adaptor.cursor


    sql2 = 'CREATE table if not exists COG_cog_names_2014 (COG_id INTEGER PRIMARY KEY AUTO_INCREMENT,' \
                                            ' COG_name varchar(100),' \
                                            ' description varchar(200))'

    sql3 = 'CREATE table if not exists COG_cog_id2cog_category (COG_id INT,' \
           ' category_id int)'


    sql4 = 'create table if not exists COG_code2category (code varchar(5),' \
                                           ' description TEXT,' \
                                           ' category_id INTEGER primary KEY AUTO_INCREMENT);'
                                           
    sql_index2 = 'create index cciccid on COG_cog_id2cog_category(COG_id)'
    sql_index3 = 'create index ccicccid on COG_cog_id2cog_category(category_id)'


    cursor.execute(sql2,)
    cursor.execute(sql3,)
    cursor.execute(sql4,)
    try:
        cursor.execute(sql_index2,)
        cursor.execute(sql_index3,)
    except:
        pass
    
    conn.commit()

def COG_id2cog_data_from_website(GOG_id):
    import urllib2
    from bs4 import BeautifulSoup
    module_file_file = 'https://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid=%s' % GOG_id
    data = urllib2.urlopen(module_file_file)
    # COG0001	H	Glutamate-1-semialdehyde aminotransferase
    soup = BeautifulSoup(data, "lxml")
    html = soup.encode('utf-8')
    title = soup.findAll("div", { "class" : "desctit" })[0]
    description = title.text.split('[')[0]
    category = title.text.split('[')[1].split(']')[0]
    return [description, category]


def load_cog_tables(biodb,
                    cognames_2014, 
                    cog_categories):
    '''

    COG names
    COG0001	H	Glutamate-1-semialdehyde aminotransferase

    COG
    158333741,Acaryochloris_marina_MBIC11017_uid58167,158333741,432,1,432,COG0001,0,
    <domain-id>, <genome-name>, <protein-id>,<protein-length>,
    <domain-start>, <domain-end>, <COG-id>, <membership-class>,

    :param cognames_2014:
    :param cog_2014:
    :return:
    '''
    
    import MySQLdb
    import os
    import re
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)
    conn = server.adaptor.conn
    cursor = server.adaptor.cursor
    
    
    try:
        conn.set_character_set('utf8')
        cursor.execute('SET NAMES utf8;')
        cursor.execute('SET CHARACTER SET utf8;')
        cursor.execute('SET character_set_connection=utf8;')
    except:
        cursor.execute('PRAGMA encoding="UTF-8";')

    sql = 'insert into COG_code2category (code, description) values (%s, %s)'
    for n,row in cog_categories.iterrows():
        cursor.execute(sql, [row.code, row.description])
    conn.commit()

    sql2 = 'select code, category_id from COG_code2category'
    cursor.execute(sql2,)
    cog_category2cog_category_id = manipulate_biosqldb.to_dict(cursor.fetchall())

    sql3 = 'select description, category_id from COG_code2category'
    cursor.execute(sql3,)
    cog_description2cog_category_id = manipulate_biosqldb.to_dict(cursor.fetchall())

    for n, row in cognames_2014.iterrows():
        # split multiple function and insert one row/function
        fonctions = list(row.functional_class)

        sql = 'insert into COG_cog_names_2014(COG_name, description) values(%s, %s)'
        cursor.execute(sql, [row.COG_id, 
                             re.sub('"','',row.COG_annotation)])
        
        cog_id = cursor.lastrowid
        
        sql = 'insert into COG_cog_id2cog_category values (%s, %s)'
        for fonction in fonctions:
            cursor.execute(sql,[cog_id,
                                cog_category2cog_category_id[fonction]])

    conn.commit()

if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    import urllib.request
    import pandas 
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--cognames_2014', type=str, help="cognames2003-2014.tab")
    parser.add_argument("-f", '--functions_2014', type=str, help="fun2003-2014.tab")
    parser.add_argument("-b", '--biodb', type=str, help="Biodb name")
    parser.add_argument("-d", '--download', action='store_true', help="Download from NCBI FTP")
    

    args = parser.parse_args()
    
    create_COG_tables(args.biodb)
    
    cognames = 'ftp://ftp.ncbi.nlm.nih.gov/pub/COG/COG2014/data/cognames2003-2014.tab'
    cog_functions = 'ftp://ftp.ncbi.nlm.nih.gov/pub/COG/COG2014/data/fun2003-2014.tab'
    
    if args.download:
        print("Downloading https://ftp.ncbi.nlm.nih.gov/pub/COG/COG2014/data/fun2003-2014.tab")
        cog_functions_data = urllib.request.urlopen(cog_functions).read().decode('utf-8').split('\n')

        print("Downloading ftp://ftp.ncbi.nlm.nih.gov/pub/COG/COG2014/data/cognames2003-2014.tab...")
        cognames_data = urllib.request.urlopen(cognames).read().decode('windows-1252').split('\n')
        
        print("done")
        print(cog_functions_data)
        load_cog_tables(args.biodb,
                        cognames_data,
                        cog_functions_data)

    else:

        cognames_2014 = pandas.read_csv(args.cognames_2014, names=["COG_id", "functional_class", "COG_annotation"], sep="\t", comment="#", encoding = "ISO-8859-1")

        functions_2014 = pandas.read_csv(args.functions_2014, names=["code", "description"], sep="\t", comment="#", encoding = "ISO-8859-1")
                  
        load_cog_tables(args.biodb,
                        cognames_2014, 
                        functions_2014)
