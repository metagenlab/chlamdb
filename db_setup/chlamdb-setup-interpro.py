#!/usr/bin/env python

def get_interpro_entry_tables(interpro_release= '60.0'):
    '''

    :param interpro_release: par exemple: 60.0
    :return:
    '''
    import urllib2
    import MySQLdb
    import os
    sqlpsw = os.environ['SQLPSW']

    conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                user="root", # your username
                                passwd=sqlpsw, # your password
                                db="interpro") # name of the data base
    cursor = conn.cursor()



    sql = 'CREATE table entry (interpro_id INT AUTO_INCREMENT, name varchar(400), description TEXT, index interpro_id(interpro_id))'

    cursor.execute(sql,)
    conn.commit()
    link = 'ftp://ftp.ebi.ac.uk/pub/databases/interpro/%s/entry.list' % interpro_release
    print (link)
    req = urllib2.Request(link)
    entry_list = urllib2.urlopen(req)
    for line in entry_list:
        if not 'IPR' in line:
            continue
        accession = line.rstrip().split(' ')[0]
        description = ' '.join(line.rstrip().split(' ')[1:])
        sql = 'insert into entry (name, description) values ("%s", "%s")' % (accession, description)
        cursor.execute(sql,)
    conn.commit()
    

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--version', type=str, help='interpro_version (e.g 60.0)')


    args = parser.parse_args()

    get_interpro_entry_tables(args.version)
