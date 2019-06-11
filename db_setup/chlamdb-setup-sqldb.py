#!/usr/bin/env python

import urllib.request
import sys
import MySQLdb
import os

url_mysql_biosql_scheme = 'https://raw.githubusercontent.com/biosql/biosql/master/sql/biosqldb-mysql.sql'

sys.stdout.write('Downloading Biosql scheme from %s ...\n' % url_mysql_biosql_scheme)
request = urllib.request.Request(url_mysql_biosql_scheme)
page = urllib.request.urlopen(request)

sys.stdout.write("Creating mysql database...\n")

sqlpsw = os.environ['SQLPSW']
conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                     user="root", # your username
                     passwd=sqlpsw) # name of the data base
cursor = conn.cursor()
sql_db = 'CREATE DATABASE IF NOT EXISTS biosqldb;
cursor.execute(sql_db,)
conn.commit()
cursor.execute("use biosqldb;",)

sys.stdout.write("Importing Biosql schema...\n")
conn.executescript(page.read().decode('unicode-escape'))

conn.commit()
