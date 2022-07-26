#!/usr/bin/env python

from chlamdb.eutils import tcdb_utils
import MySQLdb
import re
import os
from chlamdb.eutils import accession2taxon_id
from chlamdb.eutils import sequence_id2scientific_classification
from Bio import SeqIO
from bs4 import BeautifulSoup
from io import StringIO
import urllib.request
from chlamdb.biosqldb import blastswiss2sqltable
import logging
def get_protein_description(acc="1.A.1.9"):
    import urllib.request
    from bs4 import UnicodeDammit
    
    url= 'http://www.tcdb.org/search/result.php?tc=%s' % acc
    
    page = urllib.request.urlopen(url).read().decode('Latin1')
    
    soup = BeautifulSoup(page, 'lxml')
    
    #table = soup.find("table", {"id": "result-cluster"})

    superfam2description = {}
    
    # retrieve fam list and description
    famlist = soup.findAll("td", {"id": "right-border"})
    for fam in famlist:
        txt = fam.text.strip()
        try:
            tcid = re.findall("[0-9]+\.[A-Z]\.[0-9]+\.[0-9]+", txt)[0]
            description = txt.split(":")[1].strip()
            superfam2description[tcid] = description
        except:
            pass
        
    rows = soup.findAll('tr')
    
    tcid2data = {}
    #print("Number of transporters:", len(rows))
    for n, row in enumerate(rows):
        cols = row.find_all('td')
        # BeautifulSoup(raw_html, "lxml").text
        cols = [BeautifulSoup(ele.text.strip(), "lxml").text.strip() for ele in cols]
        #print (len(cols))
        #print(cols)
        if len(cols) > 1:
            tcid = re.sub("\*", "", cols[0])
            tcid2data[tcid] = {}
            try:
                tcid2data[tcid]["name"] = cols[1].replace("\r\n", "").replace("\xa0", "")
                tcid2data[tcid]["domain"] = cols[2]
                tcid2data[tcid]["superkingdom"] = cols[3]
            except:
                logging.warning("problem with row", cols)
                continue
    return tcid2data, superfam2description

def get_tcid2substrate():
    import urllib.request
    
    url = 'https://tcdb.org/cgi-bin/substrates/getSubstrates.py'

    data = urllib.request.urlopen(url).read().decode('utf-8').split('\n')
    
    tcid2substrate = {}
    
    for row in data:
        if len(row) == 0:
            continue
        row_data = [i for i in row.split("\t") if len(i) != 0]
        if len(row_data) != 2:
            logging.warning("incomplete row: %s" % row)
        else:
            tcid, substrate = row_data
        
        tcid2substrate[tcid.strip()] = substrate.strip()
    
    return tcid2substrate

tcid2substrate = get_tcid2substrate()
print(tcid2substrate)