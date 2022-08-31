#!/usr/bin/env python

from chlamdb.biosqldb import manipulate_biosqldb
import sqlite3 
from Bio import SeqIO
from Bio import Medline
from Bio import Entrez
Entrez.email = 'trestan.pillone@chuv.ch'
import pandas
from sqlalchemy import create_engine
from pandas.core.series import Series
class StringPMID():
    def __init__(self,
                 string_sqlite,
                 blast_results,
                 db_name,
                 hash2locus_df,
                 query_fasta_file,
                 string_sqlite_fa,
                 n_best_hits=20):
        import pandas 
        import os
        sqlpsw = os.environ['SQLPSW']
    
        engine = create_engine(f"mysql://root:{sqlpsw}@127.0.0.1/{db_name}")
        self.conn_mysql = engine.connect()
        
        self.string_conn = sqlite3.connect(string_sqlite)
        self.string_cursor = self.string_conn.cursor()
        
        self.string_conn_fa = sqlite3.connect(string_sqlite_fa)
        self.string_cursor_fa = self.string_conn_fa.cursor()

        self.server, self.db = manipulate_biosqldb.load_db(db_name)
        self.blast_results = blast_results
        self.hash2locus_df = hash2locus_df
        self.db_name = db_name
        logging.info("parse fasta")
        self.query2len = self.parse_fasta(query_fasta_file)
        logging.info("done")
        self.missing_pmid_data =[] 

        self.n_best_hits = n_best_hits
        # create tables if not exists
        self.create_tables()

        logging.info("retrieve string data from db")
        # retrieve data from string db
        sql = 'select publication_id,publication_date,publication_source,linkout_url,authors,title from publications;'       
        # pmid, title, journal, year, linkout_url, authors
        self.article_data_df = pandas.read_sql(sql, self.string_conn) 
        self.article_data_df = self.article_data_df.rename(columns={"publication_date":"year","publication_source":"journal"})
        self.article_data_df["pmid"] = self.article_data_df.apply (lambda row: int(row.publication_id.split(":")[1]), axis=1)
        self.article_data_df = self.article_data_df.drop("publication_id", axis="columns").set_index("pmid")
        
        logging.info(self.article_data_df.head())
        sql_species = 'select species_id,compact_name as organism from species'
        self.species_id2species_name_df = pandas.read_sql(sql_species, self.string_conn).set_index("species_id")
        logging.info(self.species_id2species_name_df.head())     
        # retrieve protein hash from biosqldb
        sql = f'select distinct hash from string_seqfeature_id2string_protein_mapping t1 inner join annotation_hash2seqfeature_id t2 on t1.seqfeature_id=t2.seqfeature_id;'
        self.hash_in_db =set([i[0] for i in self.server.adaptor.execute_and_fetchall(sql,)])
        
        # check if some data were already inserted
        self.get_pmid_in_db()
        logging.info(f"{len(self.pmid_in_db)} PMID already in db")
        self.get_string_proteins2string_protein_id()


        self.blast_header = [ "qseqid",
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
        logging.info("done")

    def get_pmid_in_db(self,):
        sql = 'select pmid from string_pmid2data_stringdb'
        self.pmid_in_db = set([i[0] for i in self.server.adaptor.execute_and_fetchall(sql,)])
    def get_string_proteins2string_protein_id(self,):
        sql = 'select accession,id as string_protein_id from string_string_protein_entry'
        self.string_proteins2string_protein_id = pandas.read_sql(sql, self.server.adaptor.conn).set_index("accession")


    def create_tables(self):
        sql1 = 'create table if not exists string_seqfeature_id2string_protein_mapping (seqfeature_id INT, string_protein_id INT, query_cov FLOAT, hit_cov FLOAT, identity FLOAT, evalue FLOAT, score FLOAT)'
        sql2 = 'create table if not exists string_string_protein2pmid (string_protein_id INT, PMID BIGINT)'
        sql3 = 'create table if not exists string_pmid2data_stringdb (pmid BIGINT, title TEXT, journal TEXT, year varchar(200), linkout_url TEXT, authors TEXT)'
        sql4 = 'create table if not exists string_string_protein_entry (id INT AUTO_INCREMENT PRIMARY KEY, accession varchar(200), organism TEXT, description TEXT, preferred_name varchar(400))'
        self.server.adaptor.execute(sql1)
        self.server.adaptor.execute(sql2)
        self.server.adaptor.execute(sql3)
        self.server.adaptor.execute(sql4)

    def index_tables(self,):
        sql1 = 'CREATE index ssspm ON string_seqfeature_id2string_protein_mapping(seqfeature_id)'
        sql2 = 'CREATE index ssspm2 ON string_seqfeature_id2string_protein_mapping(string_protein_id)'
        sql3 = 'CREATE index sspm ON string_string_protein2pmid(string_protein_id)'
        sql4 = 'CREATE index spds ON string_pmid2data_stringdb(pmid)'
        self.server.adaptor.execute(sql1)
        self.server.adaptor.execute(sql2)
        self.server.adaptor.execute(sql3)
        self.server.adaptor.execute(sql4)

    def parse_fasta(self, fasta_file):
        accession2len ={}
        if fasta_file.endswith("gz"):
            import gzip
            fagz = gzip.open(fasta_file, 'rt')
            records = SeqIO.parse(fagz, "fasta")
        else:
            records = SeqIO.parse(fasta_file, "fasta")
        for record in records:
            accession2len[record.id] = len(record.seq) 
        return accession2len


    def insert_new_string_entry(self, protein_external_id, annotation, preferred_name, nr_pmid_list):
        #logging.info(protein_external_id, nr_pmid_list)
        species_id = protein_external_id.split(".")[0]
        species_name = self.species_id2species_name[species_id]
        
        sql = f'insert into string_string_protein_entry (accession, organism, description, preferred_name) values (%s, %s, %s, %s)'
        self.server.adaptor.execute(sql, (protein_external_id, species_name, annotation, preferred_name))
        
        # add new id to dictionnary
        string_protein_id = self.server.adaptor.last_id("string_string_protein_entry")
        self.string_proteins2string_protein_id[protein_external_id] = string_protein_id
        
        # insert data if not already present in db
        nr_pmid_list = [pmid for pmid in nr_pmid_list if pmid not in self.pmid_in_db]
        
        for pmid in nr_pmid_list:
            # PMID:11825770
            #pmid = pmid.split(":")[1]
            sql = f'insert into string_string_protein2pmid (string_protein_id, pmid) values ({string_protein_id}, {pmid})'
            self.server.adaptor.execute(sql,)
            article_data = self.pmid2article_data[pmid]
            publication_date, publication_source, linkout_url, authors, title = article_data
            #sql2 = f'insert into string_pmid2data_stringdb (pmid, title, journal, year, linkout_url, authors) values (%s, %s, %s, %s, %s, %s)'
            #self.server.adaptor.execute(sql2, (pmid, title, publication_source, publication_date, linkout_url, authors))
            #self.pmid_in_db.add(pmid)

    def load_string_data(self):
        
        db_name = self.db_name
                
        sql = f'select locus_tag, seqfeature_id from custom_tables_locus2seqfeature_id'
        locus_tag2seqfeature_id_df = pandas.read_sql(sql, self.conn_mysql).set_index("locus_tag")

        for blast_file in self.blast_results:
            logging.info("Reading blast file...")
            df = pandas.read_csv(blast_file, sep="\t", names=self.blast_header).set_index("qseqid")
            df_hits = df.groupby(["qseqid"]).head(self.n_best_hits)
            logging.info("exclude hash already in db")
            df_hits = df_hits.loc[df.index.difference(self.hash_in_db), :]
            if len(df_hits) == 0:
                logging.info("all qseqid already in db!")
                continue
            logging.info("create temp table with string protein ids to join on pmid data")            
            ids = df_hits[["sseqid"]].drop_duplicates()
            ids.to_sql("stringbbh", self.string_conn, index=False, if_exists="replace")
            ids.to_sql("stringbbh", self.string_conn_fa, index=False, if_exists="replace")
            self.string_conn.commit()
            self.string_conn_fa.commit()
            logging.info("done, joining")
            sql_protein_data = 'select t1.*,t2.protein_external_id, t2.annotation as description, t2.preferred_name, t3.publication_id from stringbbh t1 inner join proteins t2 on t1.sseqid=t2.protein_external_id inner join protein_id2pmid t3 on t2.protein_id=t3.protein_id' 
            # skip already inserted entries
            df_string = pandas.read_sql(sql_protein_data, self.string_conn).set_index("protein_external_id")
            df_string["pmid"] = df_string.apply (lambda row: int(row.publication_id.split(":")[1]), axis=1)
            nr_pmid_list = set(df_string["pmid"].to_list())
            pmid_not_in_db = nr_pmid_list.difference(self.pmid_in_db) 
            logging.info(f"Adding {len(pmid_not_in_db)} missing from db")
            df_missing = self.article_data_df.loc[pmid_not_in_db]
            df_missing.to_sql("string_pmid2data_stringdb", self.conn_mysql, if_exists="append")
            logging.info(f"Added {len(df_missing)} pmid missing from db")
            logging.info(f"Loading string protein entries")
            # species_id = protein_external_id.split(".")[0]
            # not sure if it works to compare set to 
            missing_string_prot_entries = df_hits.set_index("sseqid").loc[set(df_hits["sseqid"].to_list()).difference(set(self.string_proteins2string_protein_id.index))]
            logging.info(f"{len(missing_string_prot_entries)} / {len(df_hits)} missing strint prot entries")
            logging.info(f"Extracting species id")
            if len(missing_string_prot_entries) > 0:
                missing_string_prot_entries["species_id"] = missing_string_prot_entries.apply (lambda row: int(row.name.split(".")[0]), axis=1)
                # 243161.TC_0001 | Chlamydia muridarum   | Catalyzes an early ... | hemB
                logging.info(f"Joining df_string to get protein description and preferred_name")
                missing_string_prot_entries = missing_string_prot_entries.join(df_string)
                missing_string_prot_entries = missing_string_prot_entries.set_index("species_id").join(self.species_id2species_name_df).reset_index()[["sseqid", "organism", "description", "preferred_name"]].drop_duplicates()
                missing_string_prot_entries.columns = ["accession", "organism", "description", "preferred_name"]
                logging.info("Insert string_string_protein_entry")
                # string_string_protein_entry (`index`, sseqid, organism, description, preferred_name)
                missing_string_prot_entries.to_sql("string_string_protein_entry", self.conn_mysql, if_exists="append", index=False)
                self.server.adaptor.commit()  
                self.get_string_proteins2string_protein_id()
                logging.info(f"Inserting string 2 pmid mapping")
                string_to_pmid = df_string.join(self.string_proteins2string_protein_id)[["string_protein_id", "pmid"]].drop_duplicates()
                string_to_pmid.to_sql("string_string_protein2pmid", self.conn_mysql, if_exists="append", index=False)          
                logging.info(f"Updating get_pmid_in_db")
                # update set
                self.get_pmid_in_db()

            # PMID:16132081
            logging.info("retrieve len")
            # qseqid	pident	length	mismatch	gapopen	qstart	qend	sstart	send	evalue	bitscore	annotation	preferred_name	publication_id
            df_hit_len = pandas.read_sql(f'select accession,sequence_length as hit_len from stringbbh t1 inner join hash_table t2 on t1.sseqid=t2.accession', self.string_conn_fa).set_index("accession")
            if isinstance(df_hit_len, Series):
                df_hit_len = df_hit_len.to_frame()
            # join dfs to add length 
            logging.info("joining on len")
            df_hits_len = df_hits.reset_index().set_index("sseqid").join(df_hit_len)
            def cal_cov(row):
                return round((row.qend - row.qstart + 1 / self.query2len[row.qseqid])*100, 2)
            
            df_hits_len = df_hits_len.join(self.string_proteins2string_protein_id)              
            df_hits_len["query_cov"] = df_hits_len.apply (lambda row: cal_cov(row), axis=1)
            df_hits_len["hit_cov"] = df_hits_len.apply (lambda row: round((row.send - row.sstart + 1 / row.hit_len)*100, 2), axis=1)
            # use apply to calculate query and hit cov
            #query_cov = round((query_end-query_start + 1) / self.query2len[row[0]]*100, 2)
            #hit_cov = round((hit_end-hit_start + 1)/hit_len*100, 2)
            
            # join with hash2locus
            logging.info("joining on hash2locus")
            df_hits_len = df_hits_len.reset_index().set_index("qseqid").join(self.hash2locus_df)
            # join locus2sefeature_id 
            
            df_hits_len = df_hits_len.reset_index().set_index("locus_tag").join(locus_tag2seqfeature_id_df)
            df_hits_len = df_hits_len.reset_index()[["seqfeature_id", "string_protein_id", "query_cov", "hit_cov", "pident", "evalue", "bitscore"]]
            df_hits_len.columns = ["seqfeature_id", "string_protein_id", "query_cov", "hit_cov", "identity", "evalue", "score"]
            # insert blast result df directly (no neet to loop)
            #sql = f'insert into string_seqfeature_id2string_protein_mapping (seqfeature_id, string_protein_id, query_cov, hit_cov, identity, evalue, score)' \
            #        f' values ({seqfeature_id}, {string_protein_id}, {query_cov}, {hit_cov}, {row.pident}, {row.evalue} ,{row.bitscore})'
            logging.info("Inserting string_seqfeature_id2string_protein_mapping")
            df_hits_len.to_sql("string_seqfeature_id2string_protein_mapping", self.conn_mysql, if_exists="append", index=False)
            
            self.server.adaptor.commit()  


if __name__ == '__main__':
    import argparse
    import chlamdb_setup_utils
    import logging 
    logging.basicConfig(filename="string_load.log", level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

    parser = argparse.ArgumentParser()
    parser.add_argument("-t", '--blast_output', type=str, help="BLAST output file(s)", required=True, nargs="+")
    parser.add_argument("-d", '--db_name', type=str, help="db name", required=True)
    parser.add_argument("-qf", '--query_fasta', type=str, help="db name", required=True)
    parser.add_argument("-dl", '--db_seq_len', type=str, help="string db  acc2len tsv file (len.tsv)", required=True)
    parser.add_argument("-p", '--string_db_name', type=str, help="string sqlite3 db_name", required=True)
    parser.add_argument("-c", '--corresp_table',  type=str, help="hash to locus correspondance table")



    args = parser.parse_args()
    logging.info("Reading hash2locus_df")
    hash2locus_df = pandas.read_csv(args.corresp_table, sep="\t", names=["locus_tag", "hash", "genome"]).set_index("hash")

    pb = StringPMID(args.string_db_name,
                    args.blast_output,
                    args.db_name,
                    hash2locus_df,
                    args.query_fasta,
                    args.db_seq_len)
    pb.load_string_data()
    try:
        pb.index_tables()
    except:
        pass
