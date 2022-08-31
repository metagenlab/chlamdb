#!/usr/bin/env python

import sys
from Bio import AlignIO
import numpy as np
import numpy
from multiprocessing import Process, Queue
import time
from multiprocessing import cpu_count
import os

from chlamdb.biosqldb import manipulate_biosqldb
import pickle
import time

class Orthogroup_Identity_DB:
    def __init__(self, database_name, ncpus):

        self.n_cpus = ncpus

        from chlamdb.biosqldb import manipulate_biosqldb

        self.server, self.db = manipulate_biosqldb.load_db(database_name)
        self.conn = self.server.adaptor.conn
        self.cursor = self.server.adaptor.cursor

        self.count = 0


        sql = 'CREATE TABLE if not exists orthology_identity (' \
                   ' id INTEGER PRIMARY KEY AUTO_INCREMENT,' \
                   ' orthogroup varchar(200),' \
                   ' locus_a VARCHAR(30) NOT NULL,' \
                   ' locus_b VARCHAR(30) NOT NULL,' \
                   ' identity FLOAT(5) ,' \
                   ' align_length INTEGER )'

        self.cursor.execute(sql)

        sql = f'select locus_tag,orthogroup_name from orthology_seqfeature_id2orthogroup t1 ' \
              f' inner join annotation_seqfeature_id2locus t2 on t1.seqfeature_id=t2.seqfeature_id ' \
              f' inner join orthology_orthogroup t3 on t1.orthogroup_id=t3.orthogroup_id'
              
        self.cursor.execute(sql,)
        self.locus_tag2orthogroup = {}
        for row in self.cursor.fetchall():
            self.locus_tag2orthogroup[row[0]] = row[1]

        sql = 'CREATE TABLE if not exists orthology_average_identity (orthogroup varchar(100), identity float(5))'
        self.server.adaptor.execute(sql)       
        self.server.commit()

    def import_alignments(self, cursor, alignment_files):
        import sys
        sys.stdout.write("getting orthogroup identity matrix from %s fasta alignments...\n" % (len(alignment_files)))
        all_matrix = self._get_group_id2identity_matrix(alignment_files, self.n_cpus)
        sys.stdout.write("Creating mysql table to store identity %s matrix...\n" % len(all_matrix))
        for one_matrix in all_matrix.keys():
            #print one_matrix
            #self._create_identity_table(cursor, all_matrix[one_matrix], one_matrix)
            self._add_identity_data(cursor, all_matrix[one_matrix], one_matrix)
            self.conn.commit()


    def index_identity_table(self):
        
        sql1 = 'create index oidla on orthology_identity(locus_a)'
        sql2 = 'create index oidlb on orthology_identity(locus_b)'
        sql3 = 'create index oidlg on orthology_identity(orthogroup)'
        
        self.cursor.execute(sql1,)
        self.cursor.execute(sql2,)
        self.cursor.execute(sql3,)
        

    def _add_identity_data(self, server, group_matrix, group_name):
        '''
        # old approach creating tables of the size n locus vs n locus
        # not adequate for big matrices: sql tables have a max number of columns
        # ==> use two columns tables: locus 1, locus 1, identity
        '''

        locus_list = [i.decode("utf-8") for i in group_matrix[1:,0]]

        for x in range(1, len(group_matrix[:,0])):
            for y in range(x, len(group_matrix[:,0])):
                if group_matrix[x,y] != 0:
                    sql = 'INSERT INTO orthology_identity (orthogroup, locus_a, locus_b, identity) VALUES ("%s", "%s", "%s", %s)'  % (group_name,
                                                                                                                                      locus_list[x-1],
                                                                                                                                      locus_list[y-1],
                                                                                                                                      float(group_matrix[x,y]))
                else:
                    sql = 'INSERT INTO orthology_identity (orthogroup, locus_a, locus_b, identity) VALUES ("%s", "%s", "%s", NULL)'  % (group_name,
                                                                                                                     locus_list[x-1],
                                                                                                                     locus_list[y-1]
                                                                                                                     )
                self.cursor.execute(sql)

    def _chunks(self, l, n):
        import random
        random.shuffle(l)
        return [l[i:i+n] for i in range(0, len(l), n)]

    def _pairewise_identity(self, seq1, seq2):
        import re
        A=list(seq1)
        B=list(seq2)
        identical_sites = 0
        aligned_sites = 0
        gaps=0
        for n in range(0, len(A)):
            if A[n] != "-" and B[n] != "-":
                aligned_sites+=1
            else:
                continue
            if A[n]==B[n]:
                identical_sites+=1

        # update 08.16
        # calculate coverage
        # chose to cut coverage at 30%
        seq1_no_gap = re.sub('-','', str(seq1.seq))
        seq2_no_gap = re.sub('-','', str(seq2.seq))

        if float(aligned_sites)/len(seq1_no_gap) < 0.3:
            return 0
        elif float(aligned_sites)/len(seq2_no_gap) < 0.3:
            return 0
        else:
            return 100*(identical_sites/float(aligned_sites))


    def _get_identity_matrix_from_multiple_alignment(self, alignment):
        identity_matrix = np.chararray((len(alignment)+1, len(alignment)+1), itemsize=30)
        identity_matrix[0, 0] = "-"
        for x in range(0, len(alignment)):
            identity_matrix[x+1, 0] = alignment[x].name
            identity_matrix[0, x+1] = alignment[x].name

            for y in range(x, len(alignment)):
                identity = self._pairewise_identity(alignment[x], alignment[y])
                identity_matrix[y+1, x+1] = round(identity, 2)
                identity_matrix[x+1, y+1] = round(identity, 2)
        return identity_matrix


    def _group_id(self, group_align_files, out_q, list_id):
        outdict = {}
        for i, align_file in enumerate(group_align_files):
            align = AlignIO.read(align_file, "fasta")
            id_matrix = self._get_identity_matrix_from_multiple_alignment(align)
            group_name = self.locus_tag2orthogroup[align[0].id] #os.path.basename(align_file).split(".")[0]
            outdict[group_name] = id_matrix
            self.count += 1
            if i % 20 == 0:
                print("List %s: %s / %s" % (list_id, i, len(group_align_files)))
        pickle.dump(outdict, open("list_%s.p" % list_id, "wb"))
        time.sleep(2)
        out_q.put("list_%s.p" % list_id)

    def _get_group_id2identity_matrix(self, alignments, n_cpus):

        out_q = Queue(n_cpus)

        n_cpu = n_cpus
        n_poc_per_list = int(numpy.ceil(len(alignments)/float(n_cpu)))
        query_lists = self._chunks(alignments, n_poc_per_list)
        #print query_lists
        procs = []
        print("starting... %s prallel jobs" % n_cpu)
        for n, one_list in enumerate(query_lists):
            print("list", n, len(one_list), "elements")
            proc = Process(target=self._group_id, args=(one_list, out_q, n))
            procs.append(proc)
            proc.start()

        # Collect all results into a single result dict. We know how many dicts
        # with results to expect.
        group_id2identity_matrix = {}
        for i in range(n_cpu):
            group_id2identity_matrix.update(pickle.load(open(out_q.get(),"rb")))
        #print "join proc"
        #print "n groups:", len(group_id2identity_matrix.keys())
        time.sleep(5)


        # Wait for all worker processes to finish
        for proc in procs:
            proc.join()

        return group_id2identity_matrix


    def _get_average_identity_from_identity_matrix(self, id_matrix):
        # get upper part of the matrix indexes
        indexes = np.triu_indices(len(id_matrix))
        # calculate mean
        # remove "False" elements (no alignments between some pairs of the orthogroup
        mask = np.invert(np.in1d(id_matrix[indexes], [False]))
        all_id = id_matrix[indexes][mask]
        return np.mean(all_id)


    def add_average_orthogroup_identity(self):

        sql = f'select orthogroup from comparative_tables_orthology'
                   
        groups = [i[0] for i in self.server.adaptor.execute_and_fetchall(sql, )]
        #print len(groups)

        for group in groups:
            id_values = self.get_orthogroup_identity_table(group)
            if len(id_values) > 0:
                av_id = round(np.mean(id_values), 2)
                
                sql = 'insert into orthology_average_identity values ("%s", %s)' % (group, av_id)

                self.server.adaptor.execute(sql)
        self.server.commit()


    def get_orthogroup_identity_table(self, orthogroup):
        import os

        sql = 'SELECT identity FROM orthology_identity where orthogroup="%s"' % orthogroup

        values = [float(i[0]) for i in self.server.adaptor.execute_and_fetchall(sql,)]
        return values



if __name__ == '__main__':
    import argparse
    from chlamdb.biosqldb import manipulate_biosqldb
    
    parser = argparse.ArgumentParser()

    parser.add_argument("-a",'--align_files', type=str, help="aliphment files", nargs='+')
    parser.add_argument("-d",'--db_name', type=str, help="database_name")
    parser.add_argument("-c",'--cpus', type=int, help="n cpus")

    args = parser.parse_args()

    tata = Orthogroup_Identity_DB(args.db_name, args.cpus)
    print("importing alignments...")
    
    tata.import_alignments(tata.cursor, args.align_files)  
    #print("Index table")
    #try:
    #    tata.index_identity_table()
    #except:
    #    pass
    print("add_average_orthogroup_identity")
    tata.add_average_orthogroup_identity()
    
    manipulate_biosqldb.update_config_table(args.db_name, "orthogroup_alignments")
    
