#!/usr/bin/python

from BioSQL import BioSeqDatabase
import sys
from Bio import AlignIO
import numpy as np
import numpy
from multiprocessing import Process, Queue
import time
from multiprocessing import cpu_count
import os
import MySQLdb
from chlamdb.biosqldb import manipulate_biosqldb

class Orthogroup_Identity_DB:
    def __init__(self, database_name):
        import os
        sqlpsw = os.environ['SQLPSW']
        try:
            self.conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                user="root", # your username
                                passwd=sqlpsw, # your password
                                db="orth_%s" % database_name) # name of the data base
            self.cursor = self.conn.cursor()
        except:
            conn = MySQLdb.connect(host="localhost", user="root", passwd=sqlpsw)
            sql = 'CREATE DATABASE orth_%s' % database_name
            cursor = conn.cursor()
            cursor.execute(sql)
            self.conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                user="root", # your username
                                passwd=sqlpsw, # your password
                                db="orth_%s" % database_name) # name of the data base
            self.cursor = self.conn.cursor()
        self.conn.commit()

    def import_alignments(self, cursor, alignment_files):

        import sys
        sys.stdout.write("getting orthogroup identity matrix from %s fasta alignments...\n" % (len(alignment_files)))
        all_matrix = self._get_group_id2identity_matrix(alignment_files)
        sys.stdout.write("Creating mysql table to store identity %s matrix...\n" % len(all_matrix))
        for one_matrix in all_matrix.keys():
            #print one_matrix
            #self._create_identity_table(cursor, all_matrix[one_matrix], one_matrix)
            self._add_identity_data(cursor, all_matrix[one_matrix], one_matrix)
            self.conn.commit()

    '''
    def _create_identity_table(self, server, group_matrix, group_name):



        # old table scheme with matrix style storage
        # unpossible for big matrices
        # ==> data now stored in the form: locus_a, locus_b, percentage_identity

        template_cols = ''

        for i in range(1, len(group_matrix[:,0])-1):
            template_cols += '`%s` FLOAT(5), ' % group_matrix[i,0]
        template_cols += '`%s` FLOAT(5)' % group_matrix[-1,0]

        sql = 'CREATE TABLE %s (' \
                   ' id INT(6) UNSIGNED AUTO_INCREMENT PRIMARY KEY,' \
                   ' locus_tag VARCHAR(30) NOT NULL,' \
                   ' %s )' % (group_name, template_cols)
        try:
            server.execute(sql)
        except:
            print group_name
            print group_matrix
            print sql
            exit()
    '''



    def _add_identity_data(self, server, group_matrix, group_name):
        '''
        # old approach creating tables of the size n locus vs n locus
        # not adequate for big matrices: sql tables have a max number of columns
        # ==> use two columns tables: locus 1, locus 1, identity

        if len(group_matrix) == 2:
            sql = 'INSERT INTO %s (locus_tag, `%s`) VALUES ("%s", "%s")' % (group_name, group_matrix[0, 1], group_matrix[0, 1], group_matrix[1, 1])
            server.execute(sql)

            return

        template_cols = 'locus_tag, '
        for i in range(1, len(group_matrix[:,0])-1):
            template_cols += '`%s`, ' % group_matrix[i,0]
        template_cols += '`%s`' % group_matrix[-1,0]



        for row in range(1, len(group_matrix)):
            template=''
            for i in range(0, len(group_matrix)-1):
                template += '"%s", ' % group_matrix[i, row]
            template += str(group_matrix[-1, row])
            sql = 'INSERT INTO %s (%s) VALUES (%s)' % (group_name, template_cols, template)
            server.execute(sql)
        '''
        sql = 'CREATE TABLE %s (' \
                   ' id INT(6) UNSIGNED AUTO_INCREMENT PRIMARY KEY,' \
                   ' locus_a VARCHAR(30) NOT NULL,' \
                   ' locus_b VARCHAR(30) NOT NULL,' \
                   ' identity FLOAT(5) )' % (group_name)
        server.execute(sql)

        locus_list = group_matrix[1:,0]
        #print locus_list
        for x in range(1, len(group_matrix[:,0])):
            for y in range(x, len(group_matrix[:,0])):
                if group_matrix[x,y] != 0:
                    sql = 'INSERT INTO %s (locus_a, locus_b, identity) VALUES ("%s", "%s", %s)'  % (group_name,
                                                                                                      locus_list[x-1],
                                                                                                      locus_list[y-1],
                                                                                                      group_matrix[x,y])
                else:
                    sql = 'INSERT INTO %s (locus_a, locus_b, identity) VALUES ("%s", "%s", NULL)'  % (group_name,
                                                                                                      locus_list[x-1],
                                                                                                      locus_list[y-1]
                                                                                                      )
                server.execute(sql)



    def _chunks(self, l, n):
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
        #if aligned_sites == 0:
        #    # return false if align length = 0
        #    return False
        seq1_no_gap = re.sub('-','', str(seq1.seq))
        seq2_no_gap = re.sub('-','', str(seq2.seq))

        if float(aligned_sites)/len(seq1_no_gap) < 0.3:
            return 0
        elif float(aligned_sites)/len(seq2_no_gap) < 0.3:
            return 0
        else:
            return 100*(identical_sites/float(aligned_sites))

    def _get_identity_matrix_from_multiple_alignment(self, alignment):
        # [len(alignment)+1, len(alignment)]
        identity_matrix = np.chararray((len(alignment)+1, len(alignment)+1), itemsize=30)
        identity_matrix[0,0] = "-"
        # first column = locus tags
        for x in range(0,len(alignment)):
            identity_matrix[x+1, 0] = alignment[x].name
            identity_matrix[0, x+1] = alignment[x].name

            for y in range(x, len(alignment)):
                identity = self._pairewise_identity(alignment[x], alignment[y])
                #print "identity", identity
                identity_matrix[y+1, x+1] = round(identity, 2)
                identity_matrix[x+1, y+1] = round(identity, 2)
        #print identity_matrix
        return identity_matrix


    def _group_id(self, group_align_files, out_q):
        outdict = {}
        for align_file in group_align_files:
            align = AlignIO.read(align_file, "fasta")

            id_matrix = self._get_identity_matrix_from_multiple_alignment(align)
            #self.mcl_group_id2identity_matrix[group] = id_matrix
            group_name = os.path.basename(align_file).split(".")[0]
            outdict[group_name] = id_matrix
        out_q.put(outdict)

    def _get_group_id2identity_matrix(self, alignments):

        out_q = Queue()

        n_cpu = 8
        n_poc_per_list = int(numpy.ceil(len(alignments)/float(n_cpu)))
        query_lists = self._chunks(alignments, n_poc_per_list)
        #print query_lists
        procs = []
        for one_list in query_lists:
            proc = Process(target=self._group_id, args=(one_list, out_q))
            procs.append(proc)
            proc.start()

        # Collect all results into a single result dict. We know how many dicts
        # with results to expect.
        group_id2identity_matrix = {}
        for i in range(n_cpu):
            group_id2identity_matrix.update(out_q.get())
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

    def _create_orthogroup_average_identity_column(self, server, biodatabase_name):


        sql = 'CREATE TABLE orth_%s.average_identity (orthogroup varchar(100), identity float(5))' % biodatabase_name
        server.adaptor.execute(sql)
        server.commit()

    def add_average_orthogroup_identity(self, biodatabase_name):

        #print 'adding average id to orthology table'
        server, db = manipulate_biosqldb.load_db(biodatabase_name)

        #print 'get orthogroups'
        sql = 'select orthogroup from comparative_tables.orthology_%s' % biodatabase_name
        try:
            #print 'adding column'
            self._create_orthogroup_average_identity_column(server, biodatabase_name)
        except:
            print ("column already created?")
        groups = [i[0] for i in server.adaptor.execute_and_fetchall(sql, )]
        #print len(groups)

        for group in groups:
            #print "group %s" % group
            try:
                id_table = np.array(get_orthogroup_identity_table(biodatabase_name, group))
                id_matrix = id_table[:,1:].astype(float)
                #print np.mean(id_matrix)
                #print self._get_average_identity_from_identity_matrix(id_matrix)
                av_id = round(np.mean(id_matrix[np.triu_indices(len(id_matrix), k=1)]), 2)
                #print av_id
            except:
                av_id = 0

            sql = 'insert into orth_%s.average_identity values ("%s", %s)' % (biodatabase_name, group, av_id)
            #print sql
            server.adaptor.execute(sql)
            server.commit()



            #from tabulate import tabulate
            #print tabulate(id_table)





def check_identity(cursor, orthogroup, locus1, locus2):
    #conn = MySQLdb.connect(host="localhost", # your host, usually localhost
    #                            user="root", # your username
    #                            passwd=, # your password
    #                            db="orth_%s" % biodb_name) # name of the data base
    #cursor = conn.cursor()

    sql = 'select identity from %s where (locus_a = "%s" and locus_b = "%s") or locus_a ="%s" and locus_b="%s"'  % (orthogroup,locus1, locus2, locus1, locus2)
    #print sql
    cursor.execute(sql)
    identity = int(cursor.fetchone()[0])
    return identity


def get_orthogroup_identity_table(biodb_name, orthogroup):
    import os
    sqlpsw = os.environ['SQLPSW']

    conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                user="root", # your username
                                passwd=sqlpsw, # your password
                                db="orth_%s" % biodb_name) # name of the data base
    cursor = conn.cursor()

    sql = 'SELECT * FROM %s' % orthogroup

    cursor.execute(sql)

    rows = cursor.fetchall()

    return [row[1:] for row in rows]




def get_identity_matrix_all_genomes(orthogroup, biodb_namd):
    identity_table = get_orthogroup_identity_table("saureus_01_15", "group_444")
    #for row in identity_table:
    #    print row




def locus_tag2identity_best_hit_all_genomes(biodb_name, locus, group_name, locus_tag2taxonomic_id_dict = ""):

    server, db = manipulate_biosqldb.load_db(biodb_name)

    genomes = manipulate_biosqldb.get_genome_taxons_list(server, biodb_name)

    #print genomes

    # create empty dictionnary
    genome2genome = {}
    for genome1 in genomes:
        if not genome1 in genome2genome.keys():
            genome2genome[str(genome1)] = {}
            for genome2 in genomes:
                genome2genome[str(genome1)][str(genome2)] = []

    try:
        identity_table = np.array(get_orthogroup_identity_table(biodb_name, group_name))
        all_locus = identity_table[:, 0]
    except:
        genome2best_hit = {}
        for genome1 in genomes:
            genome2best_hit[genome1] = 0
        return genome2best_hit

    #print "identity table", len(identity_table)

    #locus_tag2genome_name = manipulate_biosqldb.locus_tag2genome_description(server, biodb_name)
    locus_tag2taxonomic_id_dict = manipulate_biosqldb.locus_tag2genome_taxon_id(server, biodb_name)

    #print "locus_tag2taxonomic_id ok"


    locus_index_ref = list(all_locus).index(locus)

    '''
    for y in range(1, len(all_locus) + 1):
        genome1 = locus_tag2taxonomic_id_dict[all_locus[y-1]]
        for z in range(0, len(all_locus)):
            genome2= locus_tag2taxonomic_id_dict[all_locus[z]]
            genome2genome[str(genome1)][str(genome2)].append(identity_table[z, y])
    '''
    genome2best_hit = {}
    for genome1 in genomes:
        genome2best_hit[genome1] = 0
    #print genome2best_hit
    for target_locus in all_locus:
        locus_index_target = list(all_locus).index(target_locus)
        identity = identity_table[locus_index_ref, locus_index_target + 1]
        target_taxon = locus_tag2taxonomic_id_dict[target_locus]
        if float(identity) > genome2best_hit[str(target_taxon)]:
            genome2best_hit[str(target_taxon)] = float(identity)
    #print genome2best_hit
    return genome2best_hit




def locus_list2identity_in_other_genomes(locus_list, biodb):
    server, db = manipulate_biosqldb.load_db(biodb)

    locus_tag2seqfeature_id = manipulate_biosqldb.locus_tag2seqfeature_id_dict(server, biodb)

    taxon_id2description = manipulate_biosqldb.taxon_id2genome_description(server, biodb)

    import re
    for i in taxon_id2description.keys():
        taxon_id2description[i] = re.sub(" subsp\. aureus", "", taxon_id2description[i])
        taxon_id2description[i] = re.sub(", complete genome\.", "", taxon_id2description[i])
        taxon_id2description[i] = re.sub(", complete sequence\.", "", taxon_id2description[i])
        taxon_id2description[i] = re.sub("strain ", "", taxon_id2description[i])
        taxon_id2description[i] = re.sub("str\. ", "", taxon_id2description[i])
        taxon_id2description[i] = re.sub(" complete genome sequence\.", "", taxon_id2description[i])
        taxon_id2description[i] = re.sub(" complete genome\.", "", taxon_id2description[i])
        taxon_id2description[i] = re.sub(" chromosome", "", taxon_id2description[i])
        taxon_id2description[i] = re.sub("Staphylococcus aureus ", "", taxon_id2description[i])


    header = 'orthogroup\t'
    dico = locus_tag2identity_best_hit_all_genomes(biodb, 'wcw_1594', 'group_417')
    for i in dico.keys():

        header += taxon_id2description[i] + '\t'

    final_out = header + '\n'

    for locus in locus_list:
        #print "locus", i
        seqfeature_id = locus_tag2seqfeature_id[locus]
        orthogroup = manipulate_biosqldb.seqfeature_id2orthogroup(server, seqfeature_id, biodb)
        #print "ortho", orthogroup
        dico = locus_tag2identity_best_hit_all_genomes(biodb, locus, orthogroup)
        #print "dico done..."
        out = '%s\t' % orthogroup
        for i in dico.keys():
            identity = dico[i]
            out += '%s\t' % identity
        final_out += out + '\n'

    return final_out


def heatmap_presence_absence(biodb_name, group_name):
    server, db = manipulate_biosqldb.load_db(biodb_name)

    genomes = manipulate_biosqldb.get_genome_taxons_list(server, biodb_name)

    template = ''
    for i in range(0, len(genomes)):
        template += '`%s`, ' % genomes[i]
    template += '`%s`' % genomes[-1]

    #print "template", template

    sql = 'select %s from orthology_%s where orthogroup = "%s"' % (template, biodb_name, group_name)

    result = [int(i) for i in server.adaptor.execute_and_fetchall(sql,)[0]]
    #print result
    taxon2presence_absence = {}
    for x, y in zip(genomes, result):
        taxon2presence_absence[x] = y
    return taxon2presence_absence


def orthogroup2identity_dico(biodb_name, orthogroup):
    #server, db = manipulate_biosqldb.load_db(biodb_name)
    identity_table = np.array(get_orthogroup_identity_table(biodb_name, orthogroup))
    locus2locus_id = {}
    for row in identity_table:
        if row[0] not in locus2locus_id:
            locus2locus_id[row[0]] = {}
            locus2locus_id[row[0]][row[1]] = row[2]
        else:
            locus2locus_id[row[0]][row[1]] = row[2]
        if row[1] not in locus2locus_id:
            locus2locus_id[row[1]] = {}
            locus2locus_id[row[1]][row[0]] = row[2]
        else:
            locus2locus_id[row[1]][row[0]] = row[2]
    return locus2locus_id
    '''
    # was made for square identity matrix. Changed to columns data (mai 2016).
    identity_table = np.array(get_orthogroup_identity_table(biodb_name, orthogroup))
    print identity_table
    all_locus = identity_table[:,0]
    print "all locus", all_locus
    locus2locus_id = {}
    for counter_1, locus_1 in enumerate(all_locus):
        for counter_2, locus_2 in enumerate(all_locus):
            if locus_1 not in locus2locus_id:
                locus2locus_id[locus_1] = {}
                locus2locus_id[locus_1][locus_2] = identity_table[counter_1, counter_2+1]
            else:
                locus2locus_id[locus_1][locus_2] = identity_table[counter_1, counter_2+1]
    return locus2locus_id
    '''

def locus_list2presence_absence_all_genomes(locus_list, biodb_name):
    server, db = manipulate_biosqldb.load_db(biodb_name)

    locus_tag2seqfeature_id = manipulate_biosqldb.locus_tag2seqfeature_id_dict(server, biodb_name)

    taxon_id2description = manipulate_biosqldb.taxon_id2genome_description(server, biodb_name)

    import re
    for i in taxon_id2description.keys():
        taxon_id2description[i] = re.sub(" subsp\. aureus", "", taxon_id2description[i])
        taxon_id2description[i] = re.sub(", complete genome\.", "", taxon_id2description[i])
        taxon_id2description[i] = re.sub(", complete sequence\.", "", taxon_id2description[i])
        taxon_id2description[i] = re.sub("strain ", "", taxon_id2description[i])
        taxon_id2description[i] = re.sub("str\. ", "", taxon_id2description[i])
        taxon_id2description[i] = re.sub(" complete genome sequence\.", "", taxon_id2description[i])
        taxon_id2description[i] = re.sub(" complete genome\.", "", taxon_id2description[i])
        taxon_id2description[i] = re.sub(" chromosome", "", taxon_id2description[i])
        taxon_id2description[i] = re.sub("Staphylococcus aureus ", "", taxon_id2description[i])


    header = 'orthogroup\t'
    genomes = manipulate_biosqldb.get_genome_taxons_list(server, biodb_name)
    for i in genomes:
        header += taxon_id2description[i] + '\t'
    final_out = header + '\n'

    for i in locus_list:
        #print "locus", i
        seqfeature_id = locus_tag2seqfeature_id[i]
        orthogroup = manipulate_biosqldb.seqfeature_id2orthogroup(server, seqfeature_id, biodb_name)
        #print "ortho", orthogroup
        dico = heatmap_presence_absence(biodb_name, orthogroup)

        #print "dico done..."
        #print dico
        out = '%s\t' % orthogroup
        for i in genomes:

            out += '%s\t' % dico[i]
        final_out +=out + '\n'

    return final_out



def orthogroup2average_identity(biodatabase_name):
    server, db = manipulate_biosqldb.load_db(biodatabase_name)
    sql = 'select orthogroup,average_identity from orth_%s.average_identity' % biodatabase_name
    result = server.adaptor.execute_and_fetchall(sql)
    return manipulate_biosqldb.to_dict(result)





if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument("-a",'--align_files', type=str, help="aliphment files", nargs='+')
    parser.add_argument("-d",'--db_name', type=str, help="database_name")

    args = parser.parse_args()

    tata = Orthogroup_Identity_DB(args.db_name)
    #print "importing alignments..."
    tata.import_alignments(tata.cursor, args.align_files)
    tata.add_average_orthogroup_identity(args.db_name)
    #check_identity("Chlamydia_12_14", "group_825", "Cav1_00733", "CT565")
