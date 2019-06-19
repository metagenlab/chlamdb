#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


def interpro2biosqlV2(server,
                      hash2locus_list,
                      interpro_accession2interpro_id,
                      seqfeature_id2locus_tag,
                      locus_tag2genome_taxon_id,
                      protein_id2genome_taxon_id,
                      locus_tag2seqfeature_id,
                      protein_id2seqfeature_id,
                      seqfeature_id2organism,
                      db_name,
                      *input_files):

    import re
    '''
    1. Protein Accession (e.g. P51587)
    2. Sequence MD5 digest (e.g. 14086411a2cdf1c4cba63020e1622579)
    3.. equence Length (e.g. 3418)
    4. Analysis (e.g. Pfam / PRINTS / Gene3D)
    5. Signature Accession (e.g. PF09103 / G3DSA:2.40.50.140)
    6. Signature Description (e.g. BRCA2 repeat profile)
    7. Start location
    8. Stop location
    9. Score - is the e-value of the match reported by member database method (e.g. 3.1E-52)
    10. Status - is the status of the match (T: true)
    11. Date - is the date of the run
    12. (InterPro annotations - accession (e.g. IPR002093) - optional column; only displayed if -iprscan option is switched on)
    13. (InterPro annotations - description (e.g. BRCA2 repeat) - optional column; only displayed if -iprscan option is switched on)
    14. (GO annotations (e.g. GO:0005515) - optional column; only displayed if --goterms option is switched on)
    15. (Pathways annotations (e.g. REACT_71) - optional column; only displayed if --pathways option is switched on)
    :param input_file:
    :return:
    '''

    sql = 'CREATE TABLE interpro.interpro_%s (accession VARCHAR(100),' \
          ' seqfeature_id INT, ' \
          ' organism VARCHAR(200),  ' \
          ' taxon_id INT,' \
          ' sequence_length INT, ' \
          ' analysis VARCHAR(100) NOT NULL, ' \
          ' signature_accession VARCHAR(100), ' \
          ' signature_description VARCHAR(1000), ' \
          ' start INT, ' \
          ' stop INT, ' \
          ' score VARCHAR(20) NOT NULL, ' \
          ' interpro_id INT, ' \
          ' GO_terms varchar(10000),' \
          ' pathways varchar(10000),' \
          ' INDEX seqfeature_id (seqfeature_id),' \
          ' INDEX interpro_id (interpro_id),' \
          ' INDEX ia (interpro_accession))' % db_name
    try:
        server.adaptor.execute(sql)
    except:
        pass
    for one_interpro_file in input_files:
        print(one_interpro_file)
        from pandas import read_csv
        with open(one_interpro_file, 'r') as f:
            tsvin = read_csv(f, sep='\t', error_bad_lines=False, names=["accession",
                                                                        "MD5",
                                                                        "sequence_length",
                                                                        "analysis",
                                                                        "signature_accession",
                                                                        "signature_description",
                                                                        "start",
                                                                        "stop",
                                                                        "score",
                                                                        "status",
                                                                        "date",
                                                                        "interpro_accession",
                                                                        "interpro_description",
                                                                        "GO_terms",
                                                                        "pathways"])
            tsvin = tsvin.fillna(0)

            for i in range(len(tsvin['accession'])):

                data = list(tsvin.loc[i,:])
                for index, item in enumerate(data):
                    if type(item) == str:
                        data[index] = item.replace('\"','')

                locus_tag_list = hash2locus_list[data[0]]

                sequence_length = data[2]
                analysis = data[3]
                signature_accession = data[4]
                signature_description = data[5]
                start = data[6]
                stop = data[7]
                score = data[8]


                interpro_accession = data[11]
                interpro_description = data[12]
                GO_terms = data[13]
                pathways = data[14]

                for locus_tag in locus_tag_list:
                    taxon_id = locus_tag2genome_taxon_id[locus_tag]
                    seqfeature_id = locus_tag2seqfeature_id[locus_tag]
                    organism = seqfeature_id2organism[str(seqfeature_id)]

                    sql = 'INSERT INTO interpro.interpro_%s(accession, locus_tag, organism, taxon_id,' \
                          ' sequence_length, analysis, signature_accession, signature_description, start, ' \
                          ' stop, score, interpro_accession, interpro_description, GO_terms, pathways) ' \
                          ' values ("%s", "%s", "%s", %s, %s, "%s", "%s", "%s", %s, %s, "%s", "%s", "%s", "%s", "%s");' % (db_name,
                                                                                                                            locus_tag,
                                                                                                                            seqfeature_id,
                                                                                                                            organism,
                                                                                                                            taxon_id,
                                                                                                                            sequence_length,
                                                                                                                            analysis,
                                                                                                                            signature_accession,
                                                                                                                            signature_description,
                                                                                                                            int(start),
                                                                                                                            int(stop),
                                                                                                                            str(score),
                                                                                                                            interpro_id,
                                                                                                                            GO_terms,
                                                                                                                            pathways)
                    try:
                        server.adaptor.execute(sql)
                        server.adaptor.commit()
                    except:
                        print(sql)
                        print(data)
                        import sys
                        sys.exit()


def update_analysis_dico(server):

    from chlamdb.biosqldb import manipulate_biosqldb

    sql = 'select analysis_name, analysis_id from interpro.analysis'

    analysis_nam2analysis_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    return analysis_nam2analysis_id

def interpro2biosql(server,
                    hash2locus,
                    locus_tag2seqfeature_id,
                    db_name,
                    *input_files):

    import MySQLdb
    '''
    1. Protein Accession (e.g. P51587)
    2. Sequence MD5 digest (e.g. 14086411a2cdf1c4cba63020e1622579)
    3.. equence Length (e.g. 3418)
    4. Analysis (e.g. Pfam / PRINTS / Gene3D)
    5. Signature Accession (e.g. PF09103 / G3DSA:2.40.50.140)
    6. Signature Description (e.g. BRCA2 repeat profile)
    7. Start location
    8. Stop location
    9. Score - is the e-value of the match reported by member database method (e.g. 3.1E-52)
    10. Status - is the status of the match (T: true)
    11. Date - is the date of the run
    12. (InterPro annotations - accession (e.g. IPR002093) - optional column; only displayed if -iprscan option is switched on)
    13. (InterPro annotations - description (e.g. BRCA2 repeat) - optional column; only displayed if -iprscan option is switched on)
    14. (GO annotations (e.g. GO:0005515) - optional column; only displayed if --goterms option is switched on)
    15. (Pathways annotations (e.g. REACT_71) - optional column; only displayed if --pathways option is switched on)
    :param input_file:
    :return:
    '''

    sql = 'CREATE TABLE if not exists interpro.analysis (analysis_id INT AUTO_INCREMENT PRIMARY KEY, ' \
          ' analysis_name varchar(400),' \
          ' index analysis_name(analysis_name))'

    server.adaptor.execute(sql,)

    sql2 = 'CREATE TABLE if not exists interpro.signature (signature_id INT AUTO_INCREMENT PRIMARY KEY, ' \
           ' signature_accession varchar(400),' \
           ' signature_description TEXT,' \
           ' analysis_id INT,' \
           ' interpro_id INT, ' \
           ' GO_terms TEXT,' \
           ' pathways TEXT,' \
           ' INDEX analysis_id (analysis_id),' \
           ' index signature_accession(signature_accession),' \
           ' index interpro_id(interpro_id))' \

    server.adaptor.execute(sql2,)

    sql3 = 'CREATE TABLE if not exists interpro.interpro_%s (seqfeature_id INT,' \
          ' sequence_length INT, ' \
          ' signature_id INT, ' \
          ' start INT, ' \
          ' stop INT, ' \
          ' score TEXT,' \
          ' INDEX signature_id (signature_id),' \
          ' INDEX seqfeature(seqfeature_id))' % db_name

    server.adaptor.execute(sql3)

    analysis2analysis_id = update_analysis_dico(server)
    sql = 'select signature_accession, signature_id from interpro.signature'
    signature2signature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
    sql = 'select name, interpro_id from interpro.entry'
    interpro_entry2interpro_entry_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    for one_interpro_file in input_files:
        print(one_interpro_file)
        from pandas import read_csv
        with open(one_interpro_file, 'r') as f:
            tsvin = read_csv(f, sep='\t', error_bad_lines=False, names=["accession",
                                                                        "MD5",
                                                                        "sequence_length",
                                                                        "analysis",
                                                                        "signature_accession",
                                                                        "signature_description",
                                                                        "start",
                                                                        "stop",
                                                                        "score",
                                                                        "status",
                                                                        "date",
                                                                        "interpro_accession",
                                                                        "interpro_description",
                                                                        "GO_terms",
                                                                        "pathways"])
            tsvin = tsvin.fillna(0)

            for i in range(len(tsvin['accession'])):

                data = list(tsvin.loc[i,:])
                for index, item in enumerate(data):
                    if type(item) == str:
                        data[index] = item.replace('\"','')

                locus_tag_list = hash2locus_list[data[0]]

                sequence_length = data[2]
                analysis = data[3]
                signature_description = data[5]
                interpro_accession = data[11]
                interpro_description = data[12]
                GO_terms = data[13]
                pathways = data[14]

                try:
                    analysis_id = analysis2analysis_id[analysis]
                except KeyError:
                    print ('New analysis:', analysis)
                    sql = 'insert into interpro.analysis (analysis_name) values ("%s")' % analysis
                    server.adaptor.execute(sql)
                    server.adaptor.commit()
                    analysis2analysis_id = update_analysis_dico(server)
                    analysis_id = analysis2analysis_id[analysis]

                signature_accession = data[4]

                #sql = 'select signature_id from interpro.signature where signature_accession="%s"' % signature_accession

                try:
                    signature_id = signature2signature_id[signature_accession]#server.adaptor.execute_and_fetchall(sql,)[0][0]
                except KeyError:
                    print ('New signature', signature_accession, signature_description)
                    #sql1 = 'select interpro_id from interpro.entry where name="%s"' % (interpro_accession)

                    try:
                        interpro_id = interpro_entry2interpro_entry_id[interpro_accession]#server.adaptor.execute_and_fetchall(sql1,)[0][0]
                    except KeyError:
                        if interpro_accession == 0:
                            print('No interpro-accession for ', signature_accession, signature_description)
                            interpro_id="NULL"
                        else:
                            print('New Interpro entry', interpro_accession, interpro_description)

                            sql1b = 'insert into interpro.entry(name, description) values("%s","%s")' % (interpro_accession,
                                                                                                         interpro_description)
                            print (sql1b)
                            server.adaptor.execute(sql1b,)
                            server.adaptor.commit()
                            sql1 = 'select interpro_id from interpro.entry where name="%s"' % (interpro_accession)
                            interpro_id = server.adaptor.execute_and_fetchall(sql1,)[0][0]
                            # update dictionnray
                            interpro_entry2interpro_entry_id[interpro_accession] = interpro_id

                    sql2 = 'insert into interpro.signature (signature_accession, signature_description, ' \
                          ' analysis_id, interpro_id, GO_terms, pathways) values ("%s", "%s", %s, %s, "%s", "%s")' % (signature_accession,
                                                                                                     signature_description,
                                                                                                     analysis_id,
                                                                                                     interpro_id,
                                                                                                     GO_terms,
                                                                                                     pathways)

                    server.adaptor.execute(sql2,)
                    server.adaptor.commit()
                    sql = 'select signature_id from interpro.signature where signature_accession="%s"' % signature_accession
                    signature_id = server.adaptor.execute_and_fetchall(sql,)[0][0]
                    # update dictionnary
                    signature2signature_id[signature_accession] = signature_id

                start = data[6]
                stop = data[7]
                score = data[8]

                if analysis in ['Phobius', 'Coils', 'SignalP-TM', 'SignalP_EUK', 'ProSitePatterns', 'SignalP_GRAM_NEGATIVE', 'SignalP_GRAM_POSITIVE']:
                    score = "NULL"

                for locus_tag in locus_tag_list:
                    seqfeature_id = locus_tag2seqfeature_id[locus_tag]

                    sql = 'INSERT INTO interpro.interpro_%s(seqfeature_id,' \
                          ' signature_id,' \
                          ' sequence_length, ' \
                          ' start, ' \
                          ' stop, ' \
                          ' score) ' \
                          ' values (%s, %s, %s, %s, %s, "%s");' % (db_name,
                                                                   seqfeature_id,
                                                                   signature_id,
                                                                   sequence_length,
                                                                   int(start),
                                                                   int(stop),
                                                                   score)
                    try:
                        server.adaptor.execute(sql)
                    except:
                        print(sql)
                        print(data)
                        import sys
                        sys.exit()
            server.adaptor.commit()


def add_TM_and_SP_columns(db_name):

    server, db = manipulate_biosqldb.load_db(db_name)

    sql = 'ALTER TABLE orthology_detail_%s ADD COLUMN SP INT AFTER seqfeature_id;' % db_name
    #server.adaptor.execute(sql,)
    sql = 'ALTER TABLE orthology_detail_%s ADD COLUMN TM BOOLEAN not null default 0 AFTER SP;' % db_name
    #server.adaptor.execute(sql,)

    sql = 'select seqfeature_id, count(*) as n from interpro.interpro_%s t1 inner join interpro.signature t2 on t1.signature_id=t2.signature_id inner join interpro.analysis t3 on t2.analysis_id=t3.analysis_id where analysis_name="Phobius" and signature_accession="TRANSMEMBRANE" group by seqfeature_id;' % db_name
    sql2 = 'select seqfeature_id, count(*) as n from interpro.interpro_%s t1 inner join interpro.signature t2 on t1.signature_id=t2.signature_id inner join interpro.analysis t3 on t2.analysis_id=t3.analysis_id where analysis_name="Phobius" and signature_accession="SIGNAL_PEPTIDE" group by seqfeature_id;' % db_name

    seqfeature_id2TM = server.adaptor.execute_and_fetchall(sql,)
    for row in seqfeature_id2TM:
        sql = 'update orthology_detail_%s set SP=%s where seqfeature_id=%s' % (db_name, 1, row[0])
        server.adaptor.execute(sql,)
    seqfeature_id2signal_peptide = server.adaptor.execute_and_fetchall(sql2,)
    for row in seqfeature_id2TM:
        sql = 'update orthology_detail_%s set TM=%s where seqfeature_id=%s' % (db_name, row[1], row[0])
        server.adaptor.execute(sql,)

    server.commit()


def interpro2biosql_legacy(server,
                           seqfeature_id2locus_tag,
                           locus_tag2genome_taxon_id,
                           protein_id2genome_taxon_id,
                           locus_tag2seqfeature_id,
                           protein_id2seqfeature_id,
                           seqfeature_id2organism,
                           db_name,
                           seqfeature_id2orthogroup,
                           hash2locus_list,
                           *input_files):
    import re
    '''
    1. Protein Accession (e.g. P51587)
    2. Sequence MD5 digest (e.g. 14086411a2cdf1c4cba63020e1622579)
    3.. equence Length (e.g. 3418)
    4. Analysis (e.g. Pfam / PRINTS / Gene3D)
    5. Signature Accession (e.g. PF09103 / G3DSA:2.40.50.140)
    6. Signature Description (e.g. BRCA2 repeat profile)
    7. Start location
    8. Stop location
    9. Score - is the e-value of the match reported by member database method (e.g. 3.1E-52)
    10. Status - is the status of the match (T: true)
    11. Date - is the date of the run
    12. (InterPro annotations - accession (e.g. IPR002093) - optional column; only displayed if -iprscan option is switched on)
    13. (InterPro annotations - description (e.g. BRCA2 repeat) - optional column; only displayed if -iprscan option is switched on)
    14. (GO annotations (e.g. GO:0005515) - optional column; only displayed if --goterms option is switched on)
    15. (Pathways annotations (e.g. REACT_71) - optional column; only displayed if --pathways option is switched on)
    :param input_file:
    :return:
    '''

    sql = 'CREATE TABLE biosqldb.interpro_%s (accession VARCHAR(100),' \
          ' locus_tag VARCHAR(200), ' \
          ' organism VARCHAR(400),  ' \
          ' taxon_id INT,' \
          ' sequence_length INT, ' \
          ' analysis VARCHAR(100) NOT NULL, ' \
          ' signature_accession VARCHAR(100), ' \
          ' signature_description TEXT, ' \
          ' start INT, ' \
          ' stop INT, ' \
          ' score VARCHAR(20) NOT NULL, ' \
          ' interpro_accession varchar(400) NOT NULL, ' \
          ' interpro_description TEXT,' \
          ' GO_terms TEXT,' \
          ' pathways TEXT,' \
          ' orthogroup varchar(400),' \
          ' seqfeature_id INT,' \
          ' INDEX loc (locus_tag),' \
          ' INDEX ana (analysis),' \
          ' INDEX tid (taxon_id),' \
          ' INDEX og (organism),' \
          ' INDEX acc (accession),' \
          ' INDEX sid (seqfeature_id),' \
          ' INDEX ia (interpro_accession))' % db_name
    print(sql)
    server.adaptor.execute(sql)
    for one_interpro_file in input_files:
        print (one_interpro_file)
        from pandas import read_csv
        with open(one_interpro_file, 'r') as f:
            tsvin = read_csv(f, sep='\t', error_bad_lines=False, names=["accession",
                                                                        "MD5",
                                                                        "sequence_length",
                                                                        "analysis",
                                                                        "signature_accession",
                                                                        "signature_description",
                                                                        "start",
                                                                        "stop",
                                                                        "score",
                                                                        "status",
                                                                        "date",
                                                                        "interpro_accession",
                                                                        "interpro_description",
                                                                        "GO_terms",
                                                                        "pathways"])
            tsvin = tsvin.fillna(0)

            for i in range(len(tsvin['accession'])):

                data= list(tsvin.loc[i,:])
                for index, item in enumerate(data):
                    if type(item) == str:
                        data[index] = item.replace('\"','')

                accession = data[0]

                sequence_length = data[2]
                analysis = data[3]
                signature_accession = data[4]
                signature_description = data[5]
                start = data[6]
                stop = data[7]
                score = data[8]


                interpro_accession = data[11]
                interpro_description = data[12]
                GO_terms = data[13]
                pathways = data[14]

                for locus_tag in hash2locus_list[accession]:

                    try:
                        taxon_id = protein_id2genome_taxon_id[locus_tag]
                        seqfeature_id = protein_id2seqfeature_id[locus_tag]
                    except KeyError:
                        taxon_id = locus_tag2genome_taxon_id[locus_tag]
                        seqfeature_id = locus_tag2seqfeature_id[locus_tag]
                    organism = seqfeature_id2organism[str(seqfeature_id)]
                    #print organism
                    locus_tag = seqfeature_id2locus_tag[str(seqfeature_id)]
                    orthogroup = seqfeature_id2orthogroup[str(seqfeature_id)]
                    sql = 'INSERT INTO biosqldb.interpro_%s(accession, locus_tag, organism, taxon_id,' \
                          ' sequence_length, analysis, signature_accession, signature_description, start, ' \
                          ' stop, score, interpro_accession, interpro_description, GO_terms, pathways, orthogroup, seqfeature_id) ' \
                          ' values ("%s", "%s", "%s", %s, %s, "%s", "%s", "%s", %s, %s, "%s", "%s", "%s", "%s", "%s", "%s", %s);' % (db_name,
                                                                                                                                     locus_tag,
                                                                                                                                     locus_tag,
                                                                                                                                     organism,
                                                                                                                                     taxon_id,
                                                                                                                                     sequence_length,
                                                                                                                                     analysis,
                                                                                                                                     signature_accession,
                                                                                                                                     signature_description,
                                                                                                                                     int(start),
                                                                                                                                     int(stop),
                                                                                                                                     str(score),
                                                                                                                                     interpro_accession,
                                                                                                                                     interpro_description,
                                                                                                                                     GO_terms,
                                                                                                                                     pathways,
                                                                                                                                     orthogroup,
                                                                                                                                     seqfeature_id)
                    try:
                        server.adaptor.execute(sql)
                        server.adaptor.commit()
                    except:
                        print (sql)
                        print (data)
                        import sys
                        sys.exit()


if __name__ == '__main__':
    import argparse
    from chlamdb.biosqldb import manipulate_biosqldb

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_interpro', type=str, help="input interpro csv file", nargs='+')
    parser.add_argument("-d", '--database_name', type=str, help="database name")
    #parser.add_argument("-v2", '--v2_table', action="store_true", help="create V2 table")
    parser.add_argument("-a", '--add_SP_TM', action="store_true", help="Add Ssignal Peptide and TM counts to orthology_detail table")
    parser.add_argument("-u", '--hash2locus_tag', type=str, help="Tab separated file with correspondance between sequence hashes and locus tags")
    parser.add_argument("-l", '--legacy_table', action="store_true", help="Create legacy table in biosqldb")

    args = parser.parse_args()

    biodb = args.database_name
    server, db = manipulate_biosqldb.load_db(biodb)

    if args.hash2locus_tag:
        import chlamdb_setup_utils
        hash2locus_list = chlamdb_setup_utils.get_hash2locus_list(args.hash2locus_tag)

    if args.input_interpro:
        if not args.legacy_table:
            print("creating locus_tag2seqfeature_id")
            sql = 'select locus_tag, seqfeature_id from annotation.seqfeature_id2locus_%s' % biodb
            locus_tag2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
            print("make table v1")
            interpro2biosql(server,
                            hash2locus_list,
                            locus_tag2seqfeature_id,
                            biodb, *args.input_interpro)

        else:
            print("creating locus_tag2seqfeature_id")
            locus_tag2seqfeature_id = manipulate_biosqldb.locus_tag2seqfeature_id_dict(server, biodb)

            print("creating protein_id2seqfeature_id")
            protein_id2seqfeature_id = manipulate_biosqldb.protein_id2seqfeature_id_dict(server, biodb)

            print("getting seqfeature_id2organism")
            seqfeature_id2organism = manipulate_biosqldb.seqfeature_id2organism_dico(server, biodb)

            print("creating locus_tag2taxon_id dictionnary...")
            locus_tag2genome_taxon_id = manipulate_biosqldb.locus_tag2genome_taxon_id(server, biodb)

            print("creating protein_id2taxon_id dictionnary...")
            protein_id2genome_taxon_id = manipulate_biosqldb.protein_id2genome_taxon_id(server, biodb)

            print("getting seqfeature_id2locus_tag")
            seqfeature_id2locus_tag = manipulate_biosqldb.seqfeature_id2locus_tag_dico(server, biodb)

            sql = 'select name,interpro_id from interpro.entry'
            interpro_accession2interpro_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

            '''
            print("make table v2")
            interpro2biosqlV2(server,
                              hash2locus,
                              interpro_accession2interpro_id,
                              seqfeature_id2locus_tag,
                              locus_tag2genome_taxon_id,
                              protein_id2genome_taxon_id,
                              locus_tag2seqfeature_id,
                              protein_id2seqfeature_id,
                              seqfeature_id2organism,
                              biodb, *args.input_interpro)
            '''

            sql = 'select seqfeature_id,orthogroup_name from orthology.seqfeature_id2orthogroup_%s t1 inner join orthology.orthogroup_%s t2 on t1.orthogroup_id=t2.orthogroup_id' % (biodb, biodb)
            seqfeature_id2orthogroup = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

            interpro2biosql_legacy(server,
                                   seqfeature_id2locus_tag,
                                   locus_tag2genome_taxon_id,
                                   protein_id2genome_taxon_id,
                                   locus_tag2seqfeature_id,
                                   protein_id2seqfeature_id,
                                   seqfeature_id2organism,
                                   biodb,
                                   seqfeature_id2orthogroup,
                                   hash2locus_list,
                                   *args.input_interpro)
    if args.add_SP_TM:
        add_TM_and_SP_columns(args.database_name)
