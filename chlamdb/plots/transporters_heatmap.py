#!/usr/bin/python

def transporter_all_superfamily_heatmap(biodb, evalue, bitscore, query_cov, hit_cov, total=False):
    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select t1.taxon_id,t3.description, count(*) as n from transporters_transporters t1 ' \
          ' inner join transporters.transporter_table t2 on t1.transporter_id=t2.transporter_id ' \
          ' inner join transporters.tc_table t3 on t2.superfamily=t3.tc_id ' \
          ' where query_cov>=%s and hit_cov>=%s and evalue<=%s and bitscore_first_hsp>=%s' \
          ' group by taxon_id,tc_id;' % (biodb, query_cov, hit_cov, evalue, bitscore)

    data = server.adaptor.execute_and_fetchall(sql,)

    code2taxon2count = {}
    if total:
        code2taxon2count['TOTAL'] = {}
        transporter_list = ['TOTAL']
    else:
        transporter_list = []
    for row in data:
        row = list(row)
        if total:
            if str(row[0]) not in code2taxon2count['TOTAL']:
                code2taxon2count['TOTAL'][str(row[0])] = int(row[2])
            else:
                code2taxon2count['TOTAL'][str(row[0])] += int(row[2])

        if row[1] not in transporter_list:
            transporter_list.append(row[1])
        if row[1] not in code2taxon2count:
            code2taxon2count[row[1]] = {}
            code2taxon2count[row[1]][str(row[0])] = int(row[2])
        else:
            code2taxon2count[row[1]][str(row[0])] = int(row[2])
    return transporter_list, code2taxon2count

def transporter_fam_heatmap(biodb, family_list, rank, evalue, bitscore, query_cov, hit_cov, total=False):
    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)

    tc_name_filter = '","'.join(family_list)

    if rank == 'family':
        rank_filter = f'and t4.tc_name in ("{tc_name_filter}")'
        col = 't4.tc_name'
    elif rank == 'subfamily':
        rank_filter = f'and t5.tc_name in ("{tc_name_filter}")'
        col = 't5.tc_name'
    elif rank == 'superfamily':
        rank_filter = f'and t6.tc_name in ("{tc_name_filter}")'
        col = 't6.tc_name'
    elif rank == 'tc_name':
        rank_filter = f'and t3.tc_name in ("{tc_name_filter}")'
        col = 't3.tc_name'

    sql = 'select t1.taxon_id,%s, count(*) as n from transporters_transporters_BBH t1 ' \
          ' inner join transporters_protein_entry2transporter t2 on t1.hit_protein_id=t2.protein_id ' \
          ' inner join transporters_transporters t3 on t2.transporter_id=t3.transporter_id ' \
          ' inner join transporters_classification t4 on t3.family=t4.tc_id ' \
          ' inner join transporters_classification t5 on t3.family=t5.tc_id ' \
          ' inner join transporters_classification t6 on t3.superfamily=t6.tc_id ' \
          ' where query_cov>=%s and hit_cov>=%s and evalue<=%s and bitscore_first_hsp>=%s %s ' \
          ' group by taxon_id,%s;' % (col, query_cov, hit_cov, evalue, bitscore, rank_filter, col)
    print (sql)

    
    data = server.adaptor.execute_and_fetchall(sql,)

    code2taxon2count = {}
    if total:
        code2taxon2count['TOTAL'] = {}
        transporter_list = ['TOTAL']
    else:
        transporter_list = []
    for row in data:
        row = list(row)
        if total:
            if str(row[0]) not in code2taxon2count['TOTAL']:
                code2taxon2count['TOTAL'][str(row[0])] = int(row[2])
            else:
                code2taxon2count['TOTAL'][str(row[0])] += int(row[2])

        if row[1] not in transporter_list:
            transporter_list.append(row[1])

        if row[1] not in code2taxon2count:
            code2taxon2count[row[1]] = {}

            code2taxon2count[row[1]][str(row[0])] = int(row[2])
        else:
            code2taxon2count[row[1]][str(row[0])] = int(row[2])

    return transporter_list, code2taxon2count


def transporter_superfamily_heatmap(biodb, family,evalue, bitscore, query_cov, hit_cov):
    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select t1.taxon_id,t3.description, count(*) as n, t3.tc_name from transporters_transporters t1 ' \
          ' inner join transporters.transporter_table t2 on t1.transporter_id=t2.transporter_id ' \
          ' inner join transporters.tc_table t3 on t2.family=t3.tc_id ' \
          ' inner join transporters.tc_table t4 on t2.superfamily=t4.tc_id ' \
          ' where query_cov>%s and hit_cov>%s and evalue<=%s and bitscore_first_hsp>=%s and t4.description="%s" ' \
          ' group by taxon_id,t3.tc_id;' % (biodb, query_cov, hit_cov, evalue, bitscore, family)
    print (sql)
    data = server.adaptor.execute_and_fetchall(sql,)

    code2taxon2count = {}
    code2taxon2count['TOTAL'] = {}
    transporter_list = ['TOTAL']
    for row in data:
        if row[3] not in row[1]:
            family_label = "%s: %s" % (row[3], row[1])
        else:
            family_label = row[1]
        row = list(row)
        if str(row[0]) not in code2taxon2count['TOTAL']:
            code2taxon2count['TOTAL'][str(row[0])] = int(row[2])
        else:
            code2taxon2count['TOTAL'][str(row[0])] += int(row[2])

        if family_label not in transporter_list:
            transporter_list.append(family_label)

        if family_label not in code2taxon2count:
            code2taxon2count[family_label] = {}

            code2taxon2count[family_label][str(row[0])] = int(row[2])
        else:
            code2taxon2count[family_label][str(row[0])] = int(row[2])

    return transporter_list, code2taxon2count
