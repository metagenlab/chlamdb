

def perform_search(request, 
                   biodb,
                   search_term,
                   redo=True):

    from chlamdb.biosqldb import manipulate_biosqldb 
    from chlamdb.views import fam, KEGG_module_map, KEGG_mapp_ko, locusx
    import re
    from django.shortcuts import render

    server, db = manipulate_biosqldb.load_db(biodb)

    # check if single or multiple search terms

    split_term = search_term.split(" ")

    if len(split_term) != 1:
        synonymous = False
    else:
        print("check exact acc")
        print(re.match("^[0-9\.]+$", search_term))
        if len(search_term) == len("PF04093") and search_term[0:2] == 'PF':
            return fam(request, search_term, 'pfam')
        elif len(search_term) == len('K03652') and search_term[0:1] == 'K':
            return fam(request, search_term, 'ko')
        elif len(search_term) == len('COG0001') and search_term[0:3] == 'COG':
            return fam(request, search_term, 'cog')
        elif len(search_term) == len('IPR000014') and search_term[0:3] == 'IPR':
            return fam(request, search_term, 'interpro')
        elif len(search_term) == len('M00406') and search_term[0:3] == 'M00':
            return KEGG_module_map(request, search_term)
        elif len(search_term) == len('map00550') and search_term[0:3] == 'map':
            return KEGG_mapp_ko(request, search_term)
        elif re.match("^[0-9\.]+$", search_term) is not None and len(search_term.split(".")) == 4:
            return fam(request, search_term, 'EC')
        else:
        
            # try to match synonymous table
            sql_sy = f'select t1.db_name,t1.accession,t2.orthogroup,t2.locus_tag,t2.start,t2.stop,t2.strand,t2.gene, t2.orthogroup_size,t2.n_genomes,t2.TM,t2.SP,t2.product,t2.organism ' \
                    f' from biosqldb.cross_references_{biodb} t1 ' \
                    f' inner join orthology_detail_{biodb} t2 on t1.seqfeature_id=t2.seqfeature_id ' \
                    f' where t1.accession="{search_term}" ' \
                    f' and db_name not in ("STRING", "KEGG", "KO", "eggNOG") group by t1.seqfeature_id limit 200;'

            raw_search = server.adaptor.execute_and_fetchall(sql_sy)
            if len(raw_search) == 0:
                synonymous = False
            elif len(raw_search) == 1:
                locus = raw_search[0][3]
                # return the locus page
                return [locus]
            else:
                synonymous = True
                n = 1
                locus_list = []
                for one_hit in raw_search:
                    
                    locus_list.append((n,) + one_hit)
                    n+=1
                print(locus_list[0:10])
                #return render(request, 'chlamdb/search_synonymous.html', locals())

    if not synonymous:

        # CREATE FULLTEXT INDEX GPF1 ON orthology_detail_2019_06_PVC(gene);
        # CREATE FULLTEXT INDEX GPF2 ON orthology_detail_2019_06_PVC(product);
        # CREATE FULLTEXT INDEX GPF3 ON orthology_detail_2019_06_PVC(organism);
        # CREATE FULLTEXT INDEX GPF4 ON orthology_detail_2019_06_PVC(gene,product,organism);
        columns = 'orthogroup, locus_tag, protein_id, start, stop, ' \
                            'strand, gene, orthogroup_size, n_genomes, TM, SP, product, organism, translation'
        sql = 'SELECT %s, ' \
                ' MATCH (gene) AGAINST ("%s") AS rel1, ' \
                ' MATCH (product) AGAINST ("%s") AS rel2, ' \
                ' MATCH (organism) AGAINST ("%s") AS rel3 FROM orthology_detail_%s t1 ' \
                ' WHERE MATCH (gene,product,organism) AGAINST ("%s" IN BOOLEAN MODE) ORDER BY (rel1)+(rel2*0.5)+(rel3*10) DESC limit 100;' % (columns,
                                                                                                                            search_term,
                                                                                                                            search_term,
                                                                                                                            search_term,
                                                                                                                            biodb,
                                                                                                                            search_term)

        raw_data_gene_raw_data_product = server.adaptor.execute_and_fetchall(sql,)

        n = 1
        locus_list = []
        for one_hit in raw_data_gene_raw_data_product:
            if one_hit[2] != '-':
                interpro_id = one_hit[2]
            else:
                interpro_id = one_hit[1]
            locus_list.append((n,) + one_hit + (interpro_id,))
            n+=1


    # CREATE FULLTEXT INDEX ezf ON enzyme.enzymes_dat(value);
    sql = 'select A.ec, A.value from (select ec,value from enzyme.enzymes_dat as t1 inner join enzyme.enzymes as t2 ' \
            ' on t1.enzyme_dat_id=enzyme_id WHERE MATCH(value) AGAINST("%s" IN NATURAL LANGUAGE MODE) group by ec) A inner join comparative_tables.EC_%s' \
            ' as B on A.ec=B.id' % (search_term, biodb)

    raw_data_EC = server.adaptor.execute_and_fetchall(sql,)
    if len(raw_data_EC) == 0:
        raw_data_EC = False

    # CREATE FULLTEXT INDEX koaf ON enzyme.ko_annotation_v1(definition);
    # filter KO absent from DB
    # CREATE FULLTEXT INDEX ko_full ON ko_annotation_v1(ko_id,name,definition,EC);
    sql = 'select A.ko_id,A.name,A.definition from (select ko_id,name,definition from enzyme.ko_annotation_v1 ' \
            'WHERE MATCH(ko_id,name,definition,EC) AGAINST("%s" IN NATURAL LANGUAGE MODE)) A inner join comparative_tables.ko_%s as B on A.ko_id=B.id' % (search_term, biodb)
    try:
        raw_data_ko = server.adaptor.execute_and_fetchall(sql,)
    except:
        raw_data_ko = []
    if len(raw_data_ko) == 0:
        raw_data_ko = False

    # CREATE FULLTEXT INDEX cogf ON COG.cog_names_2014(description);
    sql = ' select COG_name,code,t3.description,t1.description from COG.cog_names_2014 t1 ' \
            ' inner join COG.cog_id2cog_category t2 on t1.COG_id=t2.COG_id ' \
            ' inner join COG.code2category t3 on t2.category_id=t3.category_id ' \
            ' WHERE MATCH(t1.description) AGAINST("%s" IN NATURAL LANGUAGE MODE);' % (search_term)
    raw_data_cog = []# :TODO server.adaptor.execute_and_fetchall(sql,)
    if len(raw_data_cog) == 0:
        raw_data_cog = False

    # CREATE FULLTEXT INDEX ipf ON interpro.signature(signature_description);
    # CREATE FULLTEXT INDEX ipf ON interpro.entry(description);
    sql = 'select analysis_name,signature_accession,signature_description,t3.name,t3.description from interpro.signature t1 ' \
            ' inner join interpro.analysis t2 on t1.analysis_id=t2.analysis_id ' \
            ' left join interpro.entry t3 on t1.interpro_id=t3.interpro_id ' \
            ' WHERE MATCH(t3.description) AGAINST("%s" IN NATURAL LANGUAGE MODE) limit 100;' % (search_term)
    raw_data_interpro = server.adaptor.execute_and_fetchall(sql,)
    if len(raw_data_interpro) == 0:
        raw_data_interpro = False

    # CREATE FULLTEXT INDEX modf ON enzyme.kegg_module_v1(description, module_sub_cat);
    sql = 'select module_name,module_sub_cat,module_sub_sub_cat,description from enzyme.kegg_module_v1 ' \
            ' where MATCH(description,module_sub_cat) AGAINST("%s" IN BOOLEAN MODE)' % (search_term)
    raw_data_module = server.adaptor.execute_and_fetchall(sql,)
    if len(raw_data_module) == 0:
        raw_data_module = False

    # CREATE FULLTEXT INDEX kegf ON enzyme.kegg_pathway(description, pathway_category);
    sql = 'select pathway_name,pathway_category,description from enzyme.kegg_pathway ' \
            ' where MATCH(description, pathway_category) AGAINST("%s" IN BOOLEAN MODE) ' % (search_term)
    raw_data_pathway = server.adaptor.execute_and_fetchall(sql,)
    if len(raw_data_pathway) == 0:
        raw_data_pathway = False
    
    # ALTER TABLE pmid2data_stringdb ADD FULLTEXT (title,journal,year,authors);;
    sql1 = 'select * from (SELECT pmid,title, MATCH (title) AGAINST ("%s" IN BOOLEAN MODE) as score ' \
            ' FROM string.pmid2data_paperblast )A where score > 0 order by score DESC limit 20;' % (search_term)
    sql2 = 'SELECT pmid,title,year,journal,authors from string.pmid2data_stringdb where MATCH(title,journal,year,authors)' \
            ' AGAINST ("%s" IN BOOLEAN MODE) > 0 limit 50;' % (search_term)
    #raw_data_pmid_paperblast2title = server.adaptor.execute_and_fetchall(sql1,)
    raw_data_pmid = list(server.adaptor.execute_and_fetchall(sql2,))
    for n,row in enumerate(raw_data_pmid):
        raw_data_pmid[n] = [str(i) for i in row]
    pmid2title = {}
    #pmid2title.update(raw_data_pmid_paperblast2title)
    #pmid2title.update(raw_data_pmid_string2title)

    pmid_data = False
    if len(raw_data_pmid) == 0:
        pmid_data = False
        pmid2n_homologs = False
    else:

        pmid_list = [i[0] for i in raw_data_pmid]

        sql = 'select A.pmid,count(*) as n from (select t1.pmid,orthogroup_id from string.pmid2data_stringdb t1 ' \
                ' inner join string.string_protein2pmid t2 on t1.pmid=t2.pmid ' \
                ' inner join string.seqfeature_id2string_protein_mapping t3 on t2.string_protein_id=t3.string_protein_id ' \
                ' inner join orthology.seqfeature_id2orthogroup_%s t4 on t3.seqfeature_id=t4.seqfeature_id ' \
                ' where t1.pmid in (%s) group by t1.pmid,orthogroup_id) A group by pmid;' % (biodb, ','.join(pmid_list))

        pmid2n_homologs = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

        pmid_data = True
        '''
        set_pmid_paperblast = set(raw_data_pmid_paperblast2title.keys())
        set_pmid_string = set(raw_data_pmid_string2title.keys())
        pmid2presence = {}
        for pmid in pmid2title:
            if pmid in set_pmid_paperblast and pmid in set_pmid_string:
                pmid2presence[pmid] = [1,1]
            if pmid in set_pmid_paperblast and pmid not in set_pmid_string: 
                pmid2presence[pmid] = [1,0]
            if pmid not in set_pmid_paperblast and pmid in set_pmid_string: 
                pmid2presence[pmid] = [0,1]    
        '''
                    
    if not raw_data_module \
            and not raw_data_pathway \
            and not raw_data_EC \
            and not raw_data_ko \
            and not raw_data_cog \
            and not raw_data_interpro \
            and not locus_list \
            and not pmid_data:
        return None 
    else:
        return [locus_list, 
                raw_data_EC,
                raw_data_ko,
                raw_data_cog,
                raw_data_interpro,
                raw_data_pathway,
                raw_data_pmid,
                pmid2n_homologs,
                pmid_data,
                synonymous]



    '''
    else:
        n = 1
        search_result = []
        for one_hit in raw_data:
            if one_hit[2] != '-':
                interpro_id = one_hit[2]
            else:
                interpro_id = one_hit[1]
            search_result.append((n,) + one_hit + (interpro_id,))
            n+=1

        return search_result
    '''



