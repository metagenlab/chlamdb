
def get_signature_taxonomy(signature_id, 
                        biodb,
                        rank='phylum'):
    import MySQLdb
    import os
    from ete3 import NCBITaxa, Tree, TextFace,TreeStyle, StackedBarFace
    ncbi = NCBITaxa()
    from chlamdb.biosqldb import manipulate_biosqldb
    
    server, db = manipulate_biosqldb.load_db(biodb)
    

    sql_domain_taxonomy = f"""select species_taxid from refseq_ref_repres_genomes_interpro_annots t1 
                              inner join refseq_ref_repres_genomes t2 on t1.assembly_id=t2.assembly_id 
                              where signature_id={signature_id} group by species_taxid;"""
    taxid_with_domain_list = [i[0] for i in server.adaptor.execute_and_fetchall(sql_domain_taxonomy,)]


    sql = f'select phylogeny from refseq_ref_repres_genomes_phylogeny where dataset="interpro" and tax_rank="{rank}"'
    tree_string = server.adaptor.execute_and_fetchall(sql)[0][0]

    tree = Tree(tree_string)

    leaves_list = [i.name for i in tree.iter_leaves()]

    taxon_id2lineage = {}
    for taxon in taxid_with_domain_list:
        taxon_id2lineage[taxon] = ncbi.get_lineage(taxon)

    leaf_taxon2n_species_with_domain = {}
    for leaf_taxon in leaves_list:
        leaf_taxon2n_species_with_domain[leaf_taxon] = 0
        for taxon in taxid_with_domain_list:
            lineage = taxon_id2lineage[taxon]
            if int(leaf_taxon) in lineage:
                leaf_taxon2n_species_with_domain[leaf_taxon]+=1

                #if taxon in taxid_with_domain_list:
                #    leaf_taxon2n_species_with_domain[leaf_taxon]+=1
    return leaf_taxon2n_species_with_domain


def get_interpro_taxonomy(interpro_id, 
                        biodb,
                        rank='phylum'):
    import MySQLdb
    import os
    from ete3 import NCBITaxa, Tree, TextFace,TreeStyle, StackedBarFace
    ncbi = NCBITaxa()
    from chlamdb.biosqldb import manipulate_biosqldb
    
    server, db = manipulate_biosqldb.load_db(biodb)
    

    sql_domain_taxonomy = f"""select distinct species_taxid from refseq_ref_repres_genomes_interpro_annots t1 
    inner join refseq_ref_repres_genomes t2 on t1.assembly_id=t2.assembly_id 
    inner join refseq_ref_repres_genomes_interpro_entries t3 on t1.signature_id=t3.signature_id 
    inner join interpro_entry t4 on t3.interpro_id=t4.interpro_id where t4.interpro_id={interpro_id};
    """
    taxid_with_domain_list = [i[0] for i in server.adaptor.execute_and_fetchall(sql_domain_taxonomy,)]


    sql = f'select phylogeny from refseq_ref_repres_genomes_phylogeny where dataset="interpro" and tax_rank="{rank}"'
    tree_string = server.adaptor.execute_and_fetchall(sql)[0][0]

    tree = Tree(tree_string)

    leaves_list = [i.name for i in tree.iter_leaves()]

    taxon_id2lineage = {}
    for taxon in taxid_with_domain_list:
        taxon_id2lineage[taxon] = ncbi.get_lineage(taxon)

    leaf_taxon2n_species_with_domain = {}
    for leaf_taxon in leaves_list:
        leaf_taxon2n_species_with_domain[leaf_taxon] = 0
        for taxon in taxid_with_domain_list:
            lineage = taxon_id2lineage[taxon]
            if int(leaf_taxon) in lineage:
                leaf_taxon2n_species_with_domain[leaf_taxon]+=1

                #if taxon in taxid_with_domain_list:
                #    leaf_taxon2n_species_with_domain[leaf_taxon]+=1
    return leaf_taxon2n_species_with_domain



def plot_phylum_counts(domain_id, 
                       biodb,
                       rank='phylum',
                       colapse_low_species_counts=4,
                       remove_unlassified=True):

    '''

    1. get phylum tree
    2. foreach species => get phylum
    3. build phylum2count dictionnary
    3. plot barchart

    # merge eukaryotes into 5 main clades
    # merge virus as a single clade


    ATTENTION: no-rank groups and no-rank species...

    '''

    import MySQLdb
    import os
    from chlamdb.biosqldb import manipulate_biosqldb
    from ete3 import NCBITaxa, Tree, TextFace,TreeStyle, StackedBarFace
    ncbi = NCBITaxa()

    from chlamdb.biosqldb import manipulate_biosqldb
    
    server, db = manipulate_biosqldb.load_db(biodb)

    
    if domain_id.startswith('IPR'):
        sql = f'''select t2.interpro_id, t2.name, t2.description from refseq_ref_repres_genomes_interpro_entries t1 
        inner join interpro_entry t2 on t1.interpro_id=t2.interpro_id
        where t2.name ="{domain_id}"
        '''
    else:
        sql = f'select signature_id, name, type from refseq_ref_repres_genomes_interpro_entries where accession like "{domain_id}%%%%";'
    print(sql)
    signature_id, signature_name, signature_type = server.adaptor.execute_and_fetchall(sql,)[0]

    sql = f'select taxon_id, n_genomes from refseq_ref_repres_genomes_leaf2n_genomes where tax_rank="{rank}" and dataset="interpro"'
    leaf_taxon2n_species = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql))

    if domain_id.startswith('IPR'):
        leaf_taxon2n_species_with_domain = get_interpro_taxonomy(signature_id, biodb, rank)
    else:

        leaf_taxon2n_species_with_domain = get_signature_taxonomy(signature_id, biodb, rank)

    sql = f'select phylogeny from refseq_ref_repres_genomes_phylogeny where dataset="interpro" and tax_rank="{rank}"'

    tree = Tree(server.adaptor.execute_and_fetchall(sql)[0][0], format=1)

    sql = f'select taxon_id,scientific_name,node_rank from refseq_ref_repres_genomes_taxid2label where tax_rank="{rank}" and dataset="interpro"'

    taxon_id2scientific_name_and_rank = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql))
    taxon_id2scientific_name_and_rank = {str(k):v for k,v in taxon_id2scientific_name_and_rank.items()}

    print(taxon_id2scientific_name_and_rank)

    tss = TreeStyle()
    tss.draw_guiding_lines = True
    tss.guiding_lines_color = "blue"

    keep = []
    for lf in tree.iter_leaves():
        # n genomes

        if remove_unlassified:
            label = taxon_id2scientific_name_and_rank[str(lf.name)][0]
            if 'unclassified' in label:
                continue

        n_genomes = int(leaf_taxon2n_species[lf.name])
        if n_genomes > colapse_low_species_counts:
            keep.append(lf.name)
    print ('number of leaves:', len(keep))

    tree.prune(keep)

    header_list = ['Rank', 'N genomes', 'N with %s' % domain_id, 'Percentage']
    for col, header in enumerate(header_list):

        n = TextFace('%s' % (header))
        n.margin_top = 0
        n.margin_right = 1
        n.margin_left = 20
        n.margin_bottom = 1
        n.rotation = 270
        n.hz_align = 2
        n.vt_align = 2
        n.inner_background.color = "white"
        n.opacity = 1.
        tss.aligned_header.add_face(n, col)

    for lf in tree.iter_leaves():
        # n genomes

        n_genomes = int(leaf_taxon2n_species[lf.name])
        if n_genomes <= colapse_low_species_counts:
            continue


        n = TextFace('  %s ' % str(leaf_taxon2n_species[lf.name]))
        n.margin_top = 1
        n.margin_right = 1
        n.margin_left = 0
        n.margin_bottom = 1
        n.fsize = 7
        n.inner_background.color = "white"
        n.opacity = 1.
        lf.add_face(n, 2, position="aligned")

        # n genomes with domain
        m = TextFace('  %s ' % str(leaf_taxon2n_species_with_domain[lf.name]))
        m.margin_top = 1
        m.margin_right = 1
        m.margin_left = 0
        m.margin_bottom = 1
        m.fsize = 7
        m.inner_background.color = "white"
        m.opacity = 1.
        lf.add_face(m, 3, position="aligned")

        # rank
        ranks = ncbi.get_rank([lf.name])
        try:
            r = ranks[max(ranks.keys())]
        except:
            r = '-'
        n = TextFace('  %s ' % r, fsize=14, fgcolor='red')
        n.margin_top = 1
        n.margin_right = 1
        n.margin_left = 0
        n.margin_bottom = 1
        n.fsize = 7
        n.inner_background.color = "white"
        n.opacity = 1.
        lf.add_face(n, 1, position="aligned")

        # percent with target domain
        percentage = (float(leaf_taxon2n_species_with_domain[lf.name])/float(leaf_taxon2n_species[lf.name]))*100

        m = TextFace('  %s ' % str(round(percentage,2)))
        m.fsize = 1
        m.margin_top = 1
        m.margin_right = 1
        m.margin_left = 0
        m.margin_bottom = 1
        m.fsize = 7
        m.inner_background.color = "white"
        m.opacity = 1.
        lf.add_face(m, 4, position="aligned")


        b = StackedBarFace([percentage,
                            100-percentage],
                            width=100, height=10, colors=["#7fc97f", "white"])
        b.rotation= 0
        b.inner_border.color = "grey"
        b.inner_border.width = 0
        b.margin_right = 15
        b.margin_left = 0
        lf.add_face(b, 5, position="aligned")

        n = TextFace('%s' % taxon_id2scientific_name_and_rank[str(lf.name)][0], fgcolor = "black", fsize = 9) # , fstyle = 'italic'

        lf.name = " %s (%s)" % (taxon_id2scientific_name_and_rank[str(lf.name)][0], str(lf.name))
        n.margin_right = 10
        lf.add_face(n, 0)

    tss.show_leaf_name = False

    for node in tree.traverse("postorder"):
        try:
            r = taxon_id2scientific_name_and_rank[str(node.name)][1]
        except:
            pass
        try:
            if r in ['phylum', 'superkingdom', 'class', 'subphylum'] or taxon_id2scientific_name_and_rank[str(node.name)][0] in ['FCB group']:

                hola = TextFace("%s" % (taxon_id2scientific_name_and_rank[str(node.name)][0]))
                node.add_face(hola, column=0, position = "branch-top")
        except:
            pass
    return tree, tss, signature_name, signature_type
    #print tree.get_ascii(attributes=["name", "rank"])