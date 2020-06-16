

from chlamdb.biosqldb import manipulate_biosqldb

def species2taxon_id(biodb):

    server, db = manipulate_biosqldb.load_db(biodb)

    sql2 = 'select distinct family,taxon_id from taxid2species_%s t1 ' \
          ' inner join species_curated_taxonomy_%s t2 on t1.species_id=t2.species_id;' % (biodb, 
                                                                                          biodb)
    
    species_id2taxon_id = {}
    data = server.adaptor.execute_and_fetchall(sql2,)
    for row in data:
        if row[0] not in species_id2taxon_id:
            species_id2taxon_id[row[0]] = [str(row[1])]
        else:
            species_id2taxon_id[row[0]].append(str(row[1]))
    
    return species_id2taxon_id


def get_species_data(server,
                     biodb):
    
    from chlamdb.biosqldb import manipulate_biosqldb

    """
    for each species, report 
    
    - number of complete genomes (1 contig)
    - number of draft genomes (>1 contig)
    - completeness range
 
    - ideally distinguish metagenomes from cultured representatives
    """

    sql1 = f'''
            select t3.species,count(*) as n_complete from biosqldb.bioentry t1
            inner join biosqldb.taxid2species_{biodb} t2 on t1.taxon_id=t2.taxon_id
            inner join biosqldb.species_curated_taxonomy_{biodb} t3 on t2.species_id=t3.species_id
            inner join biosqldb.genomes_info_{biodb} t4 on t1.accession=t4.ACCESSION where t4.n_contigs=1 and t1.description not like "%%%%plasmid%%%%" group by species;
            '''
    sql2 = f'''
            select t3.species,count(*) as n_incomplete from biosqldb.bioentry t1
            inner join biosqldb.taxid2species_{biodb} t2 on t1.taxon_id=t2.taxon_id
            inner join biosqldb.species_curated_taxonomy_{biodb} t3 on t2.species_id=t3.species_id
            inner join biosqldb.genomes_info_{biodb} t4 on t1.accession=t4.ACCESSION where t4.n_contigs>1 and t1.description not like "%%%%plasmid%%%%" group by species;
            '''

    sql3 = f'''
            select t3.species, t5.completeness,t4.n_contigs from biosqldb.bioentry t1
            inner join biosqldb.taxid2species_2019_06_PVC t2 on t1.taxon_id=t2.taxon_id
            inner join biosqldb.species_curated_taxonomy_2019_06_PVC t3 on t2.species_id=t3.species_id
            inner join biosqldb.genomes_info_2019_06_PVC t4 on t1.accession=t4.ACCESSION
            inner join custom_tables.checkm_2019_06_PVC t5 on t1.taxon_id=t5.taxon_id where t1.description not like "%%%%plasmid%%%%";
            '''

    species2n_complete_genomes = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql1,))
    species2n_draft_genomes = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2,))
    
    species2completeness = {}
    
    for row in server.adaptor.execute_and_fetchall(sql3,):
        species, completeness, n_contigs = row 
        completeness = float(completeness)
        if int(n_contigs) == 1:
            continue
        if species not in species2completeness:
            species2completeness[species] = [completeness, completeness]
        else:
            if completeness < species2completeness[species][0]:
                print("smaller!")
                species2completeness[species][0] = completeness
            if completeness > species2completeness[species][1]:
                print("Larger")
                species2completeness[species][1] = completeness
    print(species2completeness)
    return species2n_complete_genomes, species2n_draft_genomes, species2completeness


def get_species_tree(biodb):
    
    from ete3 import Tree,TreeStyle
    
    server, db = manipulate_biosqldb.load_db(biodb)
    
    species2n_complete_genomes, species2n_draft_genomes, species2completeness = get_species_data(server,
                                                                                                 biodb)

    
    sql_tree = 'select tree from reference_phylogeny t1 inner join biodatabase t2 on t1.biodatabase_id=t2.biodatabase_id ' \
               ' where t2.name="%s";' % biodb
               
    server, db = manipulate_biosqldb.load_db(biodb)
    complete_tree = Tree(server.adaptor.execute_and_fetchall(sql_tree,)[0][0])
    R = complete_tree.get_midpoint_outgroup()
    complete_tree.set_outgroup(R)

    sql = 'select distinct taxon_id,species from taxid2species_%s t1 ' \
          ' inner join species_curated_taxonomy_%s t2 on t1.species_id=t2.species_id;' % (biodb, 
                                                                                          biodb)
          
    taxon_id2species_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
    
    
    # changing taxon id to species id
    for leaf in complete_tree.iter_leaves():
        #print '%s --> %s' % (leaf.name, str(taxon_id2species_id[str(leaf.name)]))
        leaf.name = "%s" % str(taxon_id2species_id[str(leaf.name)])

    # attributing unique id to each node
    # if all node descendant have the same name, use that name as node name
    n = 0
    for node in complete_tree.traverse():
        if node.name=='':
            desc_list = list(set([i.name for i in node.iter_descendants()]))
            try:
                desc_list.remove('')
            except ValueError:
                pass
            if len(desc_list) != 1:
                node.name = '%sbb' % n
            else:
                node.name = desc_list[0]
            n+=1
 
    # Collapsing nodes while traversing
    # http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#collapsing-nodes-while-traversing-custom-is-leaf-definition
    node2labels = complete_tree.get_cached_content(store_attr="name")
    
    def collapsed_leaf(node):
        if len(node2labels[node]) == 1:
            return True
        else:
            return False

    species_tree = Tree(complete_tree.write(is_leaf_fn=collapsed_leaf))
    
    
    for lf_count, lf in enumerate(species_tree.iter_leaves()):
        
        try:
            n_complete_genomes = species2n_complete_genomes[lf.name]
        except:
            n_complete_genomes = False
        try:
            n_draft_genomes = species2n_draft_genomes[lf.name]
        except:
            n_draft_genomes = False   

        if n_draft_genomes:
            c1 = round(species2completeness[lf.name][0])
            c2 = round(species2completeness[lf.name][1])
            if c1 == c2:
                completeness = "%s%%" % c1
            else:
                completeness = "%s-%s%%" % (c1, c2)
        if n_complete_genomes and n_draft_genomes:

            lf.name = "%s (%sc/%sd, %s)" % (lf.name,
                                       n_complete_genomes,
                                       n_draft_genomes,
                                       completeness)

        if n_complete_genomes and not n_draft_genomes:
            lf.name = "%s (%sc)" % (lf.name,
                                    n_complete_genomes)
        if not n_complete_genomes and n_draft_genomes:
            lf.name = "%s (%sd, %s)" % (lf.name,
                                    n_draft_genomes,
                                    completeness)
    
    return complete_tree, species_tree