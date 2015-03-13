#!/usr/bin/python

import heatmap
import numpy as np
import manipulate_biosqldb
import parse_newick_tree
from Bio import Phylo

# heatmap Chlamydiales pan-genome

def biodb2heatmap(bio_db_name):

    server, db = manipulate_biosqldb.load_db(bio_db_name)
    matrix = np.array(manipulate_biosqldb.get_orthology_table(server, bio_db_name))

    taxon_id2description = manipulate_biosqldb.taxon_id2genome_description(server, bio_db_name)
    #print taxon_id2description

    group_names = matrix[:,0]
    taxons_ids = manipulate_biosqldb.get_taxon_id_list(server, bio_db_name)
    print taxons_ids
    print len(taxons_ids)
    taxons_ids = [taxon_id2description[str(i)] for i in taxons_ids]


    M = matrix.astype(float) # [:, 1:]
    M = heatmap.randomize_table(M)
    print len(M[:,1])
    print len(M[1,:])
    heatmap.heatmap(M, output="heatmap.pdf", breaks="-0.5, 0.5, 1.5, 2.5", rows=None, columns=taxons_ids, format="pdf", orderCols=True)


def write_ortho_matrix(bio_db_name):
    server, db = manipulate_biosqldb.load_db(bio_db_name)
    matrix = np.array(manipulate_biosqldb.get_orthology_table(server, bio_db_name))

    taxon_id2description = manipulate_biosqldb.taxon_id2genome_description(server, bio_db_name)


    group_names = matrix[:,0]
    taxons_ids = manipulate_biosqldb.get_taxon_id_list(server, bio_db_name)

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
        taxon_id2description[i] = re.sub(" DNA", "", taxon_id2description[i])

    taxons_ids = [taxon_id2description[str(i)] for i in taxons_ids]
    #print taxon_id2description
    f = open("ortho_matrix.tab", "w")

    f.write("orthogroup\t" + "\t".join(taxons_ids) + "\n")

    for row in matrix:
        f.write("\t".join(row) + "\n")
    f.close()




#accession2description = manipulate_biosqldb.accession2description(server, "saureus1")

"""
new_tree = parse_newick_tree.convert_terminal_node_names("tree.ML.tre", accession2description)
print new_tree[0]

Phylo.write(new_tree, 'renamed_tree.ML.tre', 'newick')
"""




def convert_tree_taxon2genome(biodb_name, input_tree, output_tree ):
    server, db = manipulate_biosqldb.load_db(biodb_name)
    taxon_id2genome_description = manipulate_biosqldb.taxon_id2genome_description(server, biodb_name)

    locus2genome = manipulate_biosqldb.locus_tag2genome_name(server, biodb_name)


    import re
    for i in taxon_id2genome_description.keys():
        taxon_id2genome_description[i] = re.sub(" subsp\. aureus", "", taxon_id2genome_description[i])
        taxon_id2genome_description[i] = re.sub(", complete genome\.", "", taxon_id2genome_description[i])
        taxon_id2genome_description[i] = re.sub(", complete sequence\.", "", taxon_id2genome_description[i])
        taxon_id2genome_description[i] = re.sub("strain ", "", taxon_id2genome_description[i])
        taxon_id2genome_description[i] = re.sub("str\. ", "", taxon_id2genome_description[i])
        taxon_id2genome_description[i] = re.sub(" complete genome sequence\.", "", taxon_id2genome_description[i])
        taxon_id2genome_description[i] = re.sub(" complete genome\.", "", taxon_id2genome_description[i])
        taxon_id2genome_description[i] = re.sub(" chromosome", "", taxon_id2genome_description[i])
        taxon_id2genome_description[i] = re.sub("Staphylococcus", "S.", taxon_id2genome_description[i])
        taxon_id2genome_description[i] = re.sub(" DNA", "S.", taxon_id2genome_description[i])
    #print taxon_id2genome_description[i]



    #print taxon_id2genome_description
    new_tree = parse_newick_tree.convert_terminal_node_names(input_tree, taxon_id2genome_description)
    #print new_tree[0]
    print "writing converted tree..."
    print output_tree
    Phylo.write(new_tree, output_tree, 'newick')

def convert_tree_accession2genome(biodb_name, input_tree, output_tree ):
    server, db = manipulate_biosqldb.load_db(biodb_name)
    #taxon_id2genome_description = manipulate_biosqldb.taxon_id2genome_description(server, biodb_name)
    accession2genome_description = manipulate_biosqldb.accession2description(server, biodb_name)
    #accession2genome_description = taxon_id2genome_description

    locus2genome = manipulate_biosqldb.locus_tag2genome_name(server, biodb_name)


    import re
    for i in accession2genome_description.keys():
        accession2genome_description[i] = re.sub(" subsp\. aureus", "", accession2genome_description[i])
        accession2genome_description[i] = re.sub(", complete genome\.", "", accession2genome_description[i])
        accession2genome_description[i] = re.sub(", complete sequence\.", "", accession2genome_description[i])
        accession2genome_description[i] = re.sub("strain ", "", accession2genome_description[i])
        accession2genome_description[i] = re.sub("str\. ", "", accession2genome_description[i])
        accession2genome_description[i] = re.sub(" complete genome sequence\.", "", accession2genome_description[i])
        accession2genome_description[i] = re.sub(" complete genome\.", "", accession2genome_description[i])
        accession2genome_description[i] = re.sub(" chromosome", "", accession2genome_description[i])
        accession2genome_description[i] = re.sub("Staphylococcus", "S.", accession2genome_description[i])
        accession2genome_description[i] = re.sub(" DNA", "S.", accession2genome_description[i])
    #print accession2genome_description[i]



    #print taxon_id2genome_description
    new_tree = parse_newick_tree.convert_terminal_node_names(input_tree, accession2genome_description)
    #print new_tree[0]

    Phylo.write(new_tree, output_tree, 'newick')



########### MAIN ####################
#####################################

if __name__ == '__main__':
    convert_tree("saureus_01_15",
             "/home/trestan/Dropbox/projets/saureus/results/MiSeq/final/trees/saureus_MLST_BioNJ_tree.phyloxml",
             '/home/trestan/Dropbox/projets/saureus/results/MiSeq/final/trees/saureus_MLST_BioNJ_tree_renamed.tree')

    #biodb2heatmap("Chlamydia_01_15")
    #write_ortho_matrix("giant_virus_12_14")


    """
    #print locus2genome.keys()

    import glob
    from Bio import Phylo
    listing = glob.glob("/home/trestan//Dropbox/dev/django/test_1/assets/Chlamydia_11_14_fasta/*phy_reroot.txt")
    from ete2 import Tree, TreeStyle
    #print "listing", listing



    import re
    import StringIO
    from ete2 import Phyloxml
    for tree in listing:
        out_name = re.sub("phy_reroot", "reroot_renamed", tree)
        out_map = re.sub("phy_reroot", "reroot_map", tree)
        out_tree1 = re.sub("phy_reroot", "reroot_python", tree)
        out_tree = re.sub("txt", "svg", out_tree1)

        out_reroot = re.sub("phy_reroot", "_reroot", tree)

        print out_map
        new_tree, map_file = parse_newick_tree.combine_terminal_node_names(tree, locus2genome)

        f = open(out_map, "w")
        f.write(map_file)
        Phylo.write(new_tree, out_name, "newick")



        str_output = StringIO.StringIO()
        Phylo.write(new_tree, "temp.phyloxml", 'phyloxml')
        p = Phyloxml()
        p.build_from_file("temp.phyloxml")
        t = p.get_phylogeny()[0]
        ts = TreeStyle()
        ts.show_leaf_name = True
        ts.mode = "c"
        ts.arc_start = -90 # 0 degrees = 3 o'clock
        ts.arc_span = 180
        t.render(out_tree, w=200, units="mm", tree_style=ts)
        R = t.get_midpoint_outgroup()
        print R
        t.set_outgroup(R)

        t.write(format=1, outfile=out_reroot)

        f = open(out_map, "w")
        f.write(map_file)
        f.close()
        Phylo.write(new_tree, out_name, 'newick')
        """
