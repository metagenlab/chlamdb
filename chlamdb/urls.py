from django.conf.urls import url
from django.conf.urls import include
from django.views.generic import TemplateView
from views import *
import views
from django.contrib.auth.views import logout

urlpatterns = [
                       url(r'^home/([a-zA-Z0-9_\.]+)$', home, name="home"),
                       url(r'^cog_barchart/([a-zA-Z0-9_]+)', cog_barchart, name="cog_barchart"),
                       url(r'^COG_phylo_heatmap/([a-zA-Z0-9_]+)/([a-zA-Z0-9_]+)', COG_phylo_heatmap, name="COG_phylo_heatmap"),
                       url(r'^interactions_genome/([a-zA-Z0-9_]+)', interactions_genome, name="interactions_genome"),
                       url(r'^module_barchart/([a-zA-Z0-9_]+)', module_barchart, name="module_barchart"),
                       url(r'^locus_int/([a-zA-Z0-9_]+)', locus_int, name="locus_int"),
                       url(r'^module2heatmap/([a-zA-Z0-9_]+)', module2heatmap, name="module2heatmap"),
                       url(r'^get_fasta/([a-zA-Z0-9_]+)', get_fasta, name="get_fasta"),
                       url(r'^ko2fasta/([a-zA-Z0-9_]+)/([a-zA-Z0-9_]+)/([a-zA-Z0-9_]+)', ko2fasta, name="ko2fasta"),
                       url(r'^ko2fasta/([a-zA-Z0-9_]+)/([a-zA-Z0-9_]+)/', ko2fasta, name="ko2fasta"),
                       url(r'^plot_heatmap/([a-zA-Z0-9_]+)/([a-zA-Z0-9_]+)', plot_heatmap, name="plot_heatmap"),
                       url(r'^get_newick_tree/([a-zA-Z0-9_]+)/([a-zA-Z0-9_]+)', get_newick_tree, name="get_newick_tree"),
                       url(r'^get_orthogroup_fasta/([a-zA-Z0-9_]+)/([a-zA-Z0-9_]+)/([a-zA-Z0-9_]+)', get_orthogroup_fasta, name="get_orthogroup_fasta"),
                       url(r'^add_comment/([a-zA-Z0-9_]+)/([a-zA-Z0-9_\.]+)', add_comment, name="add_comment"),
                       url(r'^add_locus_int/([a-zA-Z0-9_]+)', add_locus_int, name="add_locus_int"),
                       url(r'^neig_interactions/([a-zA-Z0-9_]+)/([a-zA-Z0-9_\.]+)', neig_interactions, name="neig_interactions"),
                       url(r'^profile_interactions/([a-zA-Z0-9_]+)/([a-zA-Z0-9_\.]+)', profile_interactions, name="profile_interactions"),
                       url(r'^interactions/([a-zA-Z0-9_]+)/([a-zA-Z0-9_\.]+)', interactions, name="interactions"),
                       url(r'^blastnr_barchart/([a-zA-Z0-9_\.]+)', blastnr_barchart, name="blastnr_barchart"),
                       url(r'^get_cog/([a-zA-Z0-9_]+)/([a-zA-Z0-9_\.]+)/([a-zA-Z0-9_\.\%]+)$', get_cog, name="get_cog"),
                       url(r'^get_cog_multiple/([a-zA-Z0-9_]+)/([A-Z]+)/$', get_cog_multiple, name="get_cog_multiple"),
                       url(r'^get_cog_multiple/([a-zA-Z0-9_]+)/([a-zA-Z]+)/([a-zA-Z]+)/$', get_cog_multiple, name="get_cog_multiple"),
                       url(r'^get_orthogroup_multiple_cog/([a-zA-Z0-9_]+)/([a-zA-Z0-9_\.\%]+)', get_orthogroup_multiple_cog, name="get_orthogroup_multiple_cog"),
                       url(r'^get_ko_multiple/([a-zA-Z0-9_]+)/([a-zA-Z0-9_]+)/([a-zA-Z0-9_\.\%\+\-]+)', get_ko_multiple, name="get_ko_multiple"),
                       url(r'^module_cat_info/([a-zA-Z0-9_]+)/([a-zA-Z0-9_\.]+)/([a-zA-Z0-9_\.\+-]+)$', module_cat_info, name="module_cat_info"),
                       url(r'^blastnr_cat_info/([a-zA-Z0-9_]+)/([a-zA-Z0-9_\.]+)/([a-zA-Z0-9_\.]+)/([a-zA-Z0-9_\.\+-]+)$', blastnr_cat_info, name="blastnr_cat_info"),
                       url(r'^cog_venn_subset/([a-zA-Z0-9_]+)/([a-zA-Z0-9_\.\%-]+)$', cog_venn_subset, name="cog_venn_subset"),
                       url(r'^ko_venn_subset/([a-zA-Z0-9_]+)/([a-zA-Z0-9_\.\+-]+)$', ko_venn_subset, name="ko_venn_subset"),
                       url(r'^homology/([a-zA-Z0-9_]+)$', homology, name="homology"),
                       url(r'^hmm/([a-zA-Z0-9_]+)$', hmm, name="hmm"),
                       url(r'^hmm2circos/([a-zA-Z0-9_]+)$', hmm2circos, name="hmm2circos"),
                       url(r'^pfam_tree/([a-zA-Z0-9_]+)/([a-zA-Z0-9_\.]+)$', pfam_tree, name="pfam_tree"),
                       url(r'^TM_tree/([a-zA-Z0-9_]+)/([a-zA-Z0-9_\.]+)$', TM_tree, name="TM_tree"),
                       url(r'^orthogroup_conservation_tree/([a-zA-Z0-9_]+)/([a-zA-Z0-9_\.]+)$', orthogroup_conservation_tree, name="orthogroup_conservation_tree"),
                       url(r'^search_taxonomy/([a-zA-Z0-9_]+)$', search_taxonomy, name="search_taxonomy"),
                       url(r'^priam_kegg/([a-zA-Z0-9_]+)$', priam_kegg, name="priam_kegg"),
                       url(r'^locusx/([a-zA-Z0-9_]+)/([a-zA-Z0-9_\.]+)/([a-zA-Z0-9_\.]+)', locusx, name="locusx"),
                       url(r'^locusx/([a-zA-Z0-9_]+)/([a-zA-Z0-9_\.]+)', locusx, name="locusx"),
                       url(r'^locusx/([a-zA-Z0-9_]+)/', locusx, name="locusx"),
                       url(r'^compare_homologs/([a-zA-Z0-9_]+)/', compare_homologs, name="compare_homologs"),
                       url(r'^pairwiseid/([a-zA-Z0-9_]+)/', pairwiseid, name="pairwiseid"),
                       url(r'^sunburst/([a-zA-Z0-9_]+)/([a-zA-Z0-9_\.]+)$', sunburst, name="sunburst"),
                       url(r'^fam/([a-zA-Z0-9_]+)/([a-zA-Z0-9_\.]+)/([a-zA-Z0-9_\.]+)$', fam, name="fam"),
                       url(r'^blastnr/([a-zA-Z0-9_]+)/([a-zA-Z0-9_\.]+)$', blastnr, name="blastnr"),
                       url(r'^blastswissprot/([a-zA-Z0-9_]+)/([a-zA-Z0-9_\.]+)$', blastswissprot, name="blastswissprot"),
                       url(r'^plot_region/([a-zA-Z0-9_]+)$', plot_region, name="plot_region"),
                       url(r'^plot_region_direct/([a-zA-Z0-9_]+)/([a-zA-Z0-9_]+)$', plot_region_direct, name="plot_region_direct"),
                       url(r'^orthogroups/$', orthogroups, name="orthogroups"),
                       url(r'^circos/([a-zA-Z0-9_]+)$', circos, name="circos"),
                       url(r'^circos_main/([a-zA-Z0-9_]+)$', circos_main, name="circos_main"),
                       url(r'^orthogroup_list_cog_barchart/([a-zA-Z0-9_]+)$', orthogroup_list_cog_barchart, name="orthogroup_list_cog_barchart"),
                       url(r'^orthogroup_list_cog_barchart/([a-zA-Z0-9_]+)/([a-zA-Z0-9_]+)/$', orthogroup_list_cog_barchart, name="orthogroup_list_cog_barchart"),
                       url(r'^cog_subset_barchart/([a-zA-Z0-9_]+)$', cog_subset_barchart, name="cog_subset_barchart"),
                       url(r'^cog_subset_barchart/([a-zA-Z0-9_]+)/([a-zA-Z0-9_]+)$', cog_subset_barchart, name="cog_subset_barchart"),
                       url(r'^ko_subset_barchart/([a-zA-Z0-9_]+)/([a-zA-Z0-9_]+)$', ko_subset_barchart, name="ko_subset_barchart"),
                       url(r'^extract_region/([a-zA-Z0-9_]+)$', extract_region, name="extract_region"),
                       url(r'^circos_homology/([a-zA-Z0-9_]+)$', circos_homology, name="circos_homology"),
                       url(r'^search/([a-zA-Z0-9_]+)$', search, name="search"),
                       url(r'^interpro/([a-zA-Z0-9_]+)$', interpro, name="interpro"),
                       url(r'^blast/([a-zA-Z0-9_]+)$$', blast, name="blast"),
                       url(r'^mummer/([a-zA-Z0-9_]+)$$', mummer, name="mummer"),
                       url(r'^extract_orthogroup/([a-zA-Z0-9_]+)$$', extract_orthogroup, name="extract_orthogroup"),
                       url(r'^extract_orthogroup/([a-zA-Z0-9_]+)/([a-zA-Z0-9_]+)$$', extract_orthogroup, name="extract_orthogroup"),
                       url(r'^venn_orthogroup/([a-zA-Z0-9_]+)$$', venn_orthogroup, name="venn_orthogroup"),
                       url(r'^extract_cog/([a-zA-Z0-9_]+)/([a-zA-Z0-9_]+)$$', extract_cog, name="extract_cog"),
                       url(r'^genome_annotation/([a-zA-Z0-9_]+)/([a-zA-Z0-9_\.]+)$$', genome_annotation, name="genome_annotation"),
                       url(r'^extract_cog/([a-zA-Z0-9_]+)$$', extract_cog, name="extract_cog"),
                       url(r'^venn_cog/([a-zA-Z0-9_]+)$$', venn_cog, name="venn_cog"),
                       url(r'^interpro_taxonomy/([a-zA-Z0-9_]+)$', interpro_taxonomy, name="interpro_taxonomy"),
                       url(r'^venn_cog/([a-zA-Z0-9_]+)/([a-zA-Z0-9_]+)$$', venn_cog, name="venn_cog"),
                       url(r'^venn_ko/([a-zA-Z0-9_]+)$$', venn_ko, name="venn_ko"),
                       url(r'^extract_interpro/([a-zA-Z0-9_]+)$$', extract_interpro, name="extract_interpro"),
                       url(r'^extract_interpro/([a-zA-Z0-9_]+)/([a-zA-Z0-9_]+)$$', extract_interpro, name="extract_interpro"),
                       url(r'^venn_interpro/([a-zA-Z0-9_]+)$$', venn_interpro, name="venn_interpro"),
                       url(r'^extract_pfam/([a-zA-Z0-9_]+)$$', extract_pfam, name="extract_pfam"),
                       url(r'^extract_pfam/([a-zA-Z0-9_]+)/([a-zA-Z0-9_]+)$$', extract_pfam, name="extract_pfam"),
                       url(r'^extract_EC/([a-zA-Z0-9_]+)$$', extract_EC, name="extract_EC"),
                       url(r'^extract_EC/([a-zA-Z0-9_]+)/([a-zA-Z0-9_]+)$$', extract_EC, name="extract_EC"),
                       url(r'^extract_ko/([a-zA-Z0-9_]+)$$', extract_ko, name="extract_ko"),
                       url(r'^extract_ko/([a-zA-Z0-9_]+)/([a-zA-Z0-9_]+)$$', extract_ko, name="extract_ko"),
                       url(r'^venn_pfam/([a-zA-Z0-9_]+)$$', venn_pfam, name="venn_pfam"),
                       url(r'^venn_EC/([a-zA-Z0-9_]+)$$', venn_EC, name="venn_EC"),
                       url(r'^blast_profile/([a-zA-Z0-9_]+)$$', blast_profile, name="blast_profile"),
                       url(r'^KEGG_mapp/([a-zA-Z0-9_]+)/([a-zA-Z0-9_\.]+)$', KEGG_mapp, name="KEGG_mapp"),
                       url(r'^KEGG_mapp_ko/([a-zA-Z0-9_]+)/([a-zA-Z0-9_\.]+)$', KEGG_mapp_ko, name="KEGG_mapp_ko"),
                       url(r'^KEGG_module_map/([a-zA-Z0-9_]+)/([a-zA-Z0-9_\.]+)$', KEGG_module_map, name="KEGG_module_map"),
                       url(r'^motif_search/([a-zA-Z0-9_]+)$$', motif_search, name="motif_search"),
                       url(r'^primer_search/([a-zA-Z0-9_]+)$$', primer_search, name="primer_search"),
                       url(r'^choose_db/$', choose_db, name="choose_db"),
                       url(r'^circos2genomes/([a-zA-Z0-9_]+)/$', circos2genomes, name="circos2genomes"),
                       url(r'^metabo_overview/([a-zA-Z0-9_]+)/$', metabo_overview, name="metabo_overview"),
                       url(r'^kegg_module/([a-zA-Z0-9_]+)/$', kegg_module, name="kegg_module"),
                       url(r'^kegg_module_subcat/([a-zA-Z0-9_]+)/$', kegg_module_subcat, name="kegg_module_subcat"),
                       url(r'^kegg_pathway_heatmap/([a-zA-Z0-9_]+)/$', kegg_pathway_heatmap, name="kegg_pathway_heatmap"),
                       url(r'^metabo_comparison/([a-zA-Z0-9_]+)/$', metabo_comparison, name="metabo_comparison"),
                       url(r'^metabo_comparison_ko/([a-zA-Z0-9_]+)/$', metabo_comparison_ko, name="metabo_comparison_ko"),
                       url(r'^module_comparison/([a-zA-Z0-9_]+)/$', module_comparison, name="module_comparison"),
                       url(r'^pfam_comparison/([a-zA-Z0-9_]+)/$', pfam_comparison, name="pfam_comparison"),
                       url(r'^ko_comparison/([a-zA-Z0-9_]+)/$', ko_comparison, name="ko_comparison"),
                       url(r'^locus_list2orthogroups/([a-zA-Z0-9_]+)/$', locus_list2orthogroups, name="locus_list2orthogroups"),
                       url(r'^orthogroup_comparison/([a-zA-Z0-9_]+)/$', orthogroup_comparison, name="orthogroup_comparison"),
                       url(r'^crossplot/$', crossplot, name="crossplot"),
                       url(r'^login/$', views.chlamdb_login, name='chlamdb_login'),
                       #url(r'^chaining/', include('smart_selects.urls')),
                       url(r'^logout/$', logout_view, name="logout_view"),
                       url(r'^orthogroup_identity/([a-zA-Z0-9_\.]+)/([a-zA-Z0-9_\.]+)/([a-zA-Z0-9_\.]+)$', orthogroup_identity, name="orthogroup_identity"),
                       url(r'^ortho_id_plot/([a-zA-Z0-9_\.]+)/', ortho_id_plot, name='ortho_id_plot'),
                       url(r'^string_page/([a-zA-Z0-9_]+)/([a-zA-Z0-9_\.]+)/([a-zA-Z0-9_\.]+)$', string_page, name="string_page"),
                       url(r'^multiple_orthogroup_heatmap/([a-zA-Z0-9_]+)/([a-zA-Z0-9_\.]+)$', multiple_orthogroup_heatmap, name="multiple_orthogroup_heatmap"),
                       url(r'^multiple_COGs_heatmap/([a-zA-Z0-9_]+)', multiple_COGs_heatmap, name="multiple_COGs_heatmap"),
                       url(r'^plot_neighborhood/([a-zA-Z0-9_\.]+)/([a-zA-Z0-9_\.]+)$', plot_neighborhood, name="plot_neighborhood"),
                       url(r'^orthogroup_annotation/([a-zA-Z0-9_]+)/([a-zA-Z0-9_]+)$', orthogroup_annotation, name="orthogroup_annotation"),
                       url(r'^locus_annotation/([a-zA-Z0-9_]+)/([a-zA-Z0-9_]+)$', locus_annotation, name="locus_annotation"),
                       url(r'^logout/$', logout, {'next_page': '/'}),

]

