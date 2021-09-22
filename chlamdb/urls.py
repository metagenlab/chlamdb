from django.conf.urls import url
from . import views
from django.contrib.auth import logout
from django.contrib.sitemaps.views import sitemap
from django.views.generic import TemplateView
from django.urls import include, path
from django.conf import settings
from django.urls import reverse
from django.contrib.sitemaps import Sitemap
from django.views import static
from django.views.generic.base import RedirectView

favicon_view = RedirectView.as_view(url='/assets/favicon.ico', permanent=True)

class ViewSitemap(Sitemap):
    """Reverse 'static' views for XML sitemap."""

    def items(self):
        # Return list of url names for views to include in sitemap
        return ['home', 'about']

    def location(self, item):
        #site = Site(domain='chlamdb.ch', name='chlamdb.ch')
        return reverse(item)


sitemaps = {'views': ViewSitemap}

# url(r'^sitemap.xml$', sitemap, {'sitemaps': sitemaps}),

urlpatterns = [        url('^robots.txt$', TemplateView.as_view(template_name='robots.txt', content_type='text/plain')),
                       url(r'^sitemap$', views.sitemap, name="sitemap"),
                       url(r'^home/$', views.home, name="home"),
                       url(r'^blastnr_euk/$', views.blastnr_euk, name="blastnr_euk"),
                       url(r'^pmid/([a-zA-Z0-9_]+)$', views.pmid, name="pmid"),
                       url(r'^cog_barchart/$', views.cog_barchart, name="cog_barchart"),
                       url(r'^blast_sets/$', views.blast_sets, name="blast_sets"),
                       url(r'^rnaseq_class/([a-zA-Z0-9_]+)/([a-zA-Z0-9_]+)', views.rnaseq_class, name="rnaseq_class"),
                       url(r'^orthogroup_KO_COG/$', views.orthogroup_KO_COG, name="orthogroup_KO_COG"),
                       url(r'^pan_genome/([a-zA-Z0-9_]+)', views.pan_genome, name="pan_genome"),
                       url(r'^blastnr_overview/$', views.blastnr_overview, name="blastnr_overview"),
                       url(r'^refseq_swissprot_tree/([a-zA-Z0-9_]+)', views.refseq_swissprot_tree, name="refseq_swissprot_tree"),
                       url(r'^COG_phylo_heatmap/([a-zA-Z0-9_\-]+)', views.COG_phylo_heatmap, name="COG_phylo_heatmap"),
                       url(r'^interactions_genome/$', views.interactions_genome, name="interactions_genome"),
                       url(r'^interactions_genome_string$', views.interactions_genome_string, name="interactions_genome_string"),
                       url(r'^module_barchart/$', views.module_barchart, name="module_barchart"),
                       url(r'^locus_int/$', views.locus_int, name="locus_int"),
                       url(r'^module2heatmap/$', views.module2heatmap, name="module2heatmap"),
                       url(r'^get_fasta/$', views.get_fasta, name="get_fasta"),
                       url(r'^effector_pred/$', views.effector_pred, name="effector_pred"),
                       url(r'^transporters/$', views.transporters, name="transporters"),
                       url(r'^paralogs/$', views.paralogs, name="paralogs"),
                       url(r'^blastnr_top_non_phylum/$', views.blastnr_top_non_phylum, name="blastnr_top_non_phylum"),
                       url(r'^transporters_list/$', views.transporters_list, name="transporters_list"),
                       url(r'^fasta/$', views.fasta, name="fasta"),
                       url(r'^species_specific_groups/$', views.species_specific_groups, name="species_specific_groups"),
                       url(r'^ko2fasta/([a-zA-Z0-9_]+)/([a-zA-Z0-9_]+)', views.ko2fasta, name="ko2fasta"),
                       url(r'^ko2fasta/([a-zA-Z0-9_]+)/', views.ko2fasta, name="ko2fasta"),
                       url(r'^pfam2fasta/([a-zA-Z0-9_]+)/', views.pfam2fasta, name="pfam2fasta"),
                       url(r'^venn_candidate_effectors/$', views.venn_candidate_effectors, name="venn_candidate_effectors"),
                       url(r'^effector_predictions/([a-zA-Z0-9_]+)', views.effector_predictions, name="effector_predictions"),
                       url(r'^plot_heatmap/([a-zA-Z0-9_\-]+)', views.plot_heatmap, name="plot_heatmap"),
                       url(r'^get_newick_tree/([a-zA-Z0-9_\-]+)/([a-zA-Z0-9_\-]+)', views.get_newick_tree, name="get_newick_tree"),
                       url(r'^get_orthogroup_fasta/([a-zA-Z0-9_\-]+)/([a-zA-Z0-9_\-]+)', views.get_orthogroup_fasta, name="get_orthogroup_fasta"),
                       url(r'^add_comment/([a-zA-Z0-9_\.]+)', views.add_comment, name="add_comment"),
                       url(r'^add_locus_int/$', views.add_locus_int, name="add_locus_int"),
                       url(r'^interpro_taxonomy_with_homologs/([a-zA-Z0-9_\-]+)/([a-zA-Z0-9_\-]+)', views.interpro_taxonomy_with_homologs, name="interpro_taxonomy_with_homologs"),
                       url(r'^pfam_taxonomy_with_homologs/([\.a-zA-Z0-9_\-]+)/([\.a-zA-Z0-9_\-]+)', views.pfam_taxonomy_with_homologs,name="pfam_taxonomy_with_homologs"),
                       url(r'^neig_interactions/([a-zA-Z0-9_\.]+)', views.neig_interactions, name="neig_interactions"),
                       url(r'^profile_interactions/([a-zA-Z0-9_\.]+)/([a-zA-Z0-9_\.]+)', views.profile_interactions, name="profile_interactions"),
                       url(r'^interactions/([a-zA-Z0-9_\.]+)', views.interactions, name="interactions"),
                       url(r'^phylogeny/([a-zA-Z0-9_\.]+)', views.phylogeny, name="phylogeny"),
                       url(r'^pfam_profile/([a-zA-Z0-9_\.]+)/([a-zA-Z0-9_\.]+)', views.pfam_profile, name="pfam_profile"),
                       url(r'^get_pfam_taxon_table/([a-zA-Z0-9_\.]+)', views.get_pfam_taxon_table, name="get_pfam_taxon_table"),
                       url(r'^get_pfam_hit_list/([a-zA-Z0-9_\.]+)/(.*)/(.*)/(.*)/(.*)/(.*)', views.get_pfam_hit_list, name="get_pfam_hit_list"),
                       url(r'^get_pfam_hit_list/([a-zA-Z0-9_\.]+)/(.*)/(.*)/(.*)/(.*)', views.get_pfam_hit_list, name="get_pfam_hit_list"),
                       url(r'^get_pfam_hit_list/([a-zA-Z0-9_\.]+)/(.*)/(.*)/(.*)', views.get_pfam_hit_list, name="get_pfam_hit_list"),
                       url(r'^get_pfam_hit_list/([a-zA-Z0-9_\.]+)/(.*)/(.*)/', views.get_pfam_hit_list, name="get_pfam_hit_list"),
                       url(r'^get_pfam_hit_list/([a-zA-Z0-9_\.]+)/(.*)/', views.get_pfam_hit_list, name="get_pfam_hit_list"),
                       url(r'^get_pfam_hit_list/([a-zA-Z0-9_\.]+)/', views.get_pfam_hit_list, name="get_pfam_hit_list"),
                       url(r'^eggnog_profile/([a-zA-Z0-9_\.]+)/([a-zA-Z0-9_\.]+)', views.eggnog_profile, name="eggnog_profile"),
                       url(r'^blastnr_barchart/$', views.blastnr_barchart, name="blastnr_barchart"),
                       url(r'^get_cog/([a-zA-Z0-9_\.]+)/([a-zA-Z0-9_\.\%]+)$', views.get_cog, name="get_cog"),
                       url(r'^get_cog_multiple/([A-Z]+)/$', views.get_cog_multiple, name="get_cog_multiple"),
                       url(r'^get_cog_multiple/([a-zA-Z]+)/([a-zA-Z]+)/$', views.get_cog_multiple, name="get_cog_multiple"),
                       url(r'^get_orthogroup_multiple_cog/([a-zA-Z0-9_\.\%]+)', views.get_orthogroup_multiple_cog, name="get_orthogroup_multiple_cog"),
                       url(r'^get_ko_multiple/([a-zA-Z0-9_]+)/([a-zA-Z0-9_\.\%\+\-]+)', views.get_ko_multiple, name="get_ko_multiple"),
                       url(r'^module_cat_info/([a-zA-Z0-9_\.]+)/([a-zA-Z0-9_\.\+-]+)$', views.module_cat_info, name="module_cat_info"),
                       url(r'^blastnr_cat_info/([a-zA-Z0-9_\.\-]+)/([a-zA-Z0-9_\.\-]+)/([a-zA-Z0-9_\.\+-]+)$', views.blastnr_cat_info, name="blastnr_cat_info"),
                       url(r'^ko_venn_subset/([a-zA-Z0-9_\.\+-]+)$', views.ko_venn_subset, name="ko_venn_subset"),
                       url(r'^homology/$', views.homology, name="homology"),
                       url(r'^api_vs_16S_identity/$', views.api_vs_16S_identity, name="api_vs_16S_identity"),
                       url(r'^hydropathy/([a-zA-Z0-9_\.\+-]+)$', views.hydropathy, name="hydropathy"),
                       url(r'^aa_comp_locus/([a-zA-Z0-9_\.\+-]+)$', views.aa_comp_locus, name="aa_comp_locus"),
                       url(r'^hmm/$', views.hmm, name="hmm"),
                       url(r'^hmm2circos/$', views.hmm2circos, name="hmm2circos"),
                       url(r'^pfam_tree/([a-zA-Z0-9_\.]+)$', views.pfam_tree, name="pfam_tree"),
                       url(r'^gc_locus/([a-zA-Z0-9_\.\-]+)$', views.gc_locus, name="gc_locus"),
                       url(r'^TM_tree/([a-zA-Z0-9_\.\-]+)$', views.TM_tree, name="TM_tree"),
                       url(r'^search_taxonomy/$', views.search_taxonomy, name="search_taxonomy"),
                       url(r'^priam_kegg/$', views.priam_kegg, name="priam_kegg"),
                       url(r'^locusx/([a-zA-Z0-9_\.\-]+)/([a-zA-Z0-9_\.\-]+)', views.locusx, name="locusx"),
                       url(r'^search_bar$', views.search_bar, name="search_bar"),
                       url(r'^search_bar/([a-zA-Z0-9_\.\-]+)/([a-zA-Z0-9_\.\-]+)', views.search_bar, name="search_bar"),
                       url(r'^locusx/([a-zA-Z0-9_\.\-]+)', views.locusx, name="locusx"),
                       url(r'^orthogroup/([a-zA-Z0-9_\.\-]+)', views.orthogroup, name="orthogroup"),
                       url(r'^locusx$', views.locusx, name="locusx"),
                       url(r'^compare_homologs$', views.compare_homologs, name="compare_homologs"),
                       url(r'^pairwiseid/$', views.pairwiseid, name="pairwiseid"),
                       url(r'^get_BBH_non_chlamydiae_taxonomy/$', views.get_BBH_non_chlamydiae_taxonomy, name="get_BBH_non_chlamydiae_taxonomy"),
                       url(r'^transporters_family/([a-zA-Z0-9_\.]+)/', views.transporters_family, name="transporters_family"),
                       url(r'^pairwiseCDS_length/$', views.pairwiseCDS_length, name="pairwiseCDS_length"),
                       url(r'^sunburst/([a-zA-Z0-9_\.]+)$', views.sunburst, name="sunburst"),
                       url(r'^fam/([a-zA-Z0-9_\.]+)/([a-zA-Z0-9_\.]+)$', views.fam, name="fam"),
                       url(r'^fam_interpro/([a-zA-Z0-9_\.\:]+)/([a-zA-Z0-9_\.]+)$', views.fam_interpro, name="fam_interpro"),
                       url(r'^blastnr/([a-zA-Z0-9_\.\-]+)$', views.blastnr, name="blastnr"),
                       url(r'^blastswissprot/([a-zA-Z0-9_\.\-]+)$', views.blastswissprot, name="blastswissprot"),
                       url(r'^plot_region/$', views.plot_region, name="plot_region"),
                       url(r'^plot_region_direct/([a-zA-Z0-9_]+)$', views.plot_region_direct, name="plot_region_direct"),
                       url(r'^orthogroups/$', views.orthogroups, name="orthogroups"),
                       url(r'^multipleGC/$', views.multipleGC, name="multipleGC"),
                       url(r'^similarity_network/([a-zA-Z0-9_]+)/([a-zA-Z0-9_]+)$', views.similarity_network, name="similarity_network"),
                       url(r'^multiple_codon_usage/$', views.multiple_codon_usage, name="multiple_codon_usage"),
                       url(r'^circos/$', views.circos, name="circos"),
                       url(r'^curated_taxonomy/$', views.curated_taxonomy, name="curated_taxonomy"),
                       url(r'^circos_blastnr/$', views.circos_blastnr, name="circos_blastnr"),
                       url(r'^circos_main/$', views.circos_main, name="circos_main"),
                       url(r'^orthogroup_list_cog_barchart/$', views.orthogroup_list_cog_barchart, name="orthogroup_list_cog_barchart"),
                       url(r'^orthogroup_list_cog_barchart/([a-zA-Z0-9_]+)/$', views.orthogroup_list_cog_barchart, name="orthogroup_list_cog_barchart"),
                       url(r'^ko_subset_barchart/([a-zA-Z0-9_]+)$', views.ko_subset_barchart, name="ko_subset_barchart"),
                       url(r'^extract_region/$', views.extract_region, name="extract_region"),
                       url(r'^circos_homology/$', views.circos_homology, name="circos_homology"),
                       url(r'^search/$', views.search, name="search"),
                       url(r'^interpro/$', views.interpro, name="interpro"),
                       url(r'^blast/$', views.blast, name="blast"),
                       url(r'^edit_species_taxonomy/([0-9]+)$', views.edit_species_taxonomy, name="edit_species_taxonomy"),
                       url(r'^mummer/$', views.mummer, name="mummer"),
                       url(r'^locus_list2circos/([0-9]+)', views.locus_list2circos, name="locus_list2circos"),
                       url(r'^prot_length_barchart/$', views.prot_length_barchart, name="prot_length_barchart"),
                       url(r'^extract_orthogroup/$', views.extract_orthogroup, name="extract_orthogroup"),
                       url(r'^extract_orthogroup/([a-zA-Z0-9_]+)$$', views.extract_orthogroup, name="extract_orthogroup"),
                       url(r'^venn_orthogroup/$', views.venn_orthogroup, name="venn_orthogroup"),
                       url(r'^extract_cog/([a-zA-Z0-9_]+)$$', views.extract_cog, name="extract_cog"),
                       url(r'^genome_annotation/([a-zA-Z0-9_\.\-]+)$$', views.genome_annotation, name="genome_annotation"),
                       url(r'^extract_cog/$', views.extract_cog, name="extract_cog"),
                       url(r'^annotation_overview/$', views.annotation_overview, name="annotation_overview"),
                       url(r'^venn_cog/$', views.venn_cog, name="venn_cog"),
                       url(r'^interpro_taxonomy/$', views.interpro_taxonomy, name="interpro_taxonomy"),
                       url(r'^cog_venn_subset/([A-Z])$', views.cog_venn_subset, name="cog_venn_subset"),
                       url(r'^venn_cog/([a-zA-Z0-9_]+)$$', views.venn_cog, name="venn_cog"),
                       url(r'^venn_ko/$', views.venn_ko, name="venn_ko"),
                       url(r'^extract_interpro/$', views.extract_interpro, name="extract_interpro"),
                       url(r'^extract_interpro/([a-zA-Z0-9_]+)$', views.extract_interpro, name="extract_interpro"),
                       url(r'^venn_interpro/$', views.venn_interpro, name="venn_interpro"),
                       url(r'^extract_pfam/$', views.extract_pfam, name="extract_pfam"),
                       url(r'^extract_pfam/([a-zA-Z0-9_]+)$', views.extract_pfam, name="extract_pfam"),
                       url(r'^extract_EC/$', views.extract_EC, name="extract_EC"),
                       url(r'^extract_EC/([a-zA-Z0-9_]+)$', views.extract_EC, name="extract_EC"),
                       url(r'^extract_ko/$', views.extract_ko, name="extract_ko"),
                       url(r'^extract_ko/([a-zA-Z0-9_]+)$', views.extract_ko, name="extract_ko"),
                       url(r'^venn_pfam/$', views.venn_pfam, name="venn_pfam"),
                       url(r'^venn_EC/$', views.venn_EC, name="venn_EC"),
                       url(r'^blast_profile/$', views.blast_profile, name="blast_profile"),
                       url(r'^KEGG_mapp/([a-zA-Z0-9_\.]+)$', views.KEGG_mapp, name="KEGG_mapp"),
                       url(r'^KEGG_mapp_ko/([a-zA-Z0-9_\.]+)$', views.KEGG_mapp_ko, name="KEGG_mapp_ko"),
                       url(r'^KEGG_mapp_ko_organism/([a-zA-Z0-9_\.]+)/([a-zA-Z0-9_\.]+)$', views.KEGG_mapp_ko_organism, name="KEGG_mapp_ko_organism"),
                       url(r'^KEGG_module_map/([a-zA-Z0-9_\.]+)$', views.KEGG_module_map, name="KEGG_module_map"),
                       url(r'^motif_search/$', views.motif_search, name="motif_search"),
                       url(r'^primer_search/$', views.primer_search, name="primer_search"),
                       url(r'^choose_db/$', views.choose_db, name="choose_db"),
                       url(r'^circos2genomes/$', views.circos2genomes, name="circos2genomes"),
                       url(r'^download_COG/$', views.download_COG, name="download_COG"),
                       url(r'^pmid_associations/([a-zA-Z0-9_\.]+)/([a-zA-Z0-9_\.]+)/([a-zA-Z0-9_\.]+)$', views.pmid_associations, name="pmid_associations"),
                       url(r'^pmid_associations_orthogroups/([a-zA-Z0-9_\.]+)/([a-zA-Z0-9_\.]+)$', views.pmid_associations_orthogroups, name="pmid_associations_orthogroups"),
                       url(r'^metabo_overview/$', views.metabo_overview, name="metabo_overview"),
                       url(r'^kegg_module_subcat$', views.kegg_module_subcat, name="kegg_module_subcat"),
                       url(r'^kegg_module/$', views.kegg_module, name="kegg_module"),
                       url(r'^kegg_pathway_heatmap/$', views.kegg_pathway_heatmap, name="kegg_pathway_heatmap"),
                       url(r'^metabo_comparison/$', views.metabo_comparison, name="metabo_comparison"),
                       url(r'^metabo_comparison_ko/$', views.metabo_comparison_ko, name="metabo_comparison_ko"),
                       url(r'^module_comparison/$', views.module_comparison, name="module_comparison"),
                       url(r'^pfam_comparison', views.pfam_comparison, name="pfam_comparison"),
                       url(r'^ko_comparison', views.ko_comparison, name="ko_comparison"),
                       url(r'^identity_heatmap', views.identity_heatmap, name="identity_heatmap"),
                       url(r'^locus_list2orthogroups', views.locus_list2orthogroups, name="locus_list2orthogroups"),
                       url(r'^orthogroup_comparison', views.orthogroup_comparison, name="orthogroup_comparison"),
                       url(r'^locus2locus/$', views.locus2locus, name="locus2locus"),
                       url(r'^crossplot/$', views.crossplot, name="crossplot"),
                       url(r'^login/$', views.chlamdb_login, name='chlamdb_login'),
                       #url(r'^chaining/', include('smart_selects.urls')),
                       url(r'^logout/$', views.logout_view, name="logout_view"),
                       url(r'^kegg_multi/([a-zA-Z0-9_]+)/([a-zA-Z0-9_]+)/$', views.kegg_multi, name="kegg_multi"),
                       url(r'^orthogroup_identity/([a-zA-Z0-9_\.]+)/([a-zA-Z0-9_\.]+)$', views.orthogroup_identity, name="orthogroup_identity"),
                       url(r'^ortho_id_plot/([a-zA-Z0-9_\.]+)$', views.ortho_id_plot, name='ortho_id_plot'),
                       url(r'^string_page/([a-zA-Z0-9_\-\.]+)/([a-zA-Z0-9_\-\.]+)$', views.string_page, name="string_page"),
                       url(r'^multiple_orthogroup_heatmap/([a-zA-Z0-9_\.]+)$', views.multiple_orthogroup_heatmap, name="multiple_orthogroup_heatmap"),
                       url(r'^multiple_COGs_heatmap', views.multiple_COGs_heatmap, name="multiple_COGs_heatmap"),
                       url(r'^plot_neighborhood/([a-zA-Z0-9_\.\-]+)$', views.plot_neighborhood, name="plot_neighborhood"),
                       url(r'^orthogroup_annotation/([a-zA-Z0-9_]+)$', views.orthogroup_annotation, name="orthogroup_annotation"),
                       url(r'^locus_annotation/([a-zA-Z0-9_]+)$', views.locus_annotation, name="locus_annotation"),
                       url(r'^logout/$', logout, {'next_page': '/'},),
                       url(r'^about$', views.about, name="about"),
                       url(r'^help', views.help, name="help"),

                       url(r'^fam_cog/(COG[0-9]+)$', views.fam_cog, name="fam_cog"), 
                       url(r'^fam_ko/(K[0-9]+)$', views.fam_ko, name="fam_ko"), 


                       url(r'^get_record/([a-zA-Z0-9_]+)/([a-zA-Z0-9_]+)', views.get_record, name="get_record"),
                       url(r'^get-task-info/', views.get_task_info),
                       url(r'^docs/(?P<path>.*)$', static.serve, {'document_root': settings.DOCS_ROOT}),
                       url(r'^favicon\.ico$', favicon_view),
                       url(r'^FAQ', views.faq, name='FAQ'),
                       url(r'^phylogeny_intro', views.phylogeny_intro, name='phylogeny_intro'),
                       url(r'^genomes_intro', views.genomes_intro, name='genomes_intro'),
                       url(r'^extract_contigs/([0-9]+)', views.extract_contigs, name='extract_contigs'),
                       url(r'^.*$', views.home, name="home"),
                       #url(r'^FAQ',TemplateView.as_view(template_name='FAQ.html')),
]


