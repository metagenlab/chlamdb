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

urlpatterns = [        
    url('^robots.txt$', TemplateView.as_view(template_name='robots.txt', content_type='text/plain')),
    url(r'^sitemap$', views.sitemap, name="sitemap"),
    url(r'^home/$', views.home, name="home"),
    url(r'^blastnr_euk/$', views.blastnr_euk, name="blastnr_euk"),
    url(r'^cog_barchart/$', views.cog_barchart, name="cog_barchart"),
    url(r'^blast_sets/$', views.blast_sets, name="blast_sets"),
    url(r'^orthogroup_KO_COG/$', views.orthogroup_KO_COG, name="orthogroup_KO_COG"),
    url(r'^pan_genome/([a-zA-Z0-9_]+)', views.pan_genome, name="pan_genome"),
    url(r'^blastnr_overview/$', views.blastnr_overview, name="blastnr_overview"),
    url(r'^refseq_swissprot_tree/([a-zA-Z0-9_]+)', views.refseq_swissprot_tree, name="refseq_swissprot_tree"),
    url(r'^COG_phylo_heatmap/([a-zA-Z0-9_\-]+)', views.COG_phylo_heatmap, name="COG_phylo_heatmap"),
    url(r'^module_barchart/$', views.module_barchart, name="module_barchart"),
    url(r'^get_fasta/$', views.get_fasta, name="get_fasta"),
    url(r'^effector_pred/$', views.effector_pred, name="effector_pred"),
    url(r'^paralogs/$', views.paralogs, name="paralogs"),
    url(r'^blastnr_top_non_phylum/$', views.blastnr_top_non_phylum, name="blastnr_top_non_phylum"),
    url(r'^fasta/$', views.fasta, name="fasta"),
    url(r'^species_specific_groups/$', views.species_specific_groups, name="species_specific_groups"),
    url(r'^ko2fasta/([a-zA-Z0-9_]+)/([a-zA-Z0-9_]+)', views.ko2fasta, name="ko2fasta"),
    url(r'^ko2fasta/([a-zA-Z0-9_]+)/', views.ko2fasta, name="ko2fasta"),
    url(r'^pfam2fasta/([a-zA-Z0-9_]+)/', views.pfam2fasta, name="pfam2fasta"),
    url(r'^plot_heatmap/([a-zA-Z0-9_\-]+)', views.plot_heatmap, name="plot_heatmap"),
    url(r'^get_newick_tree/([a-zA-Z0-9_\-]+)/([a-zA-Z0-9_\-]+)', views.get_newick_tree, name="get_newick_tree"),
    url(r'^get_orthogroup_fasta/([a-zA-Z0-9_\-]+)/([a-zA-Z0-9_\-]+)', views.get_orthogroup_fasta, name="get_orthogroup_fasta"),
    url(r'^phylogeny/([a-zA-Z0-9_\.]+)', views.phylogeny, name="phylogeny"),
    url(r'^blastnr_barchart/$', views.blastnr_barchart, name="blastnr_barchart"),
    url(r'^get_cog/([a-zA-Z0-9_\.]+)/([a-zA-Z0-9_\.\%]+)$', views.get_cog, name="get_cog"),
    url(r'^module_cat_info/([a-zA-Z0-9_\.]+)/([a-zA-Z0-9_\.\+-]+)$', views.module_cat_info, name="module_cat_info"),
    url(r'^blastnr_cat_info/([a-zA-Z0-9_\.\-]+)/([a-zA-Z0-9_\.\-]+)/([a-zA-Z0-9_\.\+-]+)$', views.blastnr_cat_info, name="blastnr_cat_info"),
    url(r'^ko_venn_subset/([a-zA-Z0-9_\.\+-]+)$', views.ko_venn_subset, name="ko_venn_subset"),
    url(r'^homology/$', views.homology, name="homology"),
    url(r'^api_vs_16S_identity/$', views.api_vs_16S_identity, name="api_vs_16S_identity"),
    url(r'^hydropathy/([a-zA-Z0-9_\.\+-]+)$', views.hydropathy, name="hydropathy"),
    url(r'^hmm2circos/$', views.hmm2circos, name="hmm2circos"),
    url(r'^pfam_tree/([a-zA-Z0-9_\.]+)$', views.pfam_tree, name="pfam_tree"),
    url(r'^TM_tree/([a-zA-Z0-9_\.\-]+)$', views.TM_tree, name="TM_tree"),
    url(r'^search_taxonomy/$', views.search_taxonomy, name="search_taxonomy"),
    url(r'^priam_kegg/$', views.priam_kegg, name="priam_kegg"),
    url(r'^locusx/([a-zA-Z0-9_\.\-]+)/([a-zA-Z0-9_\.\-]+)', views.locusx, name="locusx"),
    url(r'^search_bar$', views.search_bar, name="search_bar"),
    url(r'^search_bar/([a-zA-Z0-9_\.\-]+)/([a-zA-Z0-9_\.\-]+)', views.search_bar, name="search_bar"),
    url(r'^locusx/([a-zA-Z0-9_\.\-]+)', views.locusx, name="locusx"),
    url(r'^orthogroup/([a-zA-Z0-9_\.\-]+)', views.orthogroup, name="orthogroup"),
    url(r'^locusx$', views.locusx, name="locusx"),
    url(r'^pairwiseid/$', views.pairwiseid, name="pairwiseid"),
    url(r'^get_BBH_non_chlamydiae_taxonomy/$', views.get_BBH_non_chlamydiae_taxonomy, name="get_BBH_non_chlamydiae_taxonomy"),
    url(r'^pairwiseCDS_length/$', views.pairwiseCDS_length, name="pairwiseCDS_length"),
    url(r'^blastnr/([a-zA-Z0-9_\.\-]+)$', views.blastnr, name="blastnr"),
    url(r'^blastswissprot/([a-zA-Z0-9_\.\-]+)$', views.blastswissprot, name="blastswissprot"),
    url(r'^plot_region/$', views.plot_region, name="plot_region"),
    url(r'^orthogroups/$', views.orthogroups, name="orthogroups"),
    url(r'^multipleGC/$', views.multipleGC, name="multipleGC"),
    url(r'^multiple_codon_usage/$', views.multiple_codon_usage, name="multiple_codon_usage"),
    url(r'^circos/$', views.circos, name="circos"),
    url(r'^circos_blastnr/$', views.circos_blastnr, name="circos_blastnr"),
    url(r'^circos_main/$', views.circos_main, name="circos_main"),
    url(r'^orthogroup_list_cog_barchart/$', views.orthogroup_list_cog_barchart, name="orthogroup_list_cog_barchart"),
    url(r'^orthogroup_list_cog_barchart/([a-zA-Z0-9_]+)/$', views.orthogroup_list_cog_barchart, name="orthogroup_list_cog_barchart"),
    url(r'^blast/$', views.blast, name="blast"),
    url(r'^prot_length_barchart/$', views.prot_length_barchart, name="prot_length_barchart"),
    url(r'^extract_orthogroup/$', views.extract_orthogroup, name="extract_orthogroup"),
    url(r'^extract_orthogroup/([a-zA-Z0-9_]+)$$', views.extract_orthogroup, name="extract_orthogroup"),
    url(r'^venn_orthogroup/$', views.venn_orthogroup, name="venn_orthogroup"),
    url(r'^extract_cog/([a-zA-Z0-9_]+)$$', views.extract_cog, name="extract_cog"),
    url(r'^genome_annotation/([a-zA-Z0-9_\.\-]+)$$', views.genome_annotation, name="genome_annotation"),
    url(r'^extract_cog/$', views.extract_cog, name="extract_cog"),
    url(r'^venn_cog/$', views.venn_cog, name="venn_cog"),
    url(r'^cog_venn_subset/([A-Z])$', views.cog_venn_subset, name="cog_venn_subset"),
    url(r'^venn_cog/([a-zA-Z0-9_]+)$$', views.venn_cog, name="venn_cog"),
    url(r'^venn_ko/$', views.venn_ko, name="venn_ko"),
    url(r'^extract_pfam/$', views.extract_pfam, name="extract_pfam"),
    url(r'^extract_pfam/([a-zA-Z0-9_]+)$', views.extract_pfam, name="extract_pfam"),
    url(r'^extract_ko/$', views.extract_ko, name="extract_ko"),
    url(r'^extract_ko/([a-zA-Z0-9_]+)$', views.extract_ko, name="extract_ko"),
    url(r'^venn_pfam/$', views.venn_pfam, name="venn_pfam"),
    url(r'^blast_profile/$', views.blast_profile, name="blast_profile"),
    url(r'^KEGG_mapp_ko/([a-zA-Z0-9_\.]+)$', views.KEGG_mapp_ko, name="KEGG_mapp_ko"),
    url(r'^KEGG_mapp_ko/([a-zA-Z0-9_\.]+)/([0-9]+)$', views.KEGG_mapp_ko, name="KEGG_mapp_ko"),
    url(r'^KEGG_module_map/([a-zA-Z0-9_\.]+)$', views.KEGG_module_map, name="KEGG_module_map"),
    url(r'^motif_search/$', views.motif_search, name="motif_search"),
    url(r'^primer_search/$', views.primer_search, name="primer_search"),
    url(r'^circos2genomes/$', views.circos2genomes, name="circos2genomes"),
    url(r'^download_COG/$', views.download_COG, name="download_COG"),
    url(r'^kegg_module_subcat$', views.kegg_module_subcat, name="kegg_module_subcat"),
    url(r'^kegg_module/$', views.kegg_module, name="kegg_module"),
    url(r'^module_comparison/$', views.module_comparison, name="module_comparison"),
    url(r'^pfam_comparison', views.pfam_comparison, name="pfam_comparison"),
    url(r'^ko_comparison', views.ko_comparison, name="ko_comparison"),
    url(r'^identity_heatmap', views.identity_heatmap, name="identity_heatmap"),
    url(r'^orthogroup_comparison', views.orthogroup_comparison, name="orthogroup_comparison"),
    url(r'^locus2locus/$', views.locus2locus, name="locus2locus"),
    url(r'^kegg_multi/([a-zA-Z0-9_]+)/([a-zA-Z0-9_]+)/$', views.kegg_multi, name="kegg_multi"),
    url(r'^orthogroup_identity/([a-zA-Z0-9_\.]+)/([a-zA-Z0-9_\.]+)$', views.orthogroup_identity, name="orthogroup_identity"),
    url(r'^ortho_id_plot/([a-zA-Z0-9_\.]+)$', views.ortho_id_plot, name='ortho_id_plot'),
    url(r'^multiple_orthogroup_heatmap/([a-zA-Z0-9_\.]+)$', views.multiple_orthogroup_heatmap, name="multiple_orthogroup_heatmap"),
    url(r'^multiple_COGs_heatmap', views.multiple_COGs_heatmap, name="multiple_COGs_heatmap"),
    url(r'^plot_neighborhood/([a-zA-Z0-9_\.\-]+)$', views.plot_neighborhood, name="plot_neighborhood"),
    url(r'^logout/$', logout, {'next_page': '/'},),
    url(r'^about$', views.about, name="about"),
    url(r'^help', views.help, name="help"),
    url(r'^genomic_locus_tag', views.genomic_locus_tag_infos, name="genomic_locus_tag_infos"),

    url(r'^fam_pfam/(PF[0-9]+)$', views.fam_pfam, name="fam_pfam"),
    url(r'^fam_cog/(COG[0-9]+)$', views.fam_cog, name="fam_cog"),
    url(r'^fam_ko/(K[0-9]+)$', views.fam_ko, name="fam_ko"),

    url(r'^get_record/([a-zA-Z0-9_]+)/([a-zA-Z0-9_]+)', views.get_record, name="get_record"),
    url(r'^get-task-info/', views.get_task_info),
    url(r'^docs/(?P<path>.*)$', static.serve, {'document_root': settings.DOCS_ROOT}),
    url(r'^favicon\.ico$', favicon_view),
    url(r'^.*$', views.home, name="home"),
]
