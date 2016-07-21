from django.conf.urls import patterns, url
from django.conf.urls import include
from django.views.generic import TemplateView

urlpatterns = patterns('chlamdb.views',
                       url(r'^home/([a-zA-Z0-9_]+)$', 'home', name="home"),
                       url(r'^homology/([a-zA-Z0-9_]+)$', 'homology', name="homology"),
                       url(r'^pfam_tree/([a-zA-Z0-9_]+)/([a-zA-Z0-9_\.]+)$', 'pfam_tree', name="pfam_tree"),
                       url(r'^orthogroup_conservation_tree/([a-zA-Z0-9_]+)/([a-zA-Z0-9_\.]+)$', 'orthogroup_conservation_tree', name="orthogroup_conservation_tree"),
                       url(r'^search_taxonomy/([a-zA-Z0-9_]+)$', 'search_taxonomy', name="search_taxonomy"),
                       url(r'^priam_kegg/([a-zA-Z0-9_]+)$', 'priam_kegg', name="priam_kegg"),
                       url(r'^locusx/([a-zA-Z0-9_]+)/([a-zA-Z0-9_\.]+)/([a-zA-Z0-9_\.]+)', 'locusx', name="locusx"),
                       url(r'^locusx/([a-zA-Z0-9_]+)/([a-zA-Z0-9_\.]+)', 'locusx', name="locusx"),
                       url(r'^locusx/([a-zA-Z0-9_]+)/', 'locusx', name="locusx"),
                       url(r'^sunburst/([a-zA-Z0-9_]+)/([a-zA-Z0-9_\.]+)$', 'sunburst', name="sunburst"),
                       url(r'^fam/([a-zA-Z0-9_]+)/([a-zA-Z0-9_\.]+)/([a-zA-Z0-9_\.]+)$', 'fam', name="fam"),
                       url(r'^blastnr/([a-zA-Z0-9_]+)/([a-zA-Z0-9_]+)$', 'blastnr', name="blastnr"),
                       url(r'^plot_region/([a-zA-Z0-9_]+)$', 'plot_region', name="plot_region"),
                       url(r'^orthogroups/$', 'orthogroups', name="orthogroups"),
                       url(r'^circos/([a-zA-Z0-9_]+)$', 'circos', name="circos"),
                       url(r'^circos_main/([a-zA-Z0-9_]+)$', 'circos_main', name="circos_main"),
                       url(r'^extract_region/([a-zA-Z0-9_]+)$', 'extract_region', name="extract_region"),
                       url(r'^circos_homology/([a-zA-Z0-9_]+)$', 'circos_homology', name="circos_homology"),
                       url(r'^search/([a-zA-Z0-9_]+)$', 'search', name="search"),
                       url(r'^interpro/([a-zA-Z0-9_]+)$', 'interpro', name="interpro"),
                       url(r'^blast/([a-zA-Z0-9_]+)$$', 'blast', name="blast"),
                       url(r'^mummer/([a-zA-Z0-9_]+)$$', 'mummer', name="mummer"),
                       url(r'^extract_orthogroup/([a-zA-Z0-9_]+)$$', 'extract_orthogroup', name="extract_orthogroup"),
                       url(r'^venn_orthogroup/([a-zA-Z0-9_]+)$$', 'venn_orthogroup', name="venn_orthogroup"),
                       url(r'^extract_cog/([a-zA-Z0-9_]+)$$', 'extract_cog', name="extract_cog"),
                       url(r'^venn_cog/([a-zA-Z0-9_]+)$$', 'venn_cog', name="venn_cog"),
                       url(r'^extract_interpro/([a-zA-Z0-9_]+)$$', 'extract_interpro', name="extract_interpro"),
                       url(r'^venn_interpro/([a-zA-Z0-9_]+)$$', 'venn_interpro', name="venn_interpro"),
                       url(r'^extract_pfam/([a-zA-Z0-9_]+)$$', 'extract_pfam', name="extract_pfam"),
                       url(r'^extract_EC/([a-zA-Z0-9_]+)$$', 'extract_EC', name="extract_EC"),
                       url(r'^venn_pfam/([a-zA-Z0-9_]+)$$', 'venn_pfam', name="venn_pfam"),
                       url(r'^venn_EC/([a-zA-Z0-9_]+)$$', 'venn_EC', name="venn_EC"),
                       url(r'^KEGG_mapp/([a-zA-Z0-9_]+)/([a-zA-Z0-9_\.]+)$', 'KEGG_mapp', name="KEGG_mapp"),
                       url(r'^motif_search/([a-zA-Z0-9_]+)$$', 'motif_search', name="motif_search"),
                       url(r'^primer_search/([a-zA-Z0-9_]+)$$', 'primer_search', name="primer_search"),
                       url(r'^choose_db/$', 'choose_db', name="choose_db"),
                       url(r'^circos2genomes/([a-zA-Z0-9_]+)/$', 'circos2genomes', name="circos2genomes"),
                       url(r'^metabo_overview/([a-zA-Z0-9_]+)/$', 'metabo_overview', name="metabo_overview"),
                       url(r'^metabo_comparison/([a-zA-Z0-9_]+)/$', 'metabo_comparison', name="metabo_comparison"),
                       url(r'^pfam_comparison/([a-zA-Z0-9_]+)/$', 'pfam_comparison', name="pfam_comparison"),
                       url(r'^crossplot/$', 'crossplot', name="crossplot"),
                       url(r'^login/$', 'chlamdb_login', name='chlamdb_login'),
                       #url(r'^chaining/', include('smart_selects.urls')),
                       url(r'^autocomplete/', include('autocomplete_light.urls')),
                       url(r'^logout/$', 'logout_view', name="logout_view"),
                       url(r'^orthogroup_identity/([a-zA-Z0-9_\.]+)/([a-zA-Z0-9_\.]+)/([a-zA-Z0-9_\.]+)$', 'orthogroup_identity', name="orthogroup_identity"),
                       url(r'^ortho_id_plot/([a-zA-Z0-9_\.]+)/', 'ortho_id_plot', name='ortho_id_plot'),
                       url(r'^string_page/([a-zA-Z0-9_]+)/([a-zA-Z0-9_\.]+)/([a-zA-Z0-9_\.]+)$', 'string_page', name="string_page"),
                       url(r'^multiple_orthogroup_heatmap/([a-zA-Z0-9_]+)/([a-zA-Z0-9_\.]+)$', 'multiple_orthogroup_heatmap', name="multiple_orthogroup_heatmap"),
                       url(r'^multiple_COGs_heatmap/([a-zA-Z0-9_]+)', 'multiple_COGs_heatmap', name="multiple_COGs_heatmap"),
                       url(r'^plot_neighborhood/([a-zA-Z0-9_]+)/([a-zA-Z0-9_]+)$', 'plot_neighborhood', name="plot_neighborhood"),
                       url(r'^orthogroup_annotation/([a-zA-Z0-9_]+)$', 'orthogroup_annotation', name="orthogroup_annotation"),

)
urlpatterns += patterns('',
  (r'^logout/$', 'django.contrib.auth.views.logout', {'next_page': '/'}),
)

