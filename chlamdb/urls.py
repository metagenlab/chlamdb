from django.urls import re_path
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
    re_path('^robots.txt$', TemplateView.as_view(template_name='robots.txt', content_type='text/plain')),
    re_path(r'^sitemap$', views.sitemap, name="sitemap"),
    re_path(r'^home/$', views.home, name="home"),
    re_path(r'^cog_barchart/$', views.cog_barchart, name="cog_barchart"),
    re_path(r'^pan_genome/([a-zA-Z0-9_]+)', views.pan_genome, name="pan_genome"),
    re_path(r'^COG_phylo_heatmap/([a-zA-Z0-9_\-]+)', views.COG_phylo_heatmap, name="COG_phylo_heatmap"),
    re_path(r'^module_barchart/$', views.module_barchart, name="module_barchart"),
    re_path(r'^plot_heatmap/([a-zA-Z0-9_\-]+)', views.plot_heatmap, name="plot_heatmap"),
    re_path(r'^get_cog/([a-zA-Z0-9_\.]+)/([a-zA-Z0-9_\.\%]+)$', views.get_cog, name="get_cog"),
    re_path(r'^module_cat_info/([a-zA-Z0-9_\.]+)/([a-zA-Z0-9_\.\+-]+)$', views.module_cat_info, name="module_cat_info"),
    re_path(r'^ko_venn_subset/([a-zA-Z0-9_\.\+-]+)$', views.ko_venn_subset, name="ko_venn_subset"),
    re_path(r'^hydropathy/([a-zA-Z0-9_\.\+-]+)$', views.hydropathy, name="hydropathy"),
    re_path(r'^priam_kegg/$', views.priam_kegg, name="priam_kegg"),
    re_path(r'^locusx/([a-zA-Z0-9_\.\-]+)/([a-zA-Z0-9_\.\-]+)', views.locusx, name="locusx"),
    re_path(r'^search_bar$', views.search_bar, name="search_bar"),
    re_path(r'^search_bar/([a-zA-Z0-9_\.\-]+)/([a-zA-Z0-9_\.\-]+)', views.search_bar, name="search_bar"),
    re_path(r'^locusx/([a-zA-Z0-9_\.\-]+)', views.locusx, name="locusx"),
    re_path(r'^orthogroup/([a-zA-Z0-9_\.\-]+)', views.orthogroup, name="orthogroup"),
    re_path(r'^locusx$', views.locusx, name="locusx"),
    re_path(r'^plot_region/$', views.plot_region, name="plot_region"),
    re_path(r'^circos/$', views.circos, name="circos"),
    re_path(r'^circos_main/$', views.circos_main, name="circos_main"),
    re_path(r'^orthogroup_list_cog_barchart/$', views.orthogroup_list_cog_barchart, name="orthogroup_list_cog_barchart"),
    re_path(r'^orthogroup_list_cog_barchart/([a-zA-Z0-9_]+)/$', views.orthogroup_list_cog_barchart, name="orthogroup_list_cog_barchart"),
    re_path(r'^blast/$', views.blast, name="blast"),
    re_path(r'^extract_orthogroup/$', views.extract_orthogroup, name="extract_orthogroup"),
    re_path(r'^extract_orthogroup/([a-zA-Z0-9_]+)$$', views.extract_orthogroup, name="extract_orthogroup"),
    re_path(r'^venn_orthogroup/$', views.venn_orthogroup, name="venn_orthogroup"),
    re_path(r'^extract_cog/([a-zA-Z0-9_]+)$$', views.extract_cog, name="extract_cog"),
    re_path(r'^extract_cog/$', views.extract_cog, name="extract_cog"),
    re_path(r'^venn_cog/$', views.venn_cog, name="venn_cog"),
    re_path(r'^cog_venn_subset/([A-Z])$', views.cog_venn_subset, name="cog_venn_subset"),
    re_path(r'^venn_cog/([a-zA-Z0-9_]+)$$', views.venn_cog, name="venn_cog"),
    re_path(r'^venn_ko/$', views.venn_ko, name="venn_ko"),
    re_path(r'^extract_pfam/$', views.extract_pfam, name="extract_pfam"),
    re_path(r'^extract_pfam/([a-zA-Z0-9_]+)$', views.extract_pfam, name="extract_pfam"),
    re_path(r'^extract_ko/$', views.extract_ko, name="extract_ko"),
    re_path(r'^extract_ko/([a-zA-Z0-9_]+)$', views.extract_ko, name="extract_ko"),
    re_path(r'^venn_pfam/$', views.venn_pfam, name="venn_pfam"),
    re_path(r'^KEGG_mapp_ko/([a-zA-Z0-9_\.]+)$', views.KEGG_mapp_ko, name="KEGG_mapp_ko"),
    re_path(r'^KEGG_mapp_ko/([a-zA-Z0-9_\.]+)/([0-9]+)$', views.KEGG_mapp_ko, name="KEGG_mapp_ko"),
    re_path(r'^KEGG_module_map/([a-zA-Z0-9_\.]+)$', views.KEGG_module_map, name="KEGG_module_map"),
    re_path(r'^kegg_module_subcat$', views.kegg_module_subcat, name="kegg_module_subcat"),
    re_path(r'^kegg_module/$', views.kegg_module, name="kegg_module"),
    re_path(r'^module_comparison/$', views.module_comparison, name="module_comparison"),
    re_path(r'^pfam_comparison', views.pfam_comparison, name="pfam_comparison"),
    re_path(r'^ko_comparison', views.ko_comparison, name="ko_comparison"),
    re_path(r'^orthogroup_comparison', views.orthogroup_comparison, name="orthogroup_comparison"),
    re_path(r'^logout/$', logout, {'next_page': '/'},),
    re_path(r'^about$', views.about, name="about"),
    re_path(r'^help', views.help, name="help"),
    re_path(r'^fam_pfam/(PF[0-9]+)$', views.fam_pfam, name="fam_pfam"),
    re_path(r'^fam_cog/(COG[0-9]+)$', views.fam_cog, name="fam_cog"),
    re_path(r'^fam_ko/(K[0-9]+)$', views.fam_ko, name="fam_ko"),
    re_path(r'^get-task-info/', views.get_task_info),
    re_path(r'^favicon\.ico$', favicon_view),
    re_path(r'^FAQ', views.faq, name='FAQ'),
    re_path(r'^phylogeny_intro', views.phylogeny_intro, name='phylogeny_intro'),
    re_path(r'^genomes_intro', views.genomes_intro, name='genomes_intro'),
    re_path(r'^extract_contigs/([0-9]+)', views.extract_contigs, name='extract_contigs'),
    re_path(r'^extract_region', views.extract_region, name='extract_region'),
    re_path(r'^.*$', views.home, name="home"),
    #re_path(r'^FAQ',TemplateView.as_view(template_name='FAQ.html')),
]


