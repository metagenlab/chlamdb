from django.conf.urls import patterns, include, url
#import autocomplete_light
#autocomplete_light.autodiscover()

from django.contrib import admin
admin.autodiscover()

urlpatterns = patterns('',
    url(r'^admin/', include(admin.site.urls)),
    url(r'^chlamdb/', include('chlamdb.urls')),
)
