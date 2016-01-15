# -*- coding: utf-8 -*-

# todo circos gc file curently written in home directory, move it to other place
# todo save temp files in temp folder


#from django.shortcuts import render
#from datetime import datetime
from django.shortcuts import render
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import circos_orthology
# from django.core.cache import cache
#import pylibmc
#from django.core.cache import cache
import os
from django.http import HttpResponseRedirect
from django.http import HttpResponse
from forms import make_contact_form
from forms import make_plot_form
from forms import SearchForm
from forms import BiodatabaseForm
from forms import make_circos_form
from forms import make_circos2genomes_form
from forms import make_mummer_form
from forms import make_blast_form
from forms import make_crossplot_form
from forms import ConnexionForm
from forms import DBForm
from forms import make_motif_form
from forms import PCRForm
from forms import make_extract_form
from forms import make_circos_orthology_form
from forms import make_interpro_from
from forms import make_extract_region_form
from forms import make_venn_from
from forms import make_priam_form

from django.contrib.auth import logout
from django.conf import settings
from django.contrib.auth import authenticate, login
from django.contrib.auth.decorators import login_required

import manipulate_biosqldb
import mysqldb_plot_genomic_feature
from django.core.cache import get_cache
from tempfile import NamedTemporaryFile
from Bio import SeqIO
from gbk2table import Record
import models
import simplejson
import string
import random

def id_generator(size=6, chars=string.ascii_uppercase + string.ascii_lowercase + string.digits):
   return ''.join(random.choice(chars) for _ in range(size))


def extract_alphanumeric(input_string):
    from string import ascii_letters, digits
    return "".join([ch for ch in input_string if ch in (ascii_letters + digits + '_-.')])


def choose_db(request):
    server = manipulate_biosqldb.load_db()
    if request.method == 'POST':  # S'il s'agit d'une requête POST

        form = BiodatabaseForm(request.POST)  # Nous reprenons les données

        if form.is_valid():  # Nous vérifions que les données envoyées sont valides

            biodb = form.cleaned_data['biodatabase']
            envoi = True
    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = BiodatabaseForm()  # Nous créons un formulaire vide

    return render(request, 'chlamdb/choose_db.html', locals())


def create_user(username, mail, password, first_name, last_name, staff=True):
    from django.contrib.auth.models import User
    user = User.objects.create_user(username, mail, password)
    user.first_name, user.last_name = first_name, last_name
    user.is_staff = staff
    user.save()



def chlamdb_login(request):
    error = False

    if request.method == "POST":
        print "longin!"
        form = ConnexionForm(request.POST)
        if form.is_valid():
            username = form.cleaned_data["username"]
            password = form.cleaned_data["password"]
            user = authenticate(username=username, password=password)  # Nous vérifions si les données sont correctes
            if user:  # Si l'objet renvoyé n'est pas None
                login(request, user)  # nous connectons l'utilisateur
                print user, "logged"
                #return HttpResponseRedirect("/chlamdb/home")
            else: # sinon une erreur sera affichée
                error = True

        if form.is_valid():  # Nous vérifions que les données envoyées sont valides

            biodb = form.cleaned_data['biodatabase']
            envoi = True
    else:
        form = ConnexionForm()

    return render(request, 'chlamdb/login.html', locals())


def logout_view(request):
    logout(request)
    return render(request, 'chlamdb/logout.html', locals())

@login_required
def home(request, biodb):
    from Bio.SeqUtils import GC

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select accession, seq as length from bioentry as t1 ' \
          ' inner join biodatabase as t2 on t1.biodatabase_id=t2.biodatabase_id ' \
          ' inner join biosequence as t3 on t1.bioentry_id=t3.bioentry_id where t2.name="%s";' % biodb
    print 'getting accession2genome_sequence...'

    accession2genome_sequence = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'select accession, organism, count(*) from orthology_detail_%s group by accession;' % biodb

    print 'getting accession2genome_data...'
    # organism name, protein encoding ORF number
    accession2genome_data = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))


    #print '<td colspan="6"><table width="800" border=0  class=table_genomes>'

    print 'preparing data...'
    genomes_data = []

    for accession in accession2genome_data:
        print accession
        one_genome_data = []
        one_genome_data.append(accession)
        try:
            one_genome_data.append(round(GC(accession2genome_sequence[accession].replace("N", "")),2))
        except:
            one_genome_data.append('-')
        one_genome_data.append(accession2genome_data[accession][1])
        try:
            one_genome_data.append(accession2genome_sequence[accession].count(200*'N') + 1)
        except:
            one_genome_data.append('-')
        try:
            size = len(accession2genome_sequence[accession]) - (200*accession2genome_sequence[accession].count(200*'N'))
            one_genome_data.append(str(size))
        except:
            one_genome_data.append('-')
        one_genome_data.append(accession2genome_data[accession][0])
        genomes_data.append(one_genome_data)
    return render(request, 'chlamdb/home.html', locals())


def substription():
    create_user('tpillone', 'trestan.pillonel@gmail.com', 'estrella3', "Trestan", "Pillonel")

#cache = pylibmc.Client(['127.0.0.1:8000'])


@login_required
def circos_homology(request, biodb):


    cache = get_cache('default')
    print "loading db..."
    server, db = manipulate_biosqldb.load_db(biodb)
    print "db loaded..."

    circos_orthology_form_class = make_circos_orthology_form(biodb)

    if request.method == 'POST':  # S'il s'agit d'une requête POST

        form = circos_orthology_form_class(request.POST)

        if form.is_valid():

            accession = form.cleaned_data['accession']

            print "accession", accession

            sql = 'select accession from bioentry' \
                  ' inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id' \
                  ' and biodatabase.name = "%s"' % biodb



            all_accession = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]
            columns = 'orthogroup, accession, start, stop'
            sql_ref = 'select %s from orthology_detail_%s where locus_tag = "%s" or protein_id = "%s" or orthogroup = "%s"' % (columns,
                                                                                                          biodb,
                                                                                                          accession,
                                                                                                          accession, accession)




            ref_record = server.adaptor.execute_and_fetchall(sql_ref,)[0]

            orthogroup = ref_record[0]

            columns = 'accession, start, stop'
            sql_targets = 'select %s from orthology_detail_%s where orthogroup ="%s"' % (columns,
                                                                                          biodb,
                                                                                          orthogroup)

            target_records = server.adaptor.execute_and_fetchall(sql_targets,)

            print "ref_record", ref_record
            print "target_records", target_records

            record_list = []
            for accession in all_accession:
                if accession == "CP001848" or accession == "BX119912":
                    continue
                print "accession", accession
                biorecord = cache.get(biodb + "_" + accession)


                if not biorecord:
                    print biodb + "_" + accession, "NOT in memory"
                    new_record = db.lookup(accession=accession)
                    biorecord = SeqRecord(Seq(new_record.seq.data, new_record.seq.alphabet),
                                                             id=new_record.id, name=new_record.name,
                                                             description=new_record.description,
                                                             dbxrefs =new_record.dbxrefs,
                                                             features=new_record.features,
                                                             annotations=new_record.annotations)
                    record_id = biorecord.id.split(".")[0]
                    cache.set(biodb + "_" + record_id, biorecord)
                    record_list.append(biorecord)
                else:
                    record_list.append(biorecord)


            path = settings.BASE_DIR + "/assets/circos"

            circos_orthology.circos_orthology(record_list, ref_record[1:], target_records, location = path)
            circos_file = "circos/circos_ortho.svg"
            envoi_circos = True

    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = circos_orthology_form_class()  # Nous créons un formulaire vide

    return render(request, 'chlamdb/circos_homology.html', locals())

@login_required
def extract_orthogroup(request, biodb):

    cache = get_cache('default')
    print "loading db..."
    server = manipulate_biosqldb.load_db()
    print "db loaded..."
    extract_form_class = make_extract_form(biodb)

    if request.method == 'POST':  # S'il s'agit d'une requête POST

        form = extract_form_class(request.POST)  # Nous reprenons les données

        #form2 = ContactForm(request.POST)
        if 'comparison' in request.POST and form.is_valid():  # Nous vérifions que les données envoyées sont valides

            print request.POST
            print form.cleaned_data.keys()
            include = form.cleaned_data['orthologs_in']
            exclude = form.cleaned_data['no_orthologs_in']

            try:
                single_copy = request.POST['button_single_copy']
                single_copy = True
            except:
                single_copy = False
            print "single_copy", single_copy


            if single_copy:
                sql_include = ''
                if len(include) > 0:
                    for i in range(0, len(include)-1):
                        sql_include += ' `%s` = 1 and ' % include[i]
                    sql_include+='`%s` = 1' % include[-1]
            else:
                sql_include = ''
                if len(include) > 0:
                    for i in range(0, len(include)-1):
                        sql_include += ' `%s` > 0 and ' % include[i]
                    sql_include+='`%s` > 0' % include[-1]

            sql_exclude = ''
            if len(exclude) > 0:
                sql_exclude = 'and '
                for i in range(0, len(exclude)-1):
                    sql_exclude += ' `%s` = 0 and ' % exclude[i]
                sql_exclude+='`%s` = 0' % exclude[-1]

            server, db = manipulate_biosqldb.load_db(biodb)

            sql ='select orthogroup from orthology_%s where %s %s' % (biodb, sql_include, sql_exclude)
            print sql


            match_groups = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]
            sum_group = len(match_groups)


            sql_include = 'taxon_id ='
            for i in range(0, len(include)-1):
                sql_include+='%s or taxon_id =' % include[i]
            sql_include+=include[-1]
            n = 1
            search_result = []

            group_filter = 'where orthogroup in ('

            for i, group in enumerate(match_groups[0:-1]):
                    group_filter += '"%s",' % group
            group_filter = group_filter[0:-1]+'"%s")' % match_groups[-1][1]


            columns = 'orthogroup, locus_tag, protein_id, start, stop, ' \
                      'strand, gene, orthogroup_size, n_genomes, TM, SP, product, organism, translation'
            sql_2 = 'select %s from orthology_detail_%s %s' % (columns, biodb, group_filter)

            raw_data = server.adaptor.execute_and_fetchall(sql_2,)
            import biosql_own_sql_tables
            orthogroup2genes = biosql_own_sql_tables.orthogroup2gene(biodb)
            orthogroup2products = biosql_own_sql_tables.orthogroup2product(biodb)
            orthogroup2cogs = biosql_own_sql_tables.orthogroup2cog(biodb, match_groups)
            orthogroup2pfam = biosql_own_sql_tables.orthogroup2pfam(biodb)

            match_groups_data = []
            group_data = ''

            for i, group in enumerate(match_groups):
                genes_data = ''
                for gene in orthogroup2genes[group]:
                    genes_data += '%s (%s)<br/>' % (gene, orthogroup2genes[group][gene])
                genes_data = genes_data[0:-5]
                product_data = ''
                for product in orthogroup2products[group]:
                    product_data += '%s (%s)<br/>' % (product, orthogroup2products[group][product])
                cog_data = ''
                for cog in orthogroup2cogs[group]:
                    if cog == None:
                        cog_data += '<p>- (%s)</p><br/>' % (orthogroup2cogs[group][cog])
                    else:
                        cog_data += '<a href="http://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid=%s">' \
                                    '%s (%s)</a><br/>' % (cog,cog, orthogroup2cogs[group][cog])
                cog_data = cog_data[0:-5]
                pfam_data = ''
                try:
                    for pfam in orthogroup2pfam[group]:

                        pfam_data += '<a href="http://pfam.xfam.org/family/%s">%s (%s)</a><br/>' % (pfam,pfam, orthogroup2pfam[group][pfam])
                    pfam_data = pfam_data[0:-5]
                except:
                    pfam_data += ' <p>-</p> '
                    pfam_data = pfam_data[0:-5]



                match_groups_data.append([i, group, genes_data, product_data, cog_data, pfam_data])

            n = 1
            extract_result = []
            for one_hit in raw_data:
                extract_result.append((n,) + one_hit)
                n+=1
                #print n
            #print extract_result

            envoi_extract = True


    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = extract_form_class()  # Nous créons un formulaire vide

    return render(request, 'chlamdb/extract_orthogroup.html', locals())




@login_required
def venn_orthogroup(request, biodb):

    cache = get_cache('default')
    print "loading db..."
    server = manipulate_biosqldb.load_db()
    print "db loaded..."
    venn_form_class = make_venn_from(biodb)
    if request.method == 'POST':  # S'il s'agit d'une requête POST

        form_venn = venn_form_class(request.POST)
        #form2 = ContactForm(request.POST)
        if 'venn' in request.POST and form_venn.is_valid():
            targets = form_venn.cleaned_data['targets']

            server, db = manipulate_biosqldb.load_db(biodb)

            all_orthogroups_list = []
            series = '['
            taxon_id2genome = manipulate_biosqldb.taxon_id2genome_description(server, biodb)
            for target in targets:
                template_serie = '{name: "%s", data: %s}'
                sql ='select orthogroup from orthology_%s where `%s` > 0' % (biodb, target)
                print sql
                orthogroups = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]
                all_orthogroups_list += orthogroups
                data = '"' + '","'.join(orthogroups) + '"'
                series+=template_serie % (taxon_id2genome[target], orthogroups) + ','
            series = series[0:-1] + ']'


            #h['Marilyn Monroe'] = 1;

            orthogroup2description = ''
            sql = 'select orthogroup, gene from orthology_detail_%s' % biodb
            data = server.adaptor.execute_and_fetchall(sql,)
            orthogroup2genes = {}
            for i in data:
                if i[0] not in orthogroup2genes:
                    orthogroup2genes[i[0]] = [i[1]]
                else:
                    if i[1] in orthogroup2genes[i[0]]:
                        pass
                    else:
                        orthogroup2genes[i[0]].append(i[1])
            sql = 'select orthogroup, product from orthology_detail_%s' % biodb
            data = server.adaptor.execute_and_fetchall(sql,)
            orthogroup2product = {}
            for i in data:
                if i[0] not in orthogroup2product:
                    orthogroup2product[i[0]] = [i[1]]
                else:
                    if i[1] in orthogroup2product[i[0]]:
                        pass
                    else:
                        orthogroup2product[i[0]].append(i[1])

            for i in orthogroup2genes:
                if i in all_orthogroups_list:
                    genes = '<br>'.join(orthogroup2genes[i])
                    products = '<br>'.join(orthogroup2product[i])
                    orthogroup2description+='h["%s"] = "%s</td><td>%s;"\n' % (i, genes, products)
                else:
                    continue
            print orthogroup2description
            #print series
            envoi_venn = True
    else:  # Si ce n'est pas du POST, c'est probablement une requête GET  # Nous créons un formulaire vide
        form_venn = venn_form_class()
    return render(request, 'chlamdb/venn_orthogroup.html', locals())


@login_required
def extract_pfam(request, biodb):

    cache = get_cache('default')
    print "loading db..."
    server = manipulate_biosqldb.load_db()
    print "db loaded..."
    extract_form_class = make_extract_form(biodb)

    if request.method == 'POST':  # S'il s'agit d'une requête POST

        form = extract_form_class(request.POST)  # Nous reprenons les données

        #form2 = ContactForm(request.POST)
        if 'comparison' in request.POST and form.is_valid():  # Nous vérifions que les données envoyées sont valides

            print request.POST
            print form.cleaned_data.keys()
            include = form.cleaned_data['orthologs_in']
            exclude = form.cleaned_data['no_orthologs_in']

            try:
                single_copy = request.POST['button_single_copy']
                single_copy = True
            except:
                single_copy = False
            print "single_copy", single_copy


            if single_copy:
                sql_include = ''
                if len(include) > 0:
                    for i in range(0, len(include)-1):
                        sql_include += ' `%s` = 1 and ' % include[i]
                    sql_include+='`%s` = 1' % include[-1]
            else:
                sql_include = ''
                if len(include) > 0:
                    for i in range(0, len(include)-1):
                        sql_include += ' `%s` > 0 and ' % include[i]
                    sql_include+='`%s` > 0' % include[-1]

            sql_exclude = ''
            if len(exclude) > 0:
                sql_exclude = 'and '
                for i in range(0, len(exclude)-1):
                    sql_exclude += ' `%s` = 0 and ' % exclude[i]
                sql_exclude+='`%s` = 0' % exclude[-1]

            server, db = manipulate_biosqldb.load_db(biodb)

            sql ='select id from comparative_tables.Pfam_%s where %s %s' % (biodb, sql_include, sql_exclude)

            print sql

            match_groups = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]
            sum_group = len(match_groups)


            sql_include = 'taxon_id ='
            for i in range(0, len(include)-1):
                sql_include+='%s or taxon_id =' % include[i]
            sql_include+=include[-1]
            n = 1
            search_result = []

            group_filter = 'where ('

            for i, group in enumerate(match_groups):
                if i == 0:
                    group_filter += 'orthogroup="%s"' % group
                else:
                    group_filter += ' or orthogroup="%s"' % group
            group_filter+=')'
            #print group_filter


            columns = 'orthogroup, locus_tag, protein_id, start, stop, ' \
                      'strand, gene, orthogroup_size, n_genomes, TM, SP, product, organism, translation'
            sql_2 = 'select %s from orthology_detail_%s %s' % (columns, biodb, group_filter)
            #print sql_2
            raw_data = server.adaptor.execute_and_fetchall(sql_2,)


            import biosql_own_sql_tables
            pfam2descr = biosql_own_sql_tables.pfam2description(biodb)
            match_groups_data = []
            for i, pfam in enumerate(match_groups):
                match_groups_data.append([i, pfam, pfam2descr[pfam]])



            n = 1
            extract_result = []
            for one_hit in raw_data:
                extract_result.append((n,) + one_hit)
                n+=1
                #print n
            #print extract_result

            envoi_extract = True


    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = extract_form_class()  # Nous créons un formulaire vide

    return render(request, 'chlamdb/extract_Pfam.html', locals())




@login_required
def venn_pfam(request, biodb):

    cache = get_cache('default')
    print "loading db..."
    server = manipulate_biosqldb.load_db()
    print "db loaded..."
    venn_form_class = make_venn_from(biodb)
    if request.method == 'POST':  # S'il s'agit d'une requête POST

        form_venn = venn_form_class(request.POST)
        #form2 = ContactForm(request.POST)
        if 'venn' in request.POST and form_venn.is_valid():
            targets = form_venn.cleaned_data['targets']

            server, db = manipulate_biosqldb.load_db(biodb)

            all_pfam_list = []
            series = '['
            taxon_id2genome = manipulate_biosqldb.taxon_id2genome_description(server, biodb)
            for target in targets:
                template_serie = '{name: "%s", data: %s}'
                sql ='select id from comparative_tables.Pfam_%s where `%s` > 0' % (biodb, target)
                print sql
                cogs = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]
                all_pfam_list += cogs
                data = '"' + '","'.join(cogs) + '"'
                series+=template_serie % (taxon_id2genome[target], cogs) + ','
            series = series[0:-1] + ']'


            pfam2description = ''
            sql = 'select signature_accession, signature_description, count(*) from interpro_%s where analysis="Pfam" group by signature_accession;' % biodb
            data = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
            for i in data:
                if i in all_pfam_list:
                    print 'ok'
                    pfam2description+='h["%s"] = "%s </td><td>%s";' % (i, data[i][0], data[i][1])
                else:
                    print 'pas ok'
            print pfam2description
            #print series

            envoi_venn = True
    else:  # Si ce n'est pas du POST, c'est probablement une requête GET  # Nous créons un formulaire vide
        form_venn = venn_form_class()
    return render(request, 'chlamdb/venn_Pfam.html', locals())

@login_required
def extract_interpro(request, biodb):

    cache = get_cache('default')
    print "loading db..."
    server = manipulate_biosqldb.load_db()
    print "db loaded..."
    extract_form_class = make_extract_form(biodb)

    if request.method == 'POST':  # S'il s'agit d'une requête POST

        form = extract_form_class(request.POST)  # Nous reprenons les données

        #form2 = ContactForm(request.POST)
        if 'comparison' in request.POST and form.is_valid():  # Nous vérifions que les données envoyées sont valides

            print request.POST
            print form.cleaned_data.keys()
            include = form.cleaned_data['orthologs_in']
            exclude = form.cleaned_data['no_orthologs_in']

            try:
                single_copy = request.POST['button_single_copy']
                single_copy = True
            except:
                single_copy = False
            print "single_copy", single_copy


            if single_copy:
                sql_include = ''
                if len(include) > 0:
                    for i in range(0, len(include)-1):
                        sql_include += ' `%s` = 1 and ' % include[i]
                    sql_include+='`%s` = 1' % include[-1]
            else:
                sql_include = ''
                if len(include) > 0:
                    for i in range(0, len(include)-1):
                        sql_include += ' `%s` > 0 and ' % include[i]
                    sql_include+='`%s` > 0' % include[-1]

            sql_exclude = ''
            if len(exclude) > 0:
                sql_exclude = 'and '
                for i in range(0, len(exclude)-1):
                    sql_exclude += ' `%s` = 0 and ' % exclude[i]
                sql_exclude+='`%s` = 0' % exclude[-1]

            server, db = manipulate_biosqldb.load_db(biodb)

            sql ='select id from comparative_tables.interpro_%s where %s %s' % (biodb, sql_include, sql_exclude)

            print sql

            match_groups = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]
            sum_group = len(match_groups)

            filter = '"' + '","'.join(match_groups) + '"'

            sql2 = 'select interpro_accession, interpro_description from interpro_%s' \
            ' where interpro_accession in (%s) group by interpro_accession;' % (biodb, filter)

            match_data = list(server.adaptor.execute_and_fetchall(sql2,))


            sql_include = 'taxon_id ='
            for i in range(0, len(include)-1):
                sql_include+='%s or taxon_id =' % include[i]
            sql_include+=include[-1]
            n = 1
            search_result = []

            group_filter = 'where ('

            for i, group in enumerate(match_groups):
                if i == 0:
                    group_filter += 'orthogroup="%s"' % group
                else:
                    group_filter += ' or orthogroup="%s"' % group
            group_filter+=')'
            #print group_filter


            columns = 'orthogroup, locus_tag, protein_id, start, stop, ' \
                      'strand, gene, orthogroup_size, n_genomes, TM, SP, product, organism, translation'
            sql_2 = 'select %s from orthology_detail_%s %s' % (columns, biodb, group_filter)
            #print sql_2
            raw_data = server.adaptor.execute_and_fetchall(sql_2,)

            n = 1
            extract_result = []
            for one_hit in raw_data:
                extract_result.append((n,) + one_hit)
                n+=1
                #print n
            #print extract_result

            envoi_extract = True


    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = extract_form_class()  # Nous créons un formulaire vide

    return render(request, 'chlamdb/extract_interpro.html', locals())




@login_required
def venn_interpro(request, biodb):

    cache = get_cache('default')
    print "loading db..."
    server = manipulate_biosqldb.load_db()
    print "db loaded..."
    venn_form_class = make_venn_from(biodb)
    if request.method == 'POST':  # S'il s'agit d'une requête POST

        form_venn = venn_form_class(request.POST)
        #form2 = ContactForm(request.POST)
        if 'venn' in request.POST and form_venn.is_valid():
            targets = form_venn.cleaned_data['targets']

            server, db = manipulate_biosqldb.load_db(biodb)

            all_pfam_list = []
            series = '['
            taxon_id2genome = manipulate_biosqldb.taxon_id2genome_description(server, biodb)
            for target in targets:
                template_serie = '{name: "%s", data: %s}'
                sql ='select id from comparative_tables.interpro_%s where `%s` > 0' % (biodb, target)
                print sql
                cogs = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]
                all_pfam_list += cogs
                data = '"' + '","'.join(cogs) + '"'
                series+=template_serie % (taxon_id2genome[target], cogs) + ','
            series = series[0:-1] + ']'


            interpro2description = ''
            sql = 'select interpro_accession,interpro_description, count(*) from interpro_%s group by interpro_accession;' % biodb
            data = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
            for i in data:
                if i in all_pfam_list:
                    print 'ok'
                    interpro2description+='h["%s"] = "%s </td><td>%s";' % (i, data[i][0], data[i][1])
                else:
                    print 'pas ok'
            print interpro2description
            #print series

            envoi_venn = True
    else:  # Si ce n'est pas du POST, c'est probablement une requête GET  # Nous créons un formulaire vide
        form_venn = venn_form_class()
    return render(request, 'chlamdb/venn_interpro.html', locals())


@login_required
def extract_cog(request, biodb):

    cache = get_cache('default')
    print "loading db..."
    server = manipulate_biosqldb.load_db()
    print "db loaded..."
    extract_form_class = make_extract_form(biodb)


    if request.method == 'POST':  # S'il s'agit d'une requête POST

        form = extract_form_class(request.POST)  # Nous reprenons les données
        #form2 = ContactForm(request.POST)
        if form.is_valid():  # Nous vérifions que les données envoyées sont valides

            include = form.cleaned_data['orthologs_in']
            exclude = form.cleaned_data['no_orthologs_in']
            try:
                single_copy = request.POST['button_single_copy']
                single_copy = True
            except:
                single_copy = False
            print "single_copy", single_copy

            print "single_copy", single_copy


            if single_copy:
                sql_include = ''
                if len(include) > 0:
                    for i in range(0, len(include)-1):
                        sql_include += ' `%s` = 1 and ' % include[i]
                    sql_include+='`%s` = 1' % include[-1]
            else:
                sql_include = ''
                if len(include) > 0:
                    for i in range(0, len(include)-1):
                        sql_include += ' `%s` > 0 and ' % include[i]
                    sql_include+='`%s` > 0' % include[-1]

            sql_exclude = ''
            if len(exclude) > 0:
                sql_exclude = 'and '
                for i in range(0, len(exclude)-1):
                    sql_exclude += ' `%s` = 0 and ' % exclude[i]
                sql_exclude+='`%s` = 0' % exclude[-1]

            server, db = manipulate_biosqldb.load_db(biodb)

            sql ='select id from comparative_tables.COG_%s where %s %s' % (biodb, sql_include, sql_exclude)

            print sql

            match_groups = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]
            sum_group = len(match_groups)

            cog_data = []
            for i in match_groups:
                sql = 'select * from COG.cog_names_2014 where COG_id ="%s"' % i
                cog_data.append(list(server.adaptor.execute_and_fetchall(sql,)[0]))
            print cog_data

            '''
            sql_include = 'taxon_id ='
            for i in range(0, len(include)-1):
                sql_include+='%s or taxon_id =' % include[i]
            sql_include+=include[-1]
            n = 1
            search_result = []

            group_filter = 'where ('

            for i, group in enumerate(match_groups):
                if i == 0:
                    group_filter += 'orthogroup="%s"' % group
                else:
                    group_filter += ' or orthogroup="%s"' % group
            group_filter+=')'
            #print group_filter


            columns = 'orthogroup, locus_tag, protein_id, start, stop, ' \
                      'strand, gene, orthogroup_size, n_genomes, TM, SP, product, organism, translation'
            sql_2 = 'select %s from orthology_detail_%s %s' % (columns, biodb, group_filter)
            #print sql_2
            raw_data = server.adaptor.execute_and_fetchall(sql_2,)

            n = 1
            extract_result = []
            for one_hit in raw_data:
                extract_result.append((n,) + one_hit)
                n+=1
                #print n
            #print extract_result
            '''
            envoi_extract = True

    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = extract_form_class()  # Nous créons un formulaire vide

    return render(request, 'chlamdb/extract_cogs.html', locals())



@login_required
def venn_cog(request, biodb):

    cache = get_cache('default')
    print "loading db..."
    server = manipulate_biosqldb.load_db()
    print "db loaded..."
    venn_form_class = make_venn_from(biodb)
    if request.method == 'POST':  # S'il s'agit d'une requête POST

        form_venn = venn_form_class(request.POST)  # Nous reprenons les données
        #form2 = ContactForm(request.POST)
        if form_venn.is_valid():  # Nous vérifions que les données envoyées sont valides

            targets = form_venn.cleaned_data['targets']

            server, db = manipulate_biosqldb.load_db(biodb)

            all_cog_list = []
            series = '['
            taxon_id2genome = manipulate_biosqldb.taxon_id2genome_description(server, biodb)
            for target in targets:
                template_serie = '{name: "%s", data: %s}'
                sql ='select id from comparative_tables.COG_%s where `%s` > 0' % (biodb, target)
                print sql
                cogs = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]
                all_cog_list += cogs
                data = '"' + '","'.join(cogs) + '"'
                series+=template_serie % (taxon_id2genome[target], cogs) + ','
            series = series[0:-1] + ']'


            #h['Marilyn Monroe'] = 1;

            cog2description = ''
            sql = 'select * from COG.cog_names_2014'
            data = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
            for i in data:
                if i in all_cog_list:
                    print 'ok'
                    cog2description+='h["%s"] = "%s </td><td>%s";' % (i, data[i][0], data[i][1])
                else:
                    print 'pas ok'
            print cog2description
            #print series
            envoi_venn = True

    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form_venn = venn_form_class()  # Nous créons un formulaire vide

    return render(request, 'chlamdb/venn_cogs.html', locals())



@login_required
def extract_region(request, biodb):
    extract_region_form_class = make_extract_region_form(biodb)
    server, db = manipulate_biosqldb.load_db(biodb)
    from Bio.Alphabet import generic_protein
    import re

    cache = get_cache('default')

    if request.method == 'POST':

        form = extract_region_form_class(request.POST)

        if form.is_valid():
            genome_accession = form.cleaned_data['genome']

            region = form.cleaned_data['region']

            start_stop = re.sub(' ', '', form.cleaned_data['region']).split(",")

            extract = form.cleaned_data['extract']

            genome_description = manipulate_biosqldb.accession2description(server, biodb)[genome_accession]


            if extract == 'annotation':
                get_annotation = True
                columns = 'orthogroup, locus_tag, protein_id, start, stop, ' \
                          'strand, gene, orthogroup_size, n_genomes, TM, SP, product, organism, translation'


                sql2 = 'select %s from orthology_detail_%s where start > %s and stop < %s and accession="%s"' % (columns,
                                                                                                                 biodb,
                                                                                                                 start_stop[0],
                                                                                                                 start_stop[1],
                                                                                                                 genome_accession)


                try:
                    raw_data = server.adaptor.execute_and_fetchall(sql2, )
                except IndexError:
                    valid_id = False
                    return render(request, 'chlamdb/extract_region.html', locals())
                if not raw_data:
                        valid_id = False
                else:
                    n = 1
                    search_result = []
                    for one_hit in raw_data:

                        sql = 'select contig from '


                        if one_hit[2] != '-':
                            interpro_id = one_hit[2]
                        else:
                            interpro_id = one_hit[1]
                        search_result.append((n,) + one_hit + (interpro_id,))
                        n+=1


            if extract == 'sequence' or extract == 'sequence_trans':
                get_sequence = True
                seq = manipulate_biosqldb.location2sequence(server, genome_accession, biodb, start_stop[0], int(start_stop[1])-int(start_stop[0]))

                record = SeqRecord(Seq(seq, generic_protein),
                   id="%s_%s_%s" % (genome_accession, start_stop[0], start_stop[1]), description=genome_description)

                temp_file = os.path.join(settings.BASE_DIR, "assets/temp/%s_region.fa" % genome_accession)
                temp_location = "temp/%s_region.fa" % genome_accession
                print "temp_file", temp_file
                with open(temp_file, 'w') as f:
                    SeqIO.write(record, f, 'fasta')


                if extract == 'sequence_trans':
                    extract_trans = True
                    from Bio.Seq import reverse_complement, translate
                    anti = reverse_complement(seq)
                    comp = anti[::-1]
                    length = len(seq)
                    frames = {}
                    for i in range(0, 3):
                        fragment_length = 3 * ((length-i) // 3)
                        frames[i+1] = translate(seq[i:i+fragment_length], 1)
                        frames[-(i+1)] = translate(anti[i:i+fragment_length], 1)[::-1]



                    frame_plus_1 = frames[1]
                    frame_plus_2 = frames[2]
                    frame_plus_3 = frames[3]
                    frame_minus_1 = frames[-1]
                    frame_minus_2 = frames[-2]
                    frame_minus_3 = frames[-3]



            envoi = True

    else:
        form = extract_region_form_class()
    return render(request, 'chlamdb/extract_region.html', locals())




@login_required
def locusx(request, biodb, locus=None, menu=False):

    print 'biodb', biodb
    print 'locus', locus
    print 'menu', menu

    cache = get_cache('default')

    #cache.clear()

    if request.method == 'GET':  # S'il s'agit d'une requête POST

        valid_id = True

        server, db = manipulate_biosqldb.load_db(biodb)

        #sql1 = 'SELECT column_name FROM information_schema.columns WHERE table_name="orthology_detail_chlamydia_03_15"'

        sql1 =   'SELECT' \
                 ' CASE' \
                 '   WHEN locus_tag = "%s" THEN "locus_tag"' \
                 '   WHEN protein_id = "%s" THEN "protein_id"' \
                 '   WHEN orthogroup = "%s" THEN "orthogroup"'\
                 ' END AS "which_column"'\
                 ' FROM' \
                 ' orthology_detail_%s where locus_tag="%s" or protein_id like "%%%%%s%%%%" or orthogroup="%s"' % (locus,
                                                                                                                   locus,
                                                                                                                   locus,
                                                                                                                   biodb,
                                                                                                                   locus,
                                                                                                                   locus,
                                                                                                                   locus)



        try:
            print sql1
            input_type = server.adaptor.execute_and_fetchall(sql1, )[0][0]
        except IndexError:
            print 'not a valid id'
            valid_id = False
            return render(request, 'chlamdb/locus.html', locals())

        else:
            columns = 'orthogroup, locus_tag, protein_id, start, stop, ' \
                      'strand, gene, orthogroup_size, n_genomes, TM, SP, product, organism, translation'
            sql2 = 'select %s from orthology_detail_%s where %s="%s"' % (columns, biodb, input_type, locus)
            if input_type == 'locus_tag':
                sql3 = 'select t2.COG_id,t2.functon,t2.name from COG.locus_tag2gi_hit_%s ' \
                       ' as t1 inner join COG.cog_names_2014 as t2 on t1.COG_id=t2.COG_id where locus_tag="%s"' % (biodb, locus)
                sql4 = 'select analysis, signature_accession, signature_description, interpro_accession, interpro_description ' \
                       ' from interpro_%s where locus_tag="%s";' % (biodb, locus)
                try:
                    cog_data = server.adaptor.execute_and_fetchall(sql3, )[0]
                except IndexError:
                    cog_data = False
                try:
                    interpro_data = server.adaptor.execute_and_fetchall(sql4, )
                except IndexError:
                    interpro_data= False

            data = list(server.adaptor.execute_and_fetchall(sql2, )[0])
            if data[2] == '-':
                data[2] = data[1]

            orthogroup = data[0]

            fasta = "%s_fasta/%s.txt" % (biodb, orthogroup)
            alignment = "%s_fasta/%s.html" % (biodb, orthogroup)
            alignment_fasta = "%s_fasta/%s.fa" % (biodb, orthogroup)
            alignment_fasta_nucl = "%s_fasta_nucl/%s_nucl.txt" % (biodb, orthogroup)
            tree_unrooted = "%s_fasta/%s_tree.svg" % (biodb, orthogroup)
            tree_rooted = "%s_fasta/%s_tree_reroot.svg" % (biodb, orthogroup)
            tree_file = "%s_fasta/%s.phy_phyml_tree.txt" % (biodb, orthogroup)

            columns = 'orthogroup, locus_tag, protein_id, start, stop, ' \
                  'strand, gene, orthogroup_size, n_genomes, TM, SP, product, organism, translation'
            sql3 = 'select %s from orthology_detail_%s where orthogroup = "%s" ' % (columns, biodb, orthogroup)

            homologues = list(server.adaptor.execute_and_fetchall(sql3, ))
            print homologues

            if len(homologues) >1:
                orthologs = True
            else:
                orthologs = False
            import orthogroup_identity_db
            if len(homologues) > 1:
                orthogroup2identity_dico = orthogroup_identity_db.orthogroup2identity_dico(biodb, orthogroup)

                print "orthologs", orthologs, len(homologues)
                for count, value in enumerate(homologues):
                    value = list(value)
                    locus_2 = value[1]
                    if value[2] != '-':
                        interpro_id = value[2]
                    else:
                        value[2] = value[1]
                    #print value + [orthogroup2identity_dico[data[1]][locus_2]]
                    homologues[count] = [count+1] + value + [orthogroup2identity_dico[data[1]][locus_2]]
                    #print homologues[count]

            else:
                homologues[0] = (1,) + homologues[0] + (100,)

        envoi = True


    return render(request, 'chlamdb/locus.html', locals())


@login_required
def fam(request, biodb, fam, type):

    cache = get_cache('default')

    #cache.clear()

    if request.method == 'GET':  # S'il s'agit d'une requête POST

        valid_id = True

        server, db = manipulate_biosqldb.load_db(biodb)

        print 'type', type

        #sql1 = 'SELECT column_name FROM information_schema.columns WHERE table_name="orthology_detail_chlamydia_03_15"'
        if type =='pfam':
            sql1 =   'select locus_tag, signature_description, count(*) from interpro_%s where signature_accession="%s" group by locus_tag' % (biodb, fam)
        elif type == 'cog':
            sql1 = 'select locus_tag from COG.locus_tag2gi_hit_%s where COG_id="%s"' % (biodb, fam)
            sql2 = 'select COG_id,functon, name from COG.cog_names_2014 where COG_id = "%s"' % (fam)
            print 'cog type', sql2
            info = server.adaptor.execute_and_fetchall(sql2, )[0]



        elif type == 'interpro':
            sql1 = 'select locus_tag, signature_description, count(*) from interpro_%s where interpro_accession="%s" group by locus_tag' % (biodb, fam)
        else:
            valid_id = False
            return render(request, 'chlamdb/fam.html', locals())

        try:
            locus_list = [i[0] for i in server.adaptor.execute_and_fetchall(sql1, )]
            print 'locus', locus_list
            locus_list_form = '"' + '","'.join(locus_list) + '"'
        except IndexError:
            valid_id = False
            return render(request, 'chlamdb/fam.html', locals())

        else:
            columns = 'orthogroup, locus_tag, protein_id, start, stop, ' \
                      'strand, gene, orthogroup_size, n_genomes, TM, SP, product, organism, translation'
            sql2 = 'select %s from orthology_detail_%s where locus_tag in (%s)' % (columns, biodb, locus_list_form)

            all_locus_raw_data = server.adaptor.execute_and_fetchall(sql2, )
            all_locus_data = []
            group_count = []
            for i in range(0, len(all_locus_raw_data)):
                 all_locus_data.append([i] + list(all_locus_raw_data[i]))
                 if all_locus_raw_data[i][0] not in group_count:
                    group_count.append(all_locus_raw_data[i][0])
        envoi = True
        menu = True


    return render(request, 'chlamdb/fam.html', locals())


@login_required
def sunburst(request, biodb, locus):

    cache = get_cache('default')

    #cache.clear()

    if request.method == 'GET':  # S'il s'agit d'une requête POST
        print 'okkkkkkkkkkkkkkk'
        valid_id = True

        server, db = manipulate_biosqldb.load_db(biodb)

        sql = 'select accession from orthology_detail_%s where locus_tag = "%s"' % (biodb, locus)
        accession = server.adaptor.execute_and_fetchall(sql,)[0][0]
        print accession
        try:
            sql1 = 'select t3.superkingdom,  t3.phylum,  t3.order,  t3.family,  t3.genus,  t3.species  from ' \
                   ' blastnr.blastnr_hits_%s_%s as t1  inner join blastnr.blastnr_hits_taxonomy_filtered_%s_%s ' \
                   ' as t2 on t1.nr_hit_id = t2.nr_hit_id  inner join blastnr.blastnr_taxonomy as t3 on ' \
                   ' t2.subject_taxon_id = t3.taxon_id inner join blastnr.blastnr_hsps_chlamydia_03_15_%s as t4 ' \
                   ' on t1.nr_hit_id=t4.nr_hit_id where t1.locus_tag="%s"' % (biodb, accession, biodb, accession, accession, locus)
            print sql
            raw_data = server.adaptor.execute_and_fetchall(sql1,)

        except:
            valid_id = False
            return render(request, 'chlamdb/sunburst.html', locals())


        dico =  {}

        y = 1
        for data in raw_data:

            if data[0] not in dico:
                dico[data[0]] = {}
            else:
                if data[1] not in dico[data[0]]:
                    dico[data[0]][data[1]] = {}

                if data[2] not in dico[data[0]][data[1]]:
                    dico[data[0]][data[1]][data[2]] = {}

                if data[3] not in dico[data[0]][data[1]][data[2]]:
                    dico[data[0]][data[1]][data[2]][data[3]] = {}

                if data[4] not in dico[data[0]][data[1]][data[2]][data[3]]:
                    dico[data[0]][data[1]][data[2]][data[3]][data[4]] = {}

                if data[5] not in dico[data[0]][data[1]][data[2]][data[3]][data[4]]:
                    dico[data[0]][data[1]][data[2]][data[3]][data[4]][data[5]] = 1

                dico[data[0]][data[1]][data[2]][data[3]][data[4]][data[5]] += 1

                y += 1
        i = 1

        tt = settings.BASE_DIR + '/assets/out.tab'
        print 'sdsdfsdf', tt
        out = open(tt, 'w')
        for superkingdom in dico:
            for phylum in dico[superkingdom]:
                for order in dico[superkingdom][phylum]:
                    for family in dico[superkingdom][phylum][order]:
                        for genus in dico[superkingdom][phylum][order][family]:
                            for species in dico[superkingdom][phylum][order][family][genus]:
                                out.write('%s, %s, %s, %s\n' % (i, 1, superkingdom, 0))
                                out.write("%s, %s, %s, %s\n" % (i, 2, phylum, 0))
                                out.write("%s, %s, %s, %s\n" % (i, 3, order, 0))
                                out.write("%s, %s, %s, %s\n" % (i, 4, family, 0))
                                out.write("%s, %s, %s, %s\n" % (i, 5, genus, 0))
                                out.write("%s, %s, %s, %s\n" % (i, 6, species, dico[superkingdom][phylum][order][family][genus][species]))
                                i+=1

                                print '\'%s, %s, %s, %s\\n\' +' % (i, 1, superkingdom[1:-1], 0)
                                print "\'%s, %s, %s, %s\\n\' +" % (i, 2, phylum[1:-1], 0)
                                print "\'%s, %s, %s, %s\\n\' +" % (i, 3, order[1:-1], 0)
                                print "\'%s, %s, %s, %s\\n\' +" % (i, 4, family[1:-1], 0)
                                print "\'%s, %s, %s, %s\\n\' +" % (i, 5, genus[1:-1], 0)
                                print "\'%s, %s, %s, %s\\n\' +" % (i, 6, species[1:-1], dico[superkingdom][phylum][order][family][genus][species])
        out.close()

        envoi = True
        menu = True


    return render(request, 'chlamdb/sunburst.html', locals())



@login_required
def blastnr(request, biodb, locus_tag):


    print biodb, locus_tag

    cache = get_cache('default')
    print "cache", cache
    #cache.clear()

    #bioentry_in_memory = cache.get("biodb")
    print "loading db..."
    server = manipulate_biosqldb.load_db()
    print "db loaded..."
    if request.method == 'GET':  # S'il s'agit d'une requête POST



        server, db = manipulate_biosqldb.load_db(biodb)

        sql = 'select accession, organism from orthology_detail_%s where locus_tag="%s"' % (biodb, locus_tag)
        data = server.adaptor.execute_and_fetchall(sql,)[0]
        accession = data[0]
        organism = data[1]

        server, db = manipulate_biosqldb.load_db(biodb)

        columns = 't1.hit_number, t1.locus_tag, t1.query_accession, t1.subject_accession, t1.subject_taxid, t1.subject_scientific_name,' \
                  't1.subject_title, t1.subject_kingdom, t2.evalue, t2.percent_identity, t2.positive, t2.gaps, t2.length,' \
                  't2.query_start, t2.query_end, t2.query_cov, t2.subject_start, t2.subject_end, t1.subject_title'
        sql = 'select %s from blastnr.blastnr_hits_%s_%s as t1 inner join blastnr.blastnr_hsps_%s_%s as t2 on' \
              ' t1.nr_hit_id = t2.nr_hit_id where  t1.locus_tag="%s" ' % (columns, biodb, accession, biodb, accession, locus_tag)
        print sql
        blast_data = list(server.adaptor.execute_and_fetchall(sql))
        #blast_data = [i for i in ]


        if len(blast_data) > 0:
            valid_id = True
            blast_query_locus = blast_data[0][1]
            blast_query_protein_id = blast_data[0][2]
            if blast_query_protein_id == blast_query_locus:
                 blast_query_protein_id = ''

            for n, one_hit in enumerate(blast_data):
                blast_data[n] = [i for i in one_hit]
                subject_taxids = one_hit[4].split(';')
                subject_scientific_names = one_hit[5].split(';')
                all_taxonomy = ''
                for taxon, name in zip(subject_taxids, subject_scientific_names):
                    all_taxonomy += '<a href="http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=%s">%s<a> ' % (taxon, name)
                blast_data[n][5] = all_taxonomy



        return render(request, 'chlamdb/blastnr.html', locals())


    return render(request, 'chlamdb/blastnr.html', locals())


@login_required
def homology(request, biodb):
    import shell_command
    cache = get_cache('default')
    print "cache", cache
    #cache.clear()

    #bioentry_in_memory = cache.get("biodb")
    print "loading db..."
    server = manipulate_biosqldb.load_db()
    print "db loaded..."
    contact_form_class = make_contact_form(server, biodb)
    if request.method == 'POST':  # S'il s'agit d'une requête POST


        form = contact_form_class(request.POST)  # Nous reprenons les données
        #form2 = ContactForm(request.POST)
        if form.is_valid():  # Nous vérifions que les données envoyées sont valides
            accession = request.POST['accession']

            envoi = True
    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = make_contact_form(server, biodb)  # Nous créons un formulaire vide
    return render(request, 'chlamdb/homology.html', locals())


@login_required
def orthogroup_identity(request, biodb, orthogroup, group=False):


    print "loading db..."
    server = manipulate_biosqldb.load_db()
    print "db loaded..."
    #if request.method == 'POST':
    import numpy
    import pandas as pd
    sql = 'SELECT * FROM orth_%s.%s;' % (biodb, orthogroup)

    try:
        data = numpy.array([list(i) for i in server.adaptor.execute_and_fetchall(sql,)])
        homologs = True
    except:
        homologs = False
        return render(request, 'chlamdb/orthogroup_identity.html', locals())
    locus_list = '"' + '","'.join(data[0:,1]) + '"'

    sql2 = 'select locus_tag, organism from orthology_detail_%s where locus_tag in (%s)' % (biodb, locus_list)
    print sql2
    locus2organism = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2,))

    print locus2organism

    #print data

    print data[0:,2:]

    columns = data[0:,1]
    rows = [i + " (%s)" % locus2organism[i] for i in columns]
    frame = pd.DataFrame(data[0:,2:], index=rows, columns=columns)
    frame = frame.astype(float)
    frame = frame/100

    for i in range(0, len(frame)):

        frame.ix[i, i] = None

    path = settings.BASE_DIR + '/assets/temp/%s.json' % orthogroup
    print 'writing to', path
    print frame
    with open(path, 'w') as f:
        f.write(frame.to_json(orient="split"))

    return render(request, 'chlamdb/orthogroup_identity.html', locals())


def ortho_id_plot(request, group):
    return render(request, 'chlamdb/orthogroup_identity_plot.html', locals())


@login_required
def plot_region(request, biodb):

    cache = get_cache('default')
    print "cache", cache

    #bioentry_in_memory = cache.get("biodb")
    print "loading db..."
    server = manipulate_biosqldb.load_db()
    print "db loaded..."
    plot_region_form_class = make_plot_form(biodb)
    if request.method == 'POST':  # S'il s'agit d'une requête POST

        form = plot_region_form_class(request.POST)  # Nous reprenons les données
        #form2 = ContactForm(request.POST)
        if form.is_valid():  # Nous vérifions que les données envoyées sont valides
            valid_id = True

            accession = extract_alphanumeric(form.cleaned_data['accession'])

            server, db = manipulate_biosqldb.load_db(biodb)

            region_size = form.cleaned_data['region_size']

            genomes = form.cleaned_data['genomes']



            sql2 = 'select orthogroup, locus_tag, protein_id, start, stop, strand, organism from orthology_detail_%s where locus_tag like "%%%%%s%%%%" or protein_id like "%%%%%s%%%%"' % (biodb, accession, accession)
            print sql2
            data = server.adaptor.execute_and_fetchall(sql2, )[0]
            print "seqfeature_id", data
            if not data:
                    valid_id = False
            if valid_id:
                orthogroup = data[0]

                select = 'and (taxon_id = %s' % genomes[0]
                if len(genomes) >1:
                    for i in range(0, len(genomes)-1):
                        select+= ' or taxon_id = %s' % genomes[i]
                    select+= ' or taxon_id = %s)' % genomes[-1]
                else:
                    select+= ')'
                sql3 = 'select locus_tag from orthology_detail_%s where orthogroup = "%s" %s' % (biodb, orthogroup, select)
                print sql3
                locus_tag_target_genomes = [i[0] for i in server.adaptor.execute_and_fetchall(sql3, )]

                if plot_region:
                    print "plotting!!!!!!!!!!!!!!"
                    print "locus_tag_list", locus_tag_target_genomes
                    home_dir = os.path.dirname(__file__)
                    print "home_dir", home_dir
                    temp_location = os.path.join(home_dir, "../assets")
                    print "temp loc", temp_location
                    temp_file = NamedTemporaryFile(delete=False, dir=temp_location, suffix=".svg")
                    print "temp file", temp_file.name
                    name = os.path.basename(temp_file.name)
                    print name.split('.')
                    name_png = name.split('.')[0] + '.png'

                    locus_tags = mysqldb_plot_genomic_feature.proteins_id2cossplot(server, db, biodb, locus_tag_target_genomes,
                                                                                      temp_file.name, int(region_size),
                                                                                      cache)

                    columns = 'orthogroup, locus_tag, protein_id, start, stop, ' \
                              'strand, gene, orthogroup_size, n_genomes, TM, SP, product, organism, translation'


                    sql_locus = 'locus_tag="%s"' % locus_tags[0]
                    for locus in range(1, len(locus_tags)):
                        sql_locus += ' or locus_tag="%s"' % locus_tags[locus]

                    sql = 'select %s from orthology_detail_%s where %s' % (columns, biodb, sql_locus)
                    print sql

                    raw_data = server.adaptor.execute_and_fetchall(sql,)

                    n = 1
                    search_result = []
                    for one_hit in raw_data:
                        search_result.append((n,) + one_hit)
                        n+=1
                        print n
                    print search_result


            envoi = True

    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = plot_region_form_class()  # Nous créons un formulaire vide

    return render(request, 'chlamdb/plot_region.html', locals())

'''
@login_required
def plot_region(request, biodb):
    plot_form_class = make_plot_form(biodb)
    cache = get_cache('default')
    print "cache", cache
    server = manipulate_biosqldb.load_db()

    if request.method == 'POST':  # S'il s'agit d'une requête POST


        form = plot_form_class(request.POST)  # Nous reprenons les données

        if form.is_valid():  # Nous vérifions que les données envoyées sont valides






            envoi = True

    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = plot_form_class()  # Nous créons un formulaire vide

    return render(request, 'chlamdb/plot_region.html', locals())
'''


@login_required
def comparative_extract():
    '''
    get features present in genomes X Y Z and not in A B C

    Se baser sur la table d'orthologie pour identifier les genes

    SQL
    select orthogroup from orthology_chlam where `taxon X` > 1 and `taxon B` = 0

    :return:
    '''


@login_required
def orthogroups(request):
    if request.method == 'POST':
        form = BiodatabaseForm(request.POST)
        if form.is_valid():
            biodb = form.cleaned_data['biodatabase']
            groups = "orthogroup_size_distrib_%s.svg" % biodb
            envoi = True
    else:
        form = BiodatabaseForm()

    return render(request, 'chlamdb/orthogroups.html', locals())


@login_required
def circos(request, biodb):

    import gbk2circos
    circos_form_class = make_circos_form(biodb)
    server, db = manipulate_biosqldb.load_db(biodb)

    cache = get_cache('default')

    if request.method == 'POST':

        form = circos_form_class(request.POST)

        if form.is_valid():
            reference_taxon = form.cleaned_data['reference']

            description2accession_dict = manipulate_biosqldb.description2accession_dict(server, biodb)

            reference_accessions = manipulate_biosqldb.taxon_id2accessions(server, reference_taxon, biodb)

            print "reference_accessions", reference_accessions
            record_list = []
            for accession in reference_accessions:

                print "reference accession", accession
                biorecord = cache.get(biodb + "_" + accession)

                if not biorecord:
                    print biodb + "_" + accession, "NOT in memory"
                    new_record = db.lookup(accession=accession)
                    biorecord = SeqRecord(Seq(new_record.seq.data, new_record.seq.alphabet),
                                                             id=new_record.id, name=new_record.name,
                                                             description=new_record.description,
                                                             dbxrefs =new_record.dbxrefs,
                                                             features=new_record.features,
                                                             annotations=new_record.annotations)
                    record_id = biorecord.id.split(".")[0]
                    cache.set(biodb + "_" + record_id, biorecord)
                    record_list.append(biorecord)
                else:
                    record_list.append(biorecord)


            if 'submit_circos' in request.POST:

                ref_name = ''
                for i in reference_accessions:
                    ref_name += i
                circos_file = "circos/%s.svg" % ref_name
                import circos
                import shell_command


                querries = manipulate_biosqldb.get_genome_accessions(server, biodb)
                target_taxons = form.cleaned_data['targets']



                target_accessions = [manipulate_biosqldb.taxon_id2accessions(server,int(i),biodb)[0] for i in target_taxons]

                print 'targets!', target_accessions

                target_accessions += reference_accessions
                print target_accessions

                draft_data = []
                for biorecord in record_list:
                    draft_data.append(gbk2circos.circos_fasta_draft_misc_features(biorecord))

                home_dir = os.path.dirname(__file__)
                print "home_dir", home_dir
                temp_location = os.path.join(home_dir, "../assets/circos/")
                myplot = circos.CircosAccession2multiplot(server,
                                          db,
                                          biodb,
                                          record_list,
                                          target_accessions,
                                          locus_highlight=[],
                                          out_directory=temp_location,
                                          draft_fasta=draft_data,
                                          href="/chlamdb/locusx/%s/" % biodb)



                original_map_file = settings.BASE_DIR + "/assets/circos/%s.html" % ref_name
                with open(original_map_file, "r") as f:
                    map_string = ''.join([line for line in f.readlines()])

                circos_html = '<!DOCTYPE html>\n' \
                              ' <html>\n' \
                              ' <body>\n' \
                              ' %s\n' \
                              ' <img src="%s.svg" usemap="#%s">' \
                              ' </body>\n' \
                              ' </html>\n' % (map_string, ref_name, ref_name)


                circos_new_file = '/assets/circos/circos_clic.html'

                with open(settings.BASE_DIR + circos_new_file, "w") as f:
                    f.write(circos_html)

                #target_map_file = settings.BASE_DIR + "/templates/circos/%s.html" % ref_name
                original_map_file_svg = settings.BASE_DIR + "/assets/circos/%s.svg" % ref_name
                #target_map_file_svg = settings.BASE_DIR + "/templates/circos/%s.svg" % ref_name
                map_file = "circos/%s.html" % ref_name
                svg_file = "circos/%s.svg" % ref_name
                #a, b, c = shell_command.shell_command("mv %s %s" % (original_map_file, target_map_file))
                #a, b, c = shell_command.shell_command("cp %s %s" % (original_map_file_svg, target_map_file_svg))
                #print a,b,c
                map_name = ref_name




                envoi_circos = True

            if 'submit_region' in request.POST:
                envoi_region = True




                record = db.lookup(accession=reference_accession)
                start_stop = form.cleaned_data['region'].split(",")
                print "tart_stop", start_stop
                reformat_record = SeqRecord(Seq(record.seq.data, record.seq.alphabet), id=record.id, name=record.name,
                                            description=record.description, dbxrefs=record.dbxrefs,
                                            features=record.features,
                                            annotations=record.annotations)
                print reformat_record
                print "start,", start_stop[0], start_stop[1]
                sub_record = reformat_record[int(start_stop[0]):int(start_stop[1])]
                print "sub"
                print sub_record
                data = Record(sub_record)
                print data.features
                header = ["contig", "type", "start", "stop", "length", "GC", "strand", "gene", "function", "inference",
                          "gi", "locus", "translation"]

                result = []
                for feature in data.features:
                    result.append(
                        [feature.contig, feature.type, feature.start, feature.stop, feature.length, feature.GC,
                         feature.strand, feature.gene, feature.product, feature.inference, feature.gi, feature.locus,
                         feature.translation])

            envoi_region = True
    else:
        form = circos_form_class()
    return render(request, 'chlamdb/circos.html', locals())

@login_required
def alignment(request, input_fasta):
    print align
    handle = open(input_fasta, "rU")
    for record in SeqIO.parse(handle, "fasta"):
        pass
    return render(request, 'chlamdb/alignment.html', locals())






def format_seqfeature_values(server, biodb, seqfeature_id):


    seqfeature_data = manipulate_biosqldb.seqfeature_id2seqfeature_qualifier_values(server, seqfeature_id, biodb)
    if not 'translation' in seqfeature_data.keys():
        # TODO add handeling of other kind of features than CDS
        return None
    try:
        sql_family_size = 'select count(*) as `n_rows` from ' \
           ' (select * from orthology_detail_chlamydia_02_15 ' \
           ' where orthogroup = "%s" group by taxon_id) a' % (seqfeature_data["orthogroup"])
        sql_n_homologues = 'select count(*) as `n_rows` from ' \
           ' (select * from orthology_detail_chlamydia_02_15 ' \
           ' where orthogroup = "%s") a' % (seqfeature_data["orthogroup"])
        n_family = int(server.adaptor.execute_and_fetchall(sql_family_size,)[0][0])
        n_homologues = int(server.adaptor.execute_and_fetchall(sql_n_homologues,)[0][0])
    except KeyError:
        n_family = '-'
        n_homologues = '-'



    # 0 seqfeature_values["description"],
    # 1 seqfeature_values["locus_tag"],
    # 2 seqfeature_values["protein_id"],
    # 3 seqfeature_values["product"],
    # 4 seqfeature_values["orthogroup"],
    # 5 seqfeature_values["gene"],
    # 6 seqfeature_values["translation"],
    # 7 n_family,
    # 8 n_homologues
    template = ["-"] * 9
    try:
        template[0] = seqfeature_data["description"]
    except KeyError:
        pass
    try:
        template[1] = seqfeature_data["locus_tag"]
    except KeyError:
        pass
    try:
        template[2] = seqfeature_data["protein_id"]
    except KeyError:
        pass
    try:
        template[3] = seqfeature_data["product"]
    except KeyError:
        pass
    try:
        template[4] = seqfeature_data["orthogroup"]
    except KeyError:
        pass
    try:
        template[5] = seqfeature_data["gene"]
    except KeyError:
        pass
    try:
        template[6] = seqfeature_data["translation"]
    except KeyError:
        pass

    template[7] = n_family
    template[8] = n_homologues

    return template
    '''
        try:
            return []
        except KeyError:
            try:
                return [seqfeature_values["description"], "-", seqfeature_values["protein_id"], seqfeature_values["product"], seqfeature_values["orthogroup"], seqfeature_values["gene"], seqfeature_values["translation"], n_family, n_homologues]
            except KeyError:
                try:
                    return [seqfeature_values["description"], seqfeature_values["locus_tag"] , "-", seqfeature_values["product"], seqfeature_values["orthogroup"], seqfeature_values["gene"], seqfeature_values["translation"], n_family, n_homologues]
                except KeyError:
                    try:
                        return [seqfeature_values["description"], "-" , "-", seqfeature_values["product"], seqfeature_values["orthogroup"], seqfeature_values["gene"], seqfeature_values["translation"], n_family, n_homologues]
                    except KeyError:
                        try:
                            return [seqfeature_values["description"], seqfeature_values["locus_tag"] , seqfeature_values["protein_id"], seqfeature_values["product"], seqfeature_values["orthogroup"], "-", seqfeature_values["translation"], n_family, n_homologues]
                        except KeyError:
                            return [seqfeature_values["description"], seqfeature_values["locus_tag"] , seqfeature_values["protein_id"], "-", seqfeature_values["orthogroup"], "-", seqfeature_values["translation"], n_family, n_homologues]
    '''




def format_search(count, seqfeature_data):
    # [y, i["description"], i["locus_tag"] + " / " + i["protein_id"], i["product"], i["orthogroup"],  i["gene"], i["translation"]]
    pass


@login_required
def search_taxonomy(request, biodb):
    from collections import Counter
    server = manipulate_biosqldb.load_db()

    if request.method == 'POST':  # S'il s'agit d'une requête POST

        print request.POST

        if len(request.POST['Phylum']) == 0:
            genome_accession = request.POST['Genome']
            superkingdom = request.POST['Superkingdom'].split('_')[-1]
            print "genome_accession", genome_accession
            print 'sdfsdf', superkingdom[0]
            if 'unclassified' in superkingdom:
                superkingdom = '-'

            if len(superkingdom) == 0:
                sql= 'select t1.locus_tag, ' \
                          ' t2.subject_taxon_id,' \
                          ' t4.nr_hit_id,' \
                          ' t1.subject_accession, ' \
                          ' t3.kingdom, ' \
                          ' t3.phylum, ' \
                          ' t3.order, ' \
                          ' t3.family, ' \
                          ' t3.genus, ' \
                          ' t3.species, ' \
                          ' t4.evalue, ' \
                          ' t4.percent_identity, ' \
                          ' t4.query_start, ' \
                          ' t4.query_end, ' \
                          ' t1.subject_title, ' \
                          ' t3.taxon_id ' \
                          ' from blastnr.blastnr_hits_%s_%s as t1 ' \
                          ' inner join blastnr.blastnr_hits_taxonomy_filtered_%s_%s as t2 on t1.nr_hit_id = t2.nr_hit_id ' \
                          ' inner join blastnr.blastnr_taxonomy as t3 on t2.subject_taxon_id = t3.taxon_id' \
                          ' inner join blastnr.blastnr_hsps_%s_%s as t4 on t1.nr_hit_id=t4.nr_hit_id' \
                          ' where t1.hit_number=%s"' % (biodb,
                                                         genome_accession,
                                                         biodb,
                                                         genome_accession,
                                                         biodb,
                                                         genome_accession,
                                                         1)

            else:
                sql= 'select t1.locus_tag, ' \
                          ' t2.subject_taxon_id,' \
                          ' t4.nr_hit_id,' \
                          ' t1.subject_accession, ' \
                          ' t3.kingdom, ' \
                          ' t3.phylum, ' \
                          ' t3.order, ' \
                          ' t3.family, ' \
                          ' t3.genus, ' \
                          ' t3.species, ' \
                          ' t4.evalue, ' \
                          ' t4.percent_identity, ' \
                          ' t4.query_start, ' \
                          ' t4.query_end, ' \
                          ' t1.subject_title, ' \
                          ' t3.taxon_id ' \
                          ' from blastnr.blastnr_hits_%s_%s as t1 ' \
                          ' inner join blastnr.blastnr_hits_taxonomy_filtered_%s_%s as t2 on t1.nr_hit_id = t2.nr_hit_id ' \
                          ' inner join blastnr.blastnr_taxonomy as t3 on t2.subject_taxon_id = t3.taxon_id' \
                          ' inner join blastnr.blastnr_hsps_%s_%s as t4 on t1.nr_hit_id=t4.nr_hit_id' \
                          ' where t1.hit_number=%s and t3.superkingdom="%s"' % (biodb,
                                                                         genome_accession,
                                                                         biodb,
                                                                         genome_accession,
                                                                         biodb,
                                                                         genome_accession,
                                                                         1,
                                                                         superkingdom)

        else:
            genome_accession = request.POST['Genome']
            superkingdom = request.POST['Superkingdom'].split('_')[-1]

            if request.POST['Phylum'] == 'Unclassified':
                phylum = '-'
            else:
                phylum = request.POST['Phylum'].split('_')[-1]

            if len(phylum) ==0:
                sql= 'select t1.locus_tag, ' \
                          ' t2.subject_taxon_id,' \
                          ' t4.nr_hit_id,' \
                          ' t1.subject_accession, ' \
                          ' t3.kingdom, ' \
                          ' t3.phylum, ' \
                          ' t3.order, ' \
                          ' t3.family, ' \
                          ' t3.genus, ' \
                          ' t3.species, ' \
                          ' t4.evalue, ' \
                          ' t4.percent_identity, ' \
                          ' t4.query_start, ' \
                          ' t4.query_end, ' \
                          ' t1.subject_title, ' \
                          ' t3.taxon_id ' \
                          ' from blastnr.blastnr_hits_%s_%s as t1 ' \
                          ' inner join blastnr.blastnr_hits_taxonomy_filtered_%s_%s as t2 on t1.nr_hit_id = t2.nr_hit_id ' \
                          ' inner join blastnr.blastnr_taxonomy as t3 on t2.subject_taxon_id = t3.taxon_id' \
                          ' inner join blastnr.blastnr_hsps_%s_%s as t4 on t1.nr_hit_id=t4.nr_hit_id' \
                          ' where t1.hit_number=%s and t3.superkingdom="%s"' % (biodb,
                                                                         genome_accession,
                                                                         biodb,
                                                                         genome_accession,
                                                                         biodb,
                                                                         genome_accession,
                                                                         1,
                                                                         superkingdom)
            else:
                sql= 'select t1.locus_tag, ' \
                          ' t2.subject_taxon_id,' \
                          ' t4.nr_hit_id,' \
                          ' t1.subject_accession, ' \
                          ' t3.kingdom, ' \
                          ' t3.phylum, ' \
                          ' t3.order, ' \
                          ' t3.family, ' \
                          ' t3.genus, ' \
                          ' t3.species, ' \
                          ' t4.evalue, ' \
                          ' t4.percent_identity, ' \
                          ' t4.query_start, ' \
                          ' t4.query_end, ' \
                          ' t1.subject_title, ' \
                          ' t3.taxon_id,' \
                          ' t3.superkingdom' \
                          ' from blastnr.blastnr_hits_%s_%s as t1 ' \
                          ' inner join blastnr.blastnr_hits_taxonomy_filtered_%s_%s as t2 on t1.nr_hit_id = t2.nr_hit_id ' \
                          ' inner join blastnr.blastnr_taxonomy as t3 on t2.subject_taxon_id = t3.taxon_id' \
                          ' inner join blastnr.blastnr_hsps_%s_%s as t4 on t1.nr_hit_id=t4.nr_hit_id' \
                          ' where t1.hit_number=%s and t3.superkingdom="%s" and t3.phylum="%s";' % (biodb,
                                                                         genome_accession,
                                                                         biodb,
                                                                         genome_accession,
                                                                         biodb,
                                                                         genome_accession,
                                                                         1,
                                                                         superkingdom,
                                                                         phylum)
        print sql
        raw_data = server.adaptor.execute_and_fetchall(sql, )

        n = 1
        data = []
        families = []
        phylum = []
        superkingdom = []
        order = []
        for i, one_hit in enumerate(raw_data):
            if n == 1:
                data.append(one_hit + (n,))
                families.append(one_hit[7])
                phylum.append(one_hit[5])
                superkingdom.append(one_hit[-1])
                n+=1
            else:
                if raw_data[i][0] == raw_data[i-1][0]:
                    data.append(one_hit + (n,))
                else:
                    n+=1
                    data.append(one_hit + (n,))
                    families.append(one_hit[7])
                    phylum.append(one_hit[5])
                    superkingdom.append(one_hit[-1])
        if not 'Phylum' in request.POST and not 'Superkingdom' in request.POST:
            classif_table = dict(Counter(superkingdom))

        elif not 'Phylum' in request.POST and 'Superkingdom' in request.POST:
            classif_table = dict(Counter(phylum))

        else:
            classif_table = dict(Counter(families))

        print classif_table

        '''
        import pandas
        frame = pandas.DataFrame(data, columns= ['n',
                                                 'locus_tag',
                                                 'subject_taxon_id',
                                                 'nr_hit_id',
                                                 'subject_accession',
                                                 'kingdom',
                                                 'phylum',
                                                 'order',
                                                 'family',
                                                 'genus',
                                                 'species',
                                                 'evalue',
                                                 'percent_identity','query_start','query_end','subject_title','taxon_id'])

        print frame
        '''
        envoi = True


    return render(request, 'chlamdb/search_taxonomy.html', locals())


@login_required
def interpro(request, biodb):
    server = manipulate_biosqldb.load_db()

    interproform = make_interpro_from(biodb)

    if request.method == 'POST':  # S'il s'agit d'une requête POST

        form = interproform(request.POST)  # Nous reprenons les données
        #form2 = ContactForm(request.POST)
        if form.is_valid():  # Nous vérifions que les données envoyées sont valides
            invalid_id = False
            # Ici nous pouvons traiter les données du formulaire
            search_type = form.cleaned_data['search_type']
            search_term = form.cleaned_data['search_term']
            taxon_ids = form.cleaned_data['targets']
            #biodb = form.cleaned_data['biodatabase']
            server, db = manipulate_biosqldb.load_db(biodb)

            columns = 'accession,' \
                      'locus_tag,' \
                      'organism, ' \
                      'analysis, ' \
                      'signature_accession, ' \
                      'signature_description, ' \
                      'interpro_accession, ' \
                      'interpro_description,' \
                      'start, ' \
                      'stop, ' \
                      'score, ' \
                      'GO_terms'

            if len(taxon_ids) == 0:
                invalid_id = True
            else:

                taxon_limit = '(taxon_id=%s' % taxon_ids[0]
                if len(taxon_ids) > 1:
                    for i in range(1, len(taxon_ids)-1):
                        taxon_limit+= ' or taxon_id=%s' % taxon_ids[i]
                    taxon_limit+=' or taxon_id=%s)' % taxon_ids[-1]
                else:
                    taxon_limit += ')'

                if search_type == "description":

                    sql = 'select %s from interpro_%s where %s and (interpro_description REGEXP "%s" or signature_description REGEXP "%s")' % (columns, biodb, taxon_limit, search_term, search_term)

                if search_type == "GO":
                    sql = 'select %s from interpro_%s where %s and (GO_terms REGEXP "%s")' % (columns, biodb, taxon_limit, search_term)


                if search_type == "EC":
                    sql = 'select %s from interpro_%s where %s and (pathways REGEXP "%s")' % (columns, biodb, taxon_limit, search_term)

                if search_type == "interpro_accession":
                    sql = 'select %s from interpro_%s where %s and (interpro_accession REGEXP "%s")' % (columns, biodb, taxon_limit, search_term)


                try:
                    raw_data = server.adaptor.execute_and_fetchall(sql, )
                except IndexError:
                    invalid_id = True


            envoi = True

    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = interproform()  # Nous créons un formulaire vide

    return render(request, 'chlamdb/interpro.html', locals())

@login_required
def search(request, biodb):
    server = manipulate_biosqldb.load_db()
    if request.method == 'POST':  # S'il s'agit d'une requête POST

        form = SearchForm(request.POST)  # Nous reprenons les données
        #form2 = ContactForm(request.POST)
        if form.is_valid():  # Nous vérifions que les données envoyées sont valides
            invalid_id = False
            # Ici nous pouvons traiter les données du formulaire
            search_type = form.cleaned_data['search_type']
            search_term = form.cleaned_data['search_term']
            #biodb = form.cleaned_data['biodatabase']
            server, db = manipulate_biosqldb.load_db(biodb)
            print "biodb", biodb

            columns = 'orthogroup, locus_tag, protein_id, start, stop, ' \
                      'strand, gene, orthogroup_size, n_genomes, TM, SP, product, organism, translation'

            if search_type == "gene":

                sql = 'select %s from orthology_detail_%s where gene REGEXP "%s"' % (columns, biodb, search_term)
                raw_data = server.adaptor.execute_and_fetchall(sql,)

            if search_type == "product":
                sql = 'select %s from orthology_detail_%s where product REGEXP "%s"' % (columns, biodb, search_term)
                raw_data = server.adaptor.execute_and_fetchall(sql,)

            if search_type == "locus_tag":
                sql = 'select %s from orthology_detail_%s where locus_tag REGEXP "%s"' % (columns, biodb, search_term)
                raw_data = server.adaptor.execute_and_fetchall(sql,)


            n = 1
            search_result = []
            for one_hit in raw_data:
                if one_hit[2] != '-':
                    interpro_id = one_hit[2]
                else:
                    interpro_id = one_hit[1]
                search_result.append((n,) + one_hit + (interpro_id,))
                n+=1
                print n
            print search_result


            '''
            seqfeature_data = []
            for one_id in seqfeature_ids:
                print one_id
                seqfeature_values = manipulate_biosqldb.seqfeature_id2seqfeature_qualifier_values(server, one_id, biodb)
                seqfeature_data.append(seqfeature_values)
            search_result = []
            y = 1
            for i in seqfeature_data:
                if not 'translation' in i.keys():
                    continue
                print i.keys()
                print i
                #  ['locus_tag', 'orthogroup', 'transl_table', 'product', 'translation', 'gene']
                search_result.append(format_search(y, i))


                y+=1
            '''
            envoi = True

    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = SearchForm()  # Nous créons un formulaire vide

    return render(request, 'chlamdb/search.html', locals())



@login_required
def primer_search(request, biodb):
    server = manipulate_biosqldb.load_db()
    if request.method == 'POST':  # S'il s'agit d'une requête POST

        form = PCRForm(request.POST)  # Nous reprenons les données

        if form.is_valid():  # Nous vérifions que les données envoyées sont valides
            from Bio.Blast.Applications import NcbiblastpCommandline
            #from StringIO import StringIO
            from tempfile import NamedTemporaryFile


            from Bio.Alphabet import IUPAC
            import os
            import shell_command
            import re
            def ExtractAlphanumeric(InputString):
                from string import ascii_letters, digits
                return "".join([ch for ch in InputString if ch in (ascii_letters + digits)])

            input_sequence = form.cleaned_data['blast_input']
            input_sequence = ExtractAlphanumeric(input_sequence)
            print input_sequence

            #biodb = form.cleaned_data['biodatabase']
            input_sequence = input_sequence.rstrip(os.linesep)
            print input_sequence
            my_record = SeqRecord(Seq(input_sequence, IUPAC.protein), id="INPUT", description="INPUT")
            print my_record
            query_file = NamedTemporaryFile()
            SeqIO.write(my_record, query_file, "fasta")
            query_file.flush()

            blastdb = settings.BASE_DIR + '/assets/blast_db/%s.faa' % biodb


            blastp_cline = NcbiblastpCommandline(query=query_file.name, db=blastdb, evalue=0.001, outfmt=0)
            stdout, stderr = blastp_cline()
            print "blast!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            print stdout
            blast_file = NamedTemporaryFile()
            blast_file.write(stdout)
            mview_cmd = 'mview -in blast -ruler on -html data -css on -coloring identity %s' % blast_file.name
            stdout, stderr, code = shell_command.shell_command(mview_cmd)
            blast_result = stdout

            #blast_result = NCBIXML.parse(StringIO(stdout))
            #print blast_result
            #blast_record = next(blast_result)
            #print blast_record



            envoi = True

    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = PCRForm()  # Nous créons un formulaire vide

    return render(request, 'chlamdb/pcr.html', locals())


@login_required
def motif_search(request, biodb):
    import shell_command
    import os
    server = manipulate_biosqldb.load_db()

    motif_form_class = make_motif_form(biodb)

    if request.method == 'POST':

        form = motif_form_class(request.POST)

        if form.is_valid():

            from Bio.Emboss.Applications import FuzznucCommandline
            from tempfile import NamedTemporaryFile

            input_pattern = form.cleaned_data['motif_input']
            n_missmatch = form.cleaned_data['n_missmatch']
            target_taxon_id = form.cleaned_data['search_in']
            print "target", target_taxon_id

            if target_taxon_id != "all":
                accessions = manipulate_biosqldb.taxon_id2accessions(server, target_taxon_id, biodb)


            #input_pattern = ExtractAlphanumeric(input_pattern)
            print input_pattern
            '''
            fuzzpro -sequence CHUV_chr_and_plasmid.faa
            -rformat gff
            -auto
            -stdout
            -pattern "K-x(2)-[LIVF]-x(4)-[LIVF]-D-x(3)-R-x(2)-L-x(5)-[LIV]-Y"
            '''

            PROJECT_ROOT = os.path.abspath(os.path.dirname(__file__))



            if target_taxon_id != "all":
                db_path = os.path.join(PROJECT_ROOT,"../assets/%s/faa/" % biodb +accessions[0] + ".faa")
            else:
                db_path = os.path.join(PROJECT_ROOT,"../assets/%s/faa/all.faa" % biodb)
            print "path", db_path

            cmd = 'fuzzpro -sequence %s -pmismatch %s -rformat seqtable -auto -stdout -pattern "%s"' % (db_path ,n_missmatch ,input_pattern)
            print cmd
            std_out, std_err, code = shell_command.shell_command(cmd)
            print std_out
            print std_err
            #fuzznuc_cline = FuzznucCommandline(sequence=genome_db, mismatch=n_missmatch, pattern=input_pattern, stdout=True)#, rformat="srspair")
            #stdout,stderr = fuzznuc_cline()


            envoi = True

    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = motif_form_class()  # empty form

    return render(request, 'chlamdb/motifs.html', locals())



@login_required
def blast(request, biodb):
    server = manipulate_biosqldb.load_db()

    blast_form_class = make_blast_form(biodb)

    if request.method == 'POST':  # S'il s'agit d'une requête POST

        form = blast_form_class(request.POST)  # Nous reprenons les données

        if form.is_valid():  # Nous vérifions que les données envoyées sont valides
            from Bio.Blast.Applications import NcbiblastpCommandline
            from Bio.Blast.Applications import NcbiblastnCommandline
            from Bio.Blast.Applications import NcbitblastnCommandline
            from Bio.Blast.Applications import NcbiblastxCommandline
            from tempfile import NamedTemporaryFile

            from Bio.Alphabet import IUPAC
            import os
            import shell_command
            import re


            input_sequence = form.cleaned_data['blast_input']

            target_accession = form.cleaned_data['target']

            blast_type = form.cleaned_data['blast']


            input_sequence = extract_alphanumeric(input_sequence)
            print "one", input_sequence

            input_sequence = input_sequence.rstrip(os.linesep)
            print "two", input_sequence

            my_record = SeqRecord(Seq(input_sequence, IUPAC.protein), id="INPUT", description="INPUT")

            query_file = NamedTemporaryFile()
            SeqIO.write(my_record, query_file, "fasta")
            query_file.flush()


            if blast_type=='blastn_ffn':
                blastdb = settings.BASE_DIR + "/assets/%s/ffn/%s.ffn" % (biodb, target_accession)
                blast_cline = NcbiblastnCommandline(query=query_file.name, db=blastdb, evalue=10, outfmt=0)
            if blast_type=='blastn_fna':
                blastdb = settings.BASE_DIR + "/assets/%s/fna/%s.fna" % (biodb, target_accession)
                blast_cline = NcbiblastnCommandline(query=query_file.name, db=blastdb, evalue=10, outfmt=0)
            if blast_type=='blastp':
                blastdb = settings.BASE_DIR + "/assets/%s/faa/%s.faa" % (biodb, target_accession)
                blast_cline = NcbiblastpCommandline(query=query_file.name, db=blastdb, evalue=10, outfmt=0)
            if blast_type=='tblastn':
                blastdb = settings.BASE_DIR + "/assets/%s/fna/%s.fna" % (biodb, target_accession)
                blast_cline = NcbitblastnCommandline(query=query_file.name, db=blastdb, evalue=10, outfmt=0)
                blast_cline2 = NcbitblastnCommandline(query=query_file.name, db=blastdb, evalue=10, outfmt=5)
            if blast_type=='blastx':
                blastdb = settings.BASE_DIR + "/assets/%s/faa/%s.faa" % (biodb, target_accession)
                blast_cline = NcbiblastxCommandline(query=query_file.name, db=blastdb, evalue=10, outfmt=0)


            print blast_cline
            blast_stdout, blast_stderr = blast_cline()

            if blast_type=='tblastn':
                from Bio.SeqUtils import six_frame_translations
                from StringIO import StringIO
                blast_stdout2, blast_stderr2 = blast_cline2()
                from Bio.Blast import NCBIXML
                blast_records = NCBIXML.parse(StringIO(blast_stdout2))
                all_data = []
                for record in blast_records:
                    for alignment in record.alignments:
                        accession = alignment.title.split(' ')[1]
                        sql = 'select description from bioentry where accession="%s" ' % accession
                        description = server.adaptor.execute_and_fetchall(sql,)[0][0]
                        for hsp in alignment.hsps:
                            start = hsp.sbjct_start
                            end = hsp.sbjct_end
                            length = end-start
                            #print 'seq for acc', accession, start, end,
                            leng = end-start
                            seq = manipulate_biosqldb.location2sequence(server, accession, biodb, start, leng)
                            print seq
                            from Bio.Seq import reverse_complement, translate
                            anti = reverse_complement(seq)
                            comp = anti[::-1]
                            length = len(seq)
                            frames = {}
                            for i in range(0, 3):
                                fragment_length = 3 * ((length-i) // 3)
                                tem1 = translate(seq[i:i+fragment_length], 1)
                                frames[i+1] = '<span style="color: #181407;">%s</span><span style="color: #bb60d5;">%s</span><span style="color: #181407;">%s</span>' % (tem1[0:100], tem1[100:len(tem1)-99], tem1[len(tem1)-99:])
                                tmp2 = translate(anti[i:i+fragment_length], 1)[::-1]
                                frames[-(i+1)] = tmp2
                            all_data.append([accession, start, end, length, frames[1], frames[2], frames[3], frames[-1], frames[-2], frames[-3], description, seq])


            no_match = re.compile('.* No hits found .*', re.DOTALL)

            if no_match.match(blast_stdout):
                print "no blast hit"
                blast_no_hits = blast_stdout
            elif len(blast_stderr) != 0:
                print "blast error"
                blast_err = blast_stderr
            else:
                print "running mview"
                blast_file = NamedTemporaryFile()

                f = open('/tmp/blast.temp', 'w')
                f.write(blast_stdout)
                f.close()



                #blast_file.write(blast_stdout)
                out, err, code = shell_command.shell_command('cat /temp/blast.temp | wc -l')
                print 'n lines', out
                if blast_type=='blastp' or blast_type=='blastn_ffn':
                    mview_cmd = 'mview -in blast -srs on -ruler on -html data -css on -coloring identity /tmp/blast.temp' #% blast_file.name
                else:
                    mview_cmd = 'mview -in blast -ruler on -html data -css on -coloring identity /tmp/blast.temp'
                stdout, stderr, code = shell_command.shell_command(mview_cmd)

                if len(stdout) == 0:
                    blast_no_hits = blast_stdout
                    blast_result = None
                else:
                    blast_result = stdout
            #blast_result = NCBIXML.parse(StringIO(stdout))
            #print blast_result
            #blast_record = next(blast_result)
            #print blast_record

            envoi = True

    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = blast_form_class()  # Nous créons un formulaire vide

    return render(request, 'chlamdb/blast.html', locals())


def get_record_from_memory(biodb, cache_obj, record_key, accession):

        biorecord = cache_obj.get(record_key)
        if not biorecord:
            print record_key, "NOT in memory"
            new_record = biodb.lookup(accession=accession)
            biorecord = SeqRecord(Seq(new_record.seq.data, new_record.seq.alphabet),
                                                             id=new_record.id, name=new_record.name,
                                                             description=new_record.description,
                                                             dbxrefs =new_record.dbxrefs,
                                                             features=new_record.features,
                                                             annotations=new_record.annotations)
            cache_obj.set(record_key, biorecord)
        return biorecord





@login_required
def mummer(request, biodb):

    server = manipulate_biosqldb.load_db()
    mummer_form_class = make_mummer_form(biodb)

    cache = get_cache('default')

    if request.method == 'POST':  # S'il s'agit d'une requête POST
        #make_circos2genomes_form
        plot = True
        form = mummer_form_class(request.POST)
        if form.is_valid():  # Nous vérifions que les données envoyées sont valides
            server, db = manipulate_biosqldb.load_db(biodb)
            reference_taxon = form.cleaned_data['reference_genome']
            query_taxon = form.cleaned_data['query_genome']

            reference_accessions = manipulate_biosqldb.taxon_id2accessions(server, reference_taxon, biodb)
            query_accessions = manipulate_biosqldb.taxon_id2accessions(server, query_taxon, biodb)


            ref_accession = manipulate_biosqldb.taxon_id2chromosome_accession(server, biodb, reference_taxon)
            query_accession = manipulate_biosqldb.taxon_id2chromosome_accession(server, biodb, query_taxon)
            #result = result[0]
            print "mummer acc", ref_accession, query_accession




            print settings.STATIC_ROOT, type(settings.STATIC_ROOT)

            reference_path = settings.BASE_DIR + '/assets/%s/fna/%s.fna' % (biodb, ref_accession)
            query_path = settings.BASE_DIR + '/assets/%s/fna/%s.fna' % (biodb, query_accession)

            rand = id_generator(5)

            out_delta = settings.BASE_DIR + '/assets/temp/promer_%s' % rand
            out_plot = settings.BASE_DIR + '/assets/temp/promer_%s' % rand




            cmd1 = 'promer -l 2 -p %s %s %s' % (out_delta, reference_path, query_path)
            cmd2 = 'mummerplot -layout -small -png -p %s %s.delta' % (out_plot, out_delta)

            print cmd1

            from shell_command import shell_command

            out, err, log = shell_command(cmd1)
            out, err, log = shell_command(cmd2)

            plot_path = 'temp/promer_%s.png' % rand

            if not os.path.exists(settings.BASE_DIR + '/assets/' + plot_path):
                plot = False

            envoi = True

    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = mummer_form_class()  # Nous créons un formulaire vide

    return render(request, 'chlamdb/mummer.html', locals())





@login_required
def circos2genomes(request, biodb):
    import circos
    import shell_command
    server = manipulate_biosqldb.load_db()
    circos2genomes_form_class = make_circos2genomes_form(biodb)

    cache = get_cache('default')

    if request.method == 'POST':  # S'il s'agit d'une requête POST
        #make_circos2genomes_form

        form = circos2genomes_form_class(request.POST)
        if form.is_valid():  # Nous vérifions que les données envoyées sont valides
            valid_id = True
            server, db = manipulate_biosqldb.load_db(biodb)
            reference_taxon = form.cleaned_data['reference_genome']
            query_taxon = form.cleaned_data['query_genome']
            import re
            protein_locus_list = re.sub(" ", "", form.cleaned_data['locus_list'])
            protein_locus_list = list(protein_locus_list.split(","))

            #reference_taxon = manipulate_biosqldb.description2taxon_id(server, reference_genome, biodb)
            #query_taxon = manipulate_biosqldb.description2taxon_id(server, query_genome, biodb)

            reference_accessions = manipulate_biosqldb.taxon_id2accessions(server, reference_taxon, biodb)
            query_accessions = manipulate_biosqldb.taxon_id2accessions(server, query_taxon, biodb)

            reference_records = []
            for accession in reference_accessions:
                reference_records.append(get_record_from_memory(db, cache, biodb + "_" + accession, accession))

            query_records = []
            for accession in query_accessions:

                query_records.append(get_record_from_memory(db, cache, biodb + "_" + accession, accession))

            print "genomes", reference_records, query_records

            orthogroup_list = []
            print "protein_locus_list", protein_locus_list[0], len(protein_locus_list[0])
            if len(protein_locus_list[0]) > 0:
                print 'okkkkkkkk'
                for protein in protein_locus_list:
                    print "protein", protein
                    sql = 'select orthogroup from orthology_detail_%s where protein_id="%s" or locus_tag="%s"' % (biodb, protein, protein)
                    print sql
                    try:
                        protein_group = server.adaptor.execute_and_fetchall(sql,)[0][0]
                        orthogroup_list.append(protein_group)
                    except IndexError:
                        valid_id = False

                    #protein_group = manipulate_biosqldb.locus_tag2orthogroup_id(server, protein, biodb)
            if valid_id:
                print "LOCUS:", orthogroup_list

                #accession2description = manipulate_biosqldb.accession2description_dict(server, biodb)

                #print "reference_genome", reference_genome
                taxon_id2description = manipulate_biosqldb.taxon_id2genome_description(server, biodb)
                reference_name = taxon_id2description[reference_taxon]
                query_name = taxon_id2description[query_taxon]

                accession2taxon_id = manipulate_biosqldb.accession2taxon_id(server, biodb)
                taxon_id_reference = accession2taxon_id[reference_accessions[0]]
                taxon_id_query = accession2taxon_id[query_accessions[0]]

                reference_n_orthogroups = manipulate_biosqldb.get_genome_number_of_orthogroups(server, biodb, taxon_id_reference)
                reference_n_proteins = manipulate_biosqldb.get_genome_number_of_proteins(server, biodb, taxon_id_reference)

                query_n_orthogroups = manipulate_biosqldb.get_genome_number_of_orthogroups(server, biodb, taxon_id_query)
                query_n_proteins = manipulate_biosqldb.get_genome_number_of_proteins(server, biodb, taxon_id_query)

                n_shared_orthogroups = manipulate_biosqldb.get_number_of_shared_orthogroups(server, biodb, taxon_id_reference, taxon_id_query)

                print "n_orthogroups", reference_n_orthogroups, query_n_orthogroups
                print "n_proteins", reference_n_proteins, query_n_proteins

                print "n_shared", n_shared_orthogroups

                print reference_name, query_name


                import circos

                path = settings.BASE_DIR + "/assets/circos"



                biplot = circos.CircosAccession2biplot(server, db, biodb, reference_records, query_records,
                                                       orthogroup_list, path)

                reference_file = "circos/%s" % biplot.reference_circos
                query_file = "circos/%s" % biplot.query_circos

            envoi = True

    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = circos2genomes_form_class()  # Nous créons un formulaire vide

    return render(request, 'chlamdb/circos2genomes.html', locals())

@login_required
def circos2genomes_main(request):
    return render(request, 'chlamdb/circos2genomes_main.html', locals())


def update_db(server):
    biodb_list = manipulate_biosqldb.get_biodatabase_list(server)
    for biodb in biodb_list:
        update_genomes_db(server, biodb)

    print

def update_genomes_db(server, biodatabase_name):
    from models import Genome
    from models import Database
    from models import GenDB

    print "updating genome Database"

    database = Database.objects.get_or_create(db_name=biodatabase_name)[0]
    genome_list = manipulate_biosqldb.get_genome_description_list(server, biodatabase_name)
    for genome_description in genome_list:
        print "updating", genome_description
        genome = Genome.objects.get_or_create(genome_name=genome_description, database=database)[0]
        print "genome", genome
        GenDB.objects.get_or_create(database=database, ref_genome=genome, query_genome=genome, genome_name=genome.genome_name, database_name=database.db_name)
        print "OKKKKKKKKKKKK"




#def get_genomes(server, biodb_name):













@login_required
def crossplot(request):

    cache = get_cache('default')
    print "cache", cache
    cache.clear()
    bioentry_in_memory = cache.get('biodb')
    print "bioentry_in_memory", bioentry_in_memory
    if not bioentry_in_memory:
        print "creating cache entry"
        cache.set("biodb", {})
    bioentry_in_memory = cache.get("biodb")
    server = manipulate_biosqldb.load_db()

    update_db(server)

    crossplot_form_class = make_crossplot_form("Chlamydia_11_14")
    if request.method == 'POST':  # S'il s'agit d'une requête POST
        #make_circos2genomes_form

        form = DBForm(request.POST) #crossplot_form_class(request.POST)


        if form.is_valid():  # Nous vérifions que les données envoyées sont valides
            server, db = manipulate_biosqldb.load_db("Chlamydia_11_14")

            update_genomes_db(server, "Chlamydia_11_14")

            print "form", form
            print dir(form)
            print form.fields
            print "data1", form.cleaned_data


            reference_genome = str(form.cleaned_data['ref_genome'])
            query_genome = str(form.cleaned_data['query_genome'])
            protein_locus = form.cleaned_data['accession']
            region_size = form.cleaned_data['region_size']

            print "asdfasfffffffffffff", reference_genome


            #description2accession = manipulate_biosqldb.description2accession(server, "saureus1")

            #reference_accession = description2accession[reference_genome]
            #query_accession = description2accession[query_genome]

            orthogroup = manipulate_biosqldb.locus_tag2orthogroup_id(server, protein_locus, "Chlamydia_11_14")
            print "orthogroup", orthogroup
            ortho_detail = list(manipulate_biosqldb.orthogroup_id2locus_tag_list(server, orthogroup, "Chlamydia_11_14"))
            locus_tag_list = []
            for i in range(0, len(ortho_detail)):
                print "ortho detail", ortho_detail[i]
                if reference_genome in ortho_detail[i] or query_genome in ortho_detail[i]:
                    locus_tag_list.append(ortho_detail[i][2])
            print "locus_lis", locus_tag_list


            if "Chlamydia_11_14" not in bioentry_in_memory.keys():
                bioentry_in_memory["Chlamydia_11_14"] = {}

            home_dir = os.path.dirname(__file__)
            print "home_dir", home_dir
            temp_location = os.path.join(home_dir, "../assets")
            print "temp loc", temp_location
            temp_file = NamedTemporaryFile(delete=False, dir=temp_location, suffix=".svg")
            print "temp file", temp_file.name
            name = os.path.basename(temp_file.name)
            bioentry_dict = mysqldb_plot_genomic_feature.proteins_id2cossplot(server, db, "Chlamydia_11_14", locus_tag_list,
                                                                                  temp_file.name, int(region_size),
                                                                                  bioentry_in_memory["Chlamydia_11_14"])




            envoi = True

    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = DBForm() #crossplot_form_class()  # Nous créons un formulaire vide
        form2 = DBForm()
    return render(request, 'chlamdb/crossplot.html', locals())





def pfam_tree(request, biodb, orthogroup):
    import ete_motifs
    import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)

    sql_locus2protein_id = 'select locus_tag, protein_id from orthology_detail_%s where orthogroup="%s"' % (biodb, orthogroup)

    locus2protein_id= manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_locus2protein_id,))

    locus2pfam_data = ete_motifs.get_pfam_data(orthogroup, biodb)

    motif_count = {}
    for data in locus2pfam_data.values():
        for motif in data:
            print data
            try:
                if motif[4] not in motif_count:
                    motif_count[motif[4]] = [1, motif[5]]
                else:
                    motif_count[motif[4]][0]+=1
            except:
                print "motif", motif

    print "motif_count", motif_count

    sql_tree = 'select phylogeny from biosqldb_phylogenies.%s where orthogroup="%s"' % (biodb, orthogroup)
    try:
        tree = server.adaptor.execute_and_fetchall(sql_tree,)[0][0]
    except IndexError:
        no_tree = True
        return render(request, 'chlamdb/pfam_tree.html', locals())
    import manipulate_biosqldb

    sql = 'select taxon_id, family from genomes_classification;'

    taxon_id2family = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    print 'tree', tree
    t, ts, leaf_number = ete_motifs.draw_pfam_tree(tree, locus2pfam_data, locus2protein_id, taxon_id2family)
    path = settings.BASE_DIR + '/assets/temp/pfam_tree.svg'
    asset_path = '/assets/temp/pfam_tree.svg'
    #print "path", path
    t.render(path, h=leaf_number*12, dpi=800, tree_style=ts)

    return render(request, 'chlamdb/pfam_tree.html', locals())

def TM_tree(request, biodb, orthogroup):
    print 'bonjour', request.method
    import manipulate_biosqldb
    import ete_motifs

    server, db = manipulate_biosqldb.load_db(biodb)

    locus2TM_data = ete_motifs.get_TM_data(orthogroup, biodb)

    sql_tree = 'select phylogeny from biosqldb_phylogenies.%s where orthogroup="%s"' % (biodb, orthogroup)
    tree = server.adaptor.execute_and_fetchall(sql_tree,)[0]
    print 'tree', tree
    t, ts, leaf_number = ete_motifs.draw_TM_tree(tree, locus2TM_data)
    path = settings.BASE_DIR + '/assets/temp/tm_tree.svg'
    #print "path", path
    t.render(path, h=leaf_number*12, dpi=800, tree_style=ts)

    return render(request, 'chlamdb/pfam_tree.html', locals())





def orthogroup_conservation_tree(request, biodb, orthogroup):

    import manipulate_biosqldb
    import shell_command

    server, db = manipulate_biosqldb.load_db(biodb)

    asset_path = '/assets/temp/phylo.svg'
    path = settings.BASE_DIR + asset_path
    a,b,c = shell_command.shell_command("rm %s" % path)
    print a, b, c
    import ete_heatmap_conservation
    sql_grp = 'select taxon_id,count(*) from  orthology_detail_%s where orthogroup="%s" group by organism;' % (biodb, orthogroup)
    print sql_grp
    taxid2n = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_grp,))
    tree_sql = 'select tree from reference_phylogeny as t1 inner join biodatabase as t2 on t1.biodatabase_id=t2.biodatabase_id where t2.name="%s"' % biodb
    tree = server.adaptor.execute_and_fetchall(tree_sql,)[0][0]

    print settings.BASE_DIR
    print path
    print taxid2n
    t1, leaf_number = ete_heatmap_conservation.plot_heat_tree(biodb, taxid2n, tree)

    shell_command.shell_command('rm %s' % path)
    print path
    t1.render(path, dpi=800, h=leaf_number*12)


    return render(request, 'chlamdb/orthogroup_conservation.html', locals())



@login_required
def priam_kegg(request, biodb):

    priam_form_class = make_priam_form(biodb)

    print 'request', request.method

    if request.method == 'POST':  # S'il s'agit d'une requête POST
        print 'request', request.method
        form = priam_form_class(request.POST)  # Nous reprenons les données
        #form2 = ContactForm(request.POST)
        print 'aaa'
        if form.is_valid():
            genome = form.cleaned_data['genome']



            envoi = True

    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = priam_form_class()  # Nous créons un formulaire vide

    return render(request, 'chlamdb/priam_kegg.html', locals())
