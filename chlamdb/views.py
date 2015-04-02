# -*- coding: utf-8 -*-

# todo circos gc file curently written in home directory, move it to other place
# todo save temp files in temp folder


#from django.shortcuts import render
#from datetime import datetime
from django.shortcuts import render
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
from forms import BlastForm
from forms import make_crossplot_form
from forms import ConnexionForm
from forms import DBForm
from forms import make_motif_form
from forms import PCRForm
from forms import make_extract_form
from Bio.SeqRecord import SeqRecord
from django.contrib.auth import logout

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

def extract_alphanumeric(input_string):
    from string import ascii_letters, digits
    return "".join([ch for ch in input_string if ch in (ascii_letters + digits + '_-')])


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
    #table = create_table_from_dict(tata)

    return render(request, 'chlamdb/home.html', locals())


def substription():
    create_user('tpillone', 'trestan.pillonel@gmail.com', 'estrella3', "Trestan", "Pillonel")

#cache = pylibmc.Client(['127.0.0.1:8000'])



@login_required
def extract(request, biodb):

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

            sql_include = ''
            for i in range(0, len(include)-1):
                sql_include += ' `%s` > 0 and ' % include[i]
            sql_include+='`%s` > 0 and ' % include[-1]

            sql_exclude = ''
            for i in range(0, len(exclude)-1):
                sql_exclude += ' `%s` = 0 and ' % exclude[i]
            sql_exclude+='`%s` = 0' % exclude[-1]

            server, db = manipulate_biosqldb.load_db(biodb)

            sql ='select orthogroup from orthology_%s where %s %s' % (biodb, sql_include, sql_exclude)


            match_groups = [ i[0] for i in server.adaptor.execute_and_fetchall(sql,)]

            sql_include = 'taxon_id ='
            for i in range(0, len(include)-1):
                sql_include+='%s or taxon_id =' % include[i]
            sql_include+=include[-1]
            n = 1
            search_result = []
            for group in match_groups:
                sql_2 = 'select seqfeature_id from orthology_detail_%s where orthogroup = "%s" and (%s)' % (biodb, group, sql_include)
                all_features = [ i[0] for i in server.adaptor.execute_and_fetchall(sql_2,)]
                for feature_id in all_features:

                    data = format_seqfeature_values(server, biodb, feature_id)
                    search_result.append([n] + data)
                    n+=1
                    print n
            print search_result
            envoi_extract = True

    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = extract_form_class()  # Nous créons un formulaire vide

    return render(request, 'chlamdb/extract_genes.html', locals())





@login_required
def homology(request, biodb):

    cache = get_cache('default')
    print "cache", cache
    #cache.clear()

    #bioentry_in_memory = cache.get("biodb")
    print "loading db..."
    server = manipulate_biosqldb.load_db()
    print "db loaded..."
    contact_form_class = make_contact_form(server, biodb)
    if request.method == 'POST':  # S'il s'agit d'une requête POST

        """
        bioentry_in_memory = cache.get('biodb')
        print "bioentry_in_memory", bioentry_in_memory
        if not bioentry_in_memory:
            print "creating cache entry"
            cache.set("biodb", {})
        """

        form = contact_form_class(request.POST)  # Nous reprenons les données
        #form2 = ContactForm(request.POST)
        if form.is_valid():  # Nous vérifions que les données envoyées sont valides
            valid_id = True
            # Ici nous pouvons traiter les données du formulaire
            accession = extract_alphanumeric(form.cleaned_data['accession'])
            #biodb = form.cleaned_data['biodatabase']
            plot_region = form.cleaned_data['plot_region']
            server, db = manipulate_biosqldb.load_db(biodb)
            print "plot_region", plot_region
            region_size = form.cleaned_data['region_size']
            print "region size", region_size
            print "biodb", biodb
            #if biodb not in bioentry_in_memory.keys():
            #    bioentry_in_memory[biodb] = {}



            locus_tag2seqfeature_id = cache.get("locus_tag2seqfeature_id_%s" % biodb)

            if locus_tag2seqfeature_id is None:
                locus_tag2seqfeature_id = manipulate_biosqldb.locus_tag2seqfeature_id_dict(server, biodb)
                cache.set("locus_tag2seqfeature_id_%s" % biodb, locus_tag2seqfeature_id)
            else:
                print "locus_tag2seqfeature_id in memory"

            protein_id2seqfeature_id = cache.get("protein_id2seqfeature_id_%s" % biodb)

            if protein_id2seqfeature_id is None:
                protein_id2seqfeature_id = manipulate_biosqldb.protein_id2seqfeature_id_dict(server, biodb)
                cache.set("protein_id2seqfeature_id_%s" % biodb, protein_id2seqfeature_id)
            else:
                print "protein_id2seqfeature_id in memory"

            try:
                seqid = locus_tag2seqfeature_id[accession]
                #seqfeature_values = manipulate_biosqldb.locus_tag2seqfeature_qualifier_values(server, accession, biodb)
                #print "seqfeature_values", seqfeature_values
            except KeyError:
                try:
                    print "not a locus tag!"
                    print "trying to convert protein id to locus tag..."
                    #print protein_id2seqfeature_id.keys()[1:10]
                    seqid = protein_id2seqfeature_id[accession]

                    accession = manipulate_biosqldb.protein_id2locus_tag(server, accession, biodb)
                    #print "locus_tag", accession
                    #seqfeature_values = manipulate_biosqldb.locus_tag2seqfeature_qualifier_values(server, accession, biodb)
                    #print "seqfeature_values", seqfeature_values
                except KeyError:
                    valid_id = False

            if valid_id:
                print "getting seqfeature values for seqfeature id %s" % seqid
                seqfeature_values = manipulate_biosqldb.seqfeature_id2seqfeature_qualifier_values(server, seqid, biodb)
                print seqfeature_values
                print "getting bioentry data"
                bioentry = manipulate_biosqldb.locus_tag2bioentry_accession_and_description(server, accession, biodb)
                print bioentry
                orthogroup = seqfeature_values["orthogroup"]
                fasta = "%s_fasta/%s.txt" % (biodb, orthogroup)

                alignment = "%s_fasta/%s.html" % (biodb, orthogroup)
                alignment_fasta = "%s_fasta/%s.fa" % (biodb, orthogroup)
                tree_unrooted = "%s_fasta/%s_tree.svg" % (biodb, orthogroup)
                tree_rooted = "%s_fasta/%s_tree_reroot.svg" % (biodb, orthogroup)
                tree_file = "%s_fasta/%s.phy_phyml_tree.txt" % (biodb, orthogroup)

                bioentry_accession = bioentry[0]
                bioentry_name = bioentry[1]

                print "getting orthologs data"
                ortho_detail = list(manipulate_biosqldb.orthogroup_id2locus_tag_list(server, orthogroup, biodb))


                print "ortho_detail", ortho_detail

                new_list = []
                locus_tag_target_genomes = []
                for i in range(0, len(ortho_detail)):

                    if ortho_detail[i][3] in form.cleaned_data['genomes']:
                        locus_tag_target_genomes.append(ortho_detail[i][2])
                    print ortho_detail[i]
                    new_list.append((i + 1,) + ortho_detail[i])

                ortho_detail = new_list

                if len(ortho_detail) >1:
                    orthologs = True
                else:
                    orthologs = False
                import orthogroup_identity_db
                if len(ortho_detail) > 1:
                    orthogroup2identity_dico = orthogroup_identity_db.orthogroup2identity_dico(biodb, orthogroup)
                    print "orthologs", orthologs, len(ortho_detail)
                    for count, value in enumerate(ortho_detail):
                        locus_2 = value[3]
                        ortho_detail[count] = value + (orthogroup2identity_dico[accession][locus_2],)
                        print ortho_detail
                        print "count", count
                        print ortho_detail[count]
                        location = manipulate_biosqldb.seqfeature_id2feature_location(server, int(ortho_detail[count][2]))
                        print "location",location
                        ortho_detail[count] = ortho_detail[count] + location
                        print ortho_detail
                else:
                    ortho_detail[0] = ortho_detail[0] + (100,)



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
                    mysqldb_plot_genomic_feature.proteins_id2cossplot(server, db, biodb, locus_tag_target_genomes,
                                                                                      temp_file.name, int(region_size),
                                                                                      cache)
                    """
                    print "updating dico"
                    bioentry_in_memory[biodb].update(bioentry_dict)

                    print "bioentry_in_memory", bioentry_in_memory
                    print "updating cache"
                    import sys
                    print "cache size",sys.getsizeof(bioentry_in_memory)


                    #import pickle
                    #dat = pickle.dumps(bioentry_in_memory)
                    #print "pickle", len(dat)
                    cache.set("biodb", bioentry_in_memory)
                    #print "updated cache:", cache.get("biodb")
                    #except:
                    #  invalid_id = True
                    # Nous pourrions ici envoyer l'e-mail grâce aux données que nous venons de récupérer

                    """

                    #try:

            envoi = True

    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = contact_form_class()  # Nous créons un formulaire vide

    return render(request, 'chlamdb/homology.html', locals())






@login_required
def plot_region(request, biodb):
    plot_form_class = make_plot_form(biodb)
    cache = get_cache('default')
    print "cache", cache
    server = manipulate_biosqldb.load_db()

    if request.method == 'POST':  # S'il s'agit d'une requête POST


        form = plot_form_class(request.POST)  # Nous reprenons les données

        if form.is_valid():  # Nous vérifions que les données envoyées sont valides
            invalid_id = False
            # Ici nous pouvons traiter les données du formulaire
            location_start = form.cleaned_data['location_start']
            location_stop = form.cleaned_data['location_stop']

            genome = form.cleaned_data['genome']

            server, db = manipulate_biosqldb.load_db(biodb)

            home_dir = os.path.dirname(__file__)
            print "home_dir", home_dir
            temp_location = os.path.join(home_dir, "../assets")
            print "temp loc", temp_location
            temp_file = NamedTemporaryFile(delete=False, dir=temp_location, suffix=".svg")
            print "temp file", temp_file.name
            name = os.path.basename(temp_file.name)

            bioentry = manipulate_biosqldb.description2accession(server, genome, biodb)
            print "bioentry", bioentry
            mysqldb_plot_genomic_feature.location2simpleplot(db, biodb, bioentry, int(location_start), int(location_stop), temp_file.name, cache)

            envoi = True

    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = plot_form_class()  # Nous créons un formulaire vide

    return render(request, 'chlamdb/plot_region.html', locals())

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
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
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

                querries = manipulate_biosqldb.get_genome_accessions(server, biodb)
                target_taxons = form.cleaned_data['targets']



                target_accessions = [manipulate_biosqldb.taxon_id2accessions(server,int(i),biodb)[0] for i in target_taxons]


                target_accessions += reference_accessions
                print target_accessions

                draft_data = []
                for biorecord in record_list:
                    temp = gbk2circos.circos_fasta_draft_misc_features(biorecord)
                    if len(temp) > 0:
                        draft_data.append(temp)
                if len(draft_data) == 0:
                    draft_data = False
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
                                          draft_fasta=draft_data)
                envoi_circos = True
            if 'submit_region' in request.POST:
                envoi_region = True
                from Bio.SeqRecord import SeqRecord
                from Bio.Seq import Seq


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

    sql_family_size = 'select count(*) as `n_rows` from ' \
           ' (select * from orthology_detail_chlamydia_02_15 ' \
           ' where orthogroup = "%s" group by taxon_id) a' % (seqfeature_data["orthogroup"])
    sql_n_homologues = 'select count(*) as `n_rows` from ' \
           ' (select * from orthology_detail_chlamydia_02_15 ' \
           ' where orthogroup = "%s") a' % (seqfeature_data["orthogroup"])

    n_family = int(server.adaptor.execute_and_fetchall(sql_family_size,)[0][0])
    n_homologues = int(server.adaptor.execute_and_fetchall(sql_n_homologues,)[0][0])

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

            if search_type == "gene":
                seqfeature_ids = manipulate_biosqldb.gene_regex2seqfeature_ids(server, biodb, search_term)
            if search_type == "product":
                seqfeature_ids = manipulate_biosqldb.product_regex2seqfeature_ids(server, biodb, search_term)

            print seqfeature_ids
            n = 1
            search_result = []
            for feature_id in seqfeature_ids:

                data = format_seqfeature_values(server, biodb, feature_id)
                if data:
                    search_result.append([n] + data)
                    n+=1
                    print n



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
            from Bio.Seq import Seq
            from Bio.SeqRecord import SeqRecord
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
            blastp_cline = NcbiblastpCommandline(query=query_file.name, db="~/Dropbox/dev/django/test_1/assets/blast_db/%s.faa" % biodb, evalue=0.001, outfmt=0)
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

            accessions = manipulate_biosqldb.taxon_id2accessions(server,target_taxon_id,biodb)
            print "accessions", accessions

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

            db_path = os.path.join(PROJECT_ROOT,"../assets/blast_db/"+accessions[0] + ".faa")
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
    if request.method == 'POST':  # S'il s'agit d'une requête POST

        form = BlastForm(request.POST)  # Nous reprenons les données

        if form.is_valid():  # Nous vérifions que les données envoyées sont valides
            from Bio.Blast.Applications import NcbiblastpCommandline
            #from StringIO import StringIO
            from tempfile import NamedTemporaryFile
            from Bio.Seq import Seq
            from Bio.SeqRecord import SeqRecord
            from Bio.Alphabet import IUPAC
            import os
            import shell_command
            import re


            input_sequence = form.cleaned_data['blast_input']
            input_sequence = extract_alphanumeric(input_sequence)
            print input_sequence

            #biodb = form.cleaned_data['biodatabase']
            input_sequence = input_sequence.rstrip(os.linesep)
            print input_sequence
            my_record = SeqRecord(Seq(input_sequence, IUPAC.protein), id="INPUT", description="INPUT")
            print my_record
            query_file = NamedTemporaryFile()
            SeqIO.write(my_record, query_file, "fasta")
            query_file.flush()
            blastp_cline = NcbiblastpCommandline(query=query_file.name, db="~/Dropbox/dev/django/test_1/assets/blast_db/%s.faa" % biodb, evalue=0.001, outfmt=0)
            stdout, stderr = blastp_cline()
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
        form = BlastForm()  # Nous créons un formulaire vide

    return render(request, 'chlamdb/blast.html', locals())


def get_record_from_memory(biodb, cache_obj, record_key, accession):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
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
            if len(protein_locus_list) > 0:
                for protein in protein_locus_list:
                    print "protein", protein
                    protein_group = manipulate_biosqldb.locus_tag2orthogroup_id(server, protein, biodb)
                    orthogroup_list.append(protein_group)
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

            biplot = circos.CircosAccession2biplot(server, db, biodb, reference_records, query_records,
                                                   orthogroup_list, "/home/trestan/Dropbox/dev/django/chlamydia/assets/circos/")

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

