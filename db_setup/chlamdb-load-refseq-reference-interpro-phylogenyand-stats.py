#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

def get_rank_summary_statistics(biodb, rank='phylum' ):
    '''

    Get phylogeny from the ncbi taxonomy database given the taxon list in the table pfam.refseq_ref_repres_genomes
    Keep rank phylogeny in the table pfam.phylogeny
    Calculate genome counts for each taxon at the specified rank. Save taxid2count in the table: pfam.<rank>_leaf2n_genomes

    :param rank:
    :return:
    '''
    from ete3 import NCBITaxa, Tree, TextFace,TreeStyle, StackedBarFace
    ncbi = NCBITaxa()
    import os
    import pandas
    sqlpsw = os.environ['SQLPSW']
    
    from sqlalchemy import create_engine
    
    # engine = create_engine("postgresql+psycopg2://scott:tiger@localhost:5432/mydatabase")
    engine = create_engine(f"mysql://root:{sqlpsw}@127.0.0.1/{biodb}")
    engine_conn = engine.connect()
    conn_mysql = engine.raw_connection()   
    cursor = conn_mysql.cursor()




    # create new tables
    sql1 = 'create table if not exists refseq_ref_repres_genomes_phylogeny (dataset varchar(400), tax_rank varchar(400), phylogeny TEXT)'
    sql2 = 'CREATE table if not exists refseq_ref_repres_genomes_leaf2n_genomes(dataset varchar(200), tax_rank varchar(200),taxon_id INT, n_genomes INT)'
    sql3 = 'CREATE table if not exists refseq_ref_repres_genomes_taxid2label(dataset varchar(200), tax_rank varchar(200), taxon_id INT, scientific_name TEXT, node_rank TEXT)'
    cursor.execute(sql1,)
    cursor.execute(sql2,)
    cursor.execute(sql3,)
    conn_mysql.commit()
    
    # retrieve taxid list
    sql_taxid_list = 'select species_taxid from refseq_ref_repres_genomes where keep_interpro=1'
    cursor.execute(sql_taxid_list,)
    
    taxid_list = [i[0] for i in cursor.fetchall()]

    tree = ncbi.get_topology(taxid_list, rank_limit=rank)

    taxon_id_list = [int(i.name) for i in tree.traverse("postorder")]
    
    taxon_id2scientific_name = ncbi.get_taxid_translator(taxon_id_list)


    taxon_id2rank = {}
    for taxon in taxon_id2scientific_name:
        ranks = ncbi.get_rank([taxon])

        try:
            r = ranks[max(ranks.keys())]
        except:
            r = '-'
        taxon_id2rank[taxon] = r

    # save mapping of taxid to phylum level ranks
    for taxon in taxon_id2scientific_name:
        sql = f'insert into refseq_ref_repres_genomes_taxid2label values("interpro", "{rank}", {taxon}, "{taxon_id2scientific_name[taxon]}", "{taxon_id2rank[taxon]}")'
        cursor.execute(sql,)
    conn_mysql.commit()

    # keep some non phylum classification
    # monstly major eukaryote clades
    collapse = ['Opisthokonta', 
                'Alveolata',
                'Amoebozoa',
                'Stramenopiles',
                'Viridiplantae',
                'Rhodophyta', 
                'Trypanosomatidae', 
                'Viruses',
                'unclassified Bacteria', 
                'Leptospiraceae', 
                'unclassified Gammaproteobacteria',
                'unclassified Alphaproteobacteria', 
                'unclassified Epsilonproteobacteria',
                'unclassified Deltaproteobacteria', 
                'unclassified Cyanobacteria (miscellaneous)',
                'unclassified Firmicutes sensu stricto', 
                'unclassified Actinobacteria (class) (miscellaneous)',
                'unclassified Tissierellia', 
                'Dehalogenimonas']

    for node in tree.traverse("postorder"):
        name =  taxon_id2scientific_name[int(node.name)]
        to_detach = []
        
        if name in collapse:
            to_detach.extend(node.children)
            print ('ok-------------------', node.name)
        for n in to_detach:
            n.detach()
    leaves_list = [i.name for i in tree.iter_leaves()]
    leaf_taxon2n_species= {}
    leaf_taxon2n_species_with_domain = {}
    for leaf_taxon in leaves_list:
        leaf_taxon2n_species[leaf_taxon] = 0
        leaf_taxon2n_species_with_domain[leaf_taxon] = 0
        for taxon in taxid_list:
            lineage = ncbi.get_lineage(taxon)
            if int(leaf_taxon) in lineage:
                leaf_taxon2n_species[leaf_taxon]+=1
                #if taxon in taxid_with_domain_list:
                #    leaf_taxon2n_species_with_domain[leaf_taxon]+=1
                
    for leaf_taxon in leaf_taxon2n_species:
        sql = f'insert into refseq_ref_repres_genomes_leaf2n_genomes values("interpro", "{rank}", {leaf_taxon}, "{leaf_taxon2n_species[leaf_taxon]}")'
        cursor.execute(sql,)
    conn_mysql.commit()

    sql = f'insert into refseq_ref_repres_genomes_phylogeny values("interpro","{rank}", "{tree.write(format=1)}")'
    cursor.execute(sql,)
    conn_mysql.commit()


def get_ip_entry_stat(biodb):
    '''
    Combine data of multiple signatures into one IP entry
    Avoid conflicting results where some signature of the same IP where extremely rare and other very frequent
    '''
        
    import os
    import pandas
    sqlpsw = os.environ['SQLPSW']
    
    from sqlalchemy import create_engine
    
    # engine = create_engine("postgresql+psycopg2://scott:tiger@localhost:5432/mydatabase")
    engine = create_engine(f"mysql://root:{sqlpsw}@127.0.0.1/{biodb}")
    engine_conn = engine.connect()
    conn_mysql = engine.raw_connection()   
    cursor = conn_mysql.cursor()
    
    sql = '''
    select distinct t1.interpro_id,t2.name from refseq_ref_repres_genomes_interpro_entries t1 inner join interpro_entry t2 on t1.interpro_id=t2.interpro_id;
    '''
    df_interpro = pandas.read_sql_query(sql, engine)
    
    sql = '''create table refseq_ref_repres_genomes_interpro_entries_freq_v3 
            (interpro_id INT, 
            eukaryote INT, bacteria INT, archae INT, virus INT, total INT, 
            p_eukaryote FLOAT, p_bacteria FLOAT, p_archae FLOAT, p_virus FLOAT,
            s_eukaryote INT, s_bacteria INT, s_archae INT, s_virus INT, 
            s_p_eukaryote FLOAT, s_p_bacteria FLOAT, s_p_archae FLOAT, s_p_virus FLOAT       
            )'''
    cursor.execute(sql)
    conn_mysql.commit()
    
    sql = '''select superkingdom,count(*) as n from (select distinct superkingdom,species from refseq_ref_repres_genomes t1 
    inner join blastnr_blastnr_taxonomy t2 on t1.species_taxid=t2.taxon_id) A group by superkingdom;'''
    
    cursor.execute(sql)
    superk2n_species = {i[0]:int(i[1]) for i in cursor.fetchall()}
    
    
    data_stats = []
    for n,interpro_entry in df_interpro.iterrows():
        '''
        if n == 0:
            continue
        if n == 5:
            break
        '''
        if n % 100 == 0:
            print (f"{n} / {len(df_interpro)}")
        
        
        #################
        # total number of hits
        #################
        counts = {"Eukaryota": 0, "Bacteria": 0, "Archaea": 0, "Viruses": 0}
        sql = f"""select A.superkingdom,count(*) from (select distinct t4.interpro_id,t1.accession,superkingdom from refseq_ref_repres_genomes_interpro_annots t1 
              inner join refseq_ref_repres_genomes_interpro_entries t4 on t1.signature_id=t4.signature_id
              inner join refseq_ref_repres_genomes t2 on t1.assembly_id=t2.assembly_id 
              inner join blastnr_blastnr_taxonomy t3 on t2.species_taxid=t3.taxon_id 
              where t4.interpro_id={interpro_entry.interpro_id} ) A group by A.superkingdom;
              """
        cursor.execute(sql)
        data = cursor.fetchall()
        for kingdom, count in data:
            counts[kingdom] = int(count)
        
        total = sum(counts.values())
        if total == 0:
            continue
        percentages = {superkingdom:count/total * 100 for superkingdom, count in counts.items()}
        #################
        # Number/percentage of species
        #################
        
        n_species = {"Eukaryota": 0, "Bacteria": 0, "Archaea": 0, "Viruses": 0}
        sql = f"""select superkingdom, count(*) as n from (select distinct t3.superkingdom,t3.species from refseq_ref_repres_genomes_interpro_annots t1 
              inner join refseq_ref_repres_genomes_interpro_entries t4 on t1.signature_id=t4.signature_id
              inner join refseq_ref_repres_genomes t2 on t1.assembly_id=t2.assembly_id 
              inner join blastnr_blastnr_taxonomy t3 on t2.species_taxid=t3.taxon_id 
              where t4.interpro_id={interpro_entry.interpro_id}) A group by A.superkingdom;
              """
        cursor.execute(sql)
        data = cursor.fetchall()
        for kingdom, count in data:
            n_species[kingdom] = int(count)
        
        # percentage of all species encoding the domain for each superingdom
        species_percentages = {superkingdom:count/superk2n_species[superkingdom] * 100 for superkingdom, count in n_species.items()}
        
        data_stats.append([interpro_entry["interpro_id"],
                           counts["Eukaryota"], 
                           counts["Bacteria"],
                           counts["Archaea"],
                           counts["Viruses"],
                           total,
                           percentages["Eukaryota"], 
                           percentages["Bacteria"],
                           percentages["Archaea"],
                           percentages["Viruses"],
                           n_species["Eukaryota"], 
                           n_species["Bacteria"],
                           n_species["Archaea"],
                           n_species["Viruses"],
                           species_percentages["Eukaryota"], 
                           species_percentages["Bacteria"],
                           species_percentages["Archaea"],
                           species_percentages["Viruses"],
                     
                     ])

        #if n == 10:
        #    break
    df = pandas.DataFrame(data_stats, columns=["interpro_id", "eukaryote", "bacteria", "archae", "virus", "total", "p_eukaryote", "p_bacteria", "p_archae", "p_virus", "s_eukaryote", 
                                               "s_bacteria", "s_archae", "s_virus", "s_p_eukaryote", "s_p_bacteria", "s_p_archae", "s_p_virus"])
    df.to_sql("refseq_ref_repres_genomes_interpro_entries_freq_v3", 
                                        engine, 
                                        index=False, 
                                        if_exists="append")
        
    



def get_signature_stats(biodb):
    import os
    import pandas
    sqlpsw = os.environ['SQLPSW']
    
    from sqlalchemy import create_engine
    
    # engine = create_engine("postgresql+psycopg2://scott:tiger@localhost:5432/mydatabase")
    engine = create_engine(f"mysql://root:{sqlpsw}@127.0.0.1/{biodb}")
    engine_conn = engine.connect()
    conn_mysql = engine.raw_connection()   
    cursor = conn_mysql.cursor()
    
    df_signatures = pandas.read_sql_query('select accession as entry_accession,signature_id from refseq_ref_repres_genomes_interpro_entries', engine)
    
    sql = '''create table refseq_ref_repres_genomes_interpro_entries_freq_v2 
            (signature_id INT, 
            eukaryote INT, bacteria INT, archae INT, virus INT, total INT, 
            p_eukaryote FLOAT, p_bacteria FLOAT, p_archae FLOAT, p_virus FLOAT,
            s_eukaryote INT, s_bacteria INT, s_archae INT, s_virus INT, 
            s_p_eukaryote FLOAT, s_p_bacteria FLOAT, s_p_archae FLOAT, s_p_virus FLOAT       
            )'''
    cursor.execute(sql)
    conn_mysql.commit()
    
    sql = 'select superkingdom,count(*) as n from (select distinct superkingdom,species from refseq_ref_repres_genomes t1 inner join blastnr_blastnr_taxonomy t2 on t1.species_taxid=t2.taxon_id) A group by superkingdom;'
    cursor.execute(sql)
    superk2n_species = {i[0]:int(i[1]) for i in cursor.fetchall()}
    
    
    data_stats = []
    for n,signature in df_signatures.iterrows():
        #if n == 0:
        #    continue
        if n % 100 == 0:
            print (f"{n} / {len(df_signatures)}")
        
        
        #################
        # total number of hits
        #################
        counts = {"Eukaryota": 0, "Bacteria": 0, "Archaea": 0, "Viruses": 0}
        sql = f"""select t3.superkingdom,count(*) from refseq_ref_repres_genomes_interpro_annots t1 
              inner join refseq_ref_repres_genomes t2 on t1.assembly_id=t2.assembly_id 
              inner join blastnr_blastnr_taxonomy t3 on t2.species_taxid=t3.taxon_id 
              where signature_id={signature.signature_id} group by t3.superkingdom;
              """
        cursor.execute(sql)
        data = cursor.fetchall()
        for kingdom, count in data:
            counts[kingdom] = int(count)
        
        total = sum(counts.values())
        if total == 0:
            continue
        percentages = {superkingdom:count/total * 100 for superkingdom, count in counts.items()}
        #################
        # Number/percentage of species
        #################
        
        n_species = {"Eukaryota": 0, "Bacteria": 0, "Archaea": 0, "Viruses": 0}
        sql = f"""select superkingdom, count(*) as n from (select distinct t3.superkingdom,t3.species from refseq_ref_repres_genomes_interpro_annots t1 
              inner join refseq_ref_repres_genomes t2 on t1.assembly_id=t2.assembly_id 
              inner join blastnr_blastnr_taxonomy t3 on t2.species_taxid=t3.taxon_id 
              where signature_id={signature.signature_id}) A group by A.superkingdom;
              """
        cursor.execute(sql)
        data = cursor.fetchall()
        for kingdom, count in data:
            n_species[kingdom] = int(count)
        
        # percentage of all species encoding the domain for each superingdom
        species_percentages = {superkingdom:count/superk2n_species[superkingdom] * 100 for superkingdom, count in n_species.items()}
        
        data_stats.append([signature["signature_id"],
                           counts["Eukaryota"], 
                           counts["Bacteria"],
                           counts["Archaea"],
                           counts["Viruses"],
                           total,
                           percentages["Eukaryota"], 
                           percentages["Bacteria"],
                           percentages["Archaea"],
                           percentages["Viruses"],
                           n_species["Eukaryota"], 
                           n_species["Bacteria"],
                           n_species["Archaea"],
                           n_species["Viruses"],
                           species_percentages["Eukaryota"], 
                           species_percentages["Bacteria"],
                           species_percentages["Archaea"],
                           species_percentages["Viruses"],
                     
                     ])

        #if n == 10:
        #    break
    df = pandas.DataFrame(data_stats, columns=["signature_id", "eukaryote", "bacteria", "archae", "virus", "total", "p_eukaryote", "p_bacteria", "p_archae", "p_virus", "s_eukaryote", 
                                               "s_bacteria", "s_archae", "s_virus", "s_p_eukaryote", "s_p_bacteria", "s_p_archae", "s_p_virus"])
    df.to_sql("refseq_ref_repres_genomes_interpro_entries_freq_v2", 
                                        engine, 
                                        index=False, 
                                        if_exists="append")
        
    
    

    

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", '--mysql_database', type=str, help="Biosql biodatabase name")

    
    args = parser.parse_args()
    
    #get_rank_summary_statistics(args.mysql_database, "phylum")
    #get_signature_stats(args.mysql_database)
    get_ip_entry_stat(args.mysql_database)