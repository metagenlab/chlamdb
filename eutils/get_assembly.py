#!/usr/bin/env python

from Bio import Entrez, SeqIO
import eutils
Entrez.email = "trestan.pillonel@unil.ch"

def gbk2faa(seq_record, outname):
    output_handle = open(outname, "w")
    for seq_feature in seq_record.features :
        if seq_feature.type=="CDS" :
            #print seq_feature
            assert len(seq_feature.qualifiers['translation'])==1
            # gi|83716028|ref|YP_443839.1| matrix protein [Avian metapneumovirus]
            try:
                output_handle.write(">gi|%s|ref|%s| %s [%s]\n%s\n" % (
                        seq_feature.qualifiers["db_xref"][0].split(":")[1],
                        seq_feature.qualifiers["protein_id"][0],
                        seq_feature.qualifiers["note"][0],
                        seq_record.description,
                        seq_feature.qualifiers['translation'][0]))
            except:
                output_handle.write(">gi|%s|ref|%s| [%s]\n%s\n" % (
                        seq_feature.qualifiers["db_xref"][0].split(":")[1],
                        seq_feature.qualifiers["protein_id"][0],
                        #seq_feature.qualifiers["note"][0],
                        seq_record.description,
                        seq_feature.qualifiers['translation'][0]))

def gbk2ffn(seq_record, outname):
    output_handle = open(outname, "w")
    for seq_feature in seq_record.features :
        if seq_feature.type=="CDS" :
            #print seq_feature
            assert len(seq_feature.qualifiers['translation'])==1
            # gi|83716028|ref|YP_443839.1| matrix protein [Avian metapneumovirus]
            try:
                output_handle.write(">gi|%s|ref|%s| %s [%s]\n%s\n" % (
                        seq_feature.qualifiers["db_xref"][0].split(":")[1],
                        seq_record.id,
                        seq_feature.qualifiers["note"][0],
                        seq_record.description,
                        seq_feature.extract(seq_record.seq)))
            except:
                output_handle.write(">gi|%s|ref|%s| [%s]\n%s\n" % (
                        seq_feature.qualifiers["db_xref"][0].split(":")[1],
                        seq_record.id,
                        seq_record.description,
                        seq_feature.extract(seq_record.seq)))



def download_one_wgs(wgs_link):
        handle = Entrez.elink(dbfrom="nuccore", db="nuccore", id=wgs_link)
        record = Entrez.read(handle)

        # get the list of link
        sublinks = [link["Id"] for link in record[0]["LinkSetDb"][1]["Link"]]
        #print "sublinks", sublinks

        if len(sublinks) > 500:
            print "More than 500 contigs, aborting download for id %s" % wgs_link
            return None
        
        # for each sublink (contig, scaffold,...), get record and append the sequence to a single file 
        output_handle = open("%s.fasta" % wgs_link, "a")
        no_sequences = False

        
        for seq_link in sublinks:
            if no_sequences:
                print "No sequences fo record %s" % wgs_link
                output_handle.close()
                
                import os
                os.remove("%s.fasta" % wgs_link)
                
                break
            handle = Entrez.efetch(db="nucleotide", id=seq_link, rettype="gb", retmode="text")
            seq_records = list(SeqIO.parse(handle, "genbank"))
            # in case multiple record for a single link (shouldn't append???)
            for record in seq_records:
                print record.name
                print record.description
                if record.seq.count("N") == len(record.seq):
                    no_sequences = True
                    break
                else:
                    SeqIO.write(record, output_handle, "fasta")
        output_handle.close()    

def get_wgs_links(one_species_link):

        # get all WGS linked to this species
        handle = Entrez.elink(dbfrom="genome", db="nuccore", id=one_species_link, term="wgs[prop]")
        record = Entrez.read(handle)
        if len(record[0]["LinkSetDb"][0]["Link"]) == 0:
            print "No WGS genome seq for %s" % one_genome_id
            return False
        else:
            linked = [link["Id"] for link in record[0]["LinkSetDb"][0]["Link"]]
            #print "WGS genome(s):", linked
            return linked
            
                

def multiple_wgs_links(ncbi_taxon):
    #handle = Entrez.esearch(db="genome", term="klebsiella+pneumoniae[orgn]")
    #txid570[Organism:exp]
    handle = Entrez.esearch(db="genome", term="txid%s[Organism:exp]" % ncbi_taxon)
    record = Entrez.read(handle)

    # get genome overview id
    genome_id_list = record["IdList"]

    # elink.fcgi?dbfrom=genome&db=nuccore&id=1076

    # get whole genomes only
    genome_record_id_list = []
    for one_genome_id in genome_id_list:
        one_genome_ids = get_wgs_links(one_genome_id)
        genome_record_id_list += one_genome_ids
        
    return genome_record_id_list
        
        
def get_complete_genomes_data(ncbi_taxon):

    import eutils

    #handle = Entrez.esearch(db="genome", term="klebsiella+pneumoniae[orgn]")
    #txid570[Organism:exp]
    handle = Entrez.esearch(db="genome", term="txid%s[Organism:exp]" % ncbi_taxon)
    record = Entrez.read(handle)
    print record
    # get genome overview id
    genome_id_list = record["IdList"]

    # elink.fcgi?dbfrom=genome&db=nuccore&id=1076

    # get whole genomes only
    genome_record_id_list = []
    for one_genome_id in genome_id_list:
        print "considering genome ID", one_genome_id
        # srcdb+ddbj/embl/genbank[prop] AND 
        handle = Entrez.elink(dbfrom="genome", db="nuccore", id=one_genome_id, term="gene+in+genomic[prop] OR gene+in+chromosome[prop]")


        record = Entrez.read(handle)
        print "all links", record[0]["LinkSetDb"][0]["Link"]
        
        if len(record[0]["LinkSetDb"][0]["Link"]) == 0:
            print "No whole genome seq for %s" % one_genome_id

        else:
            print "Whole genome data for %s " % one_genome_id
            print record


            
            """
            handle_plasmids = Entrez.elink(dbfrom="genome", db="nuccore", id=one_genome_id, term="srcdb+ddbj/embl/genbank[prop] AND gene+in+plasmid[prop]")
            record_plasmids = Entrez.read(handle_plasmids)
                

            if len(record_plasmids[0]["LinkSetDb"][0]["Link"]) == 0:
                print "No plasmid seq for %s" % one_genome_id
            else:
                linked_plasmids = [link["Id"] for link in record_plasmids[0]["LinkSetDb"][0]["Link"]]
                print "Plasmid(s):", linked_plasmids
                genome_record_id_list += linked_plasmids
            """
            
            linked = [link["Id"] for link in record[0]["LinkSetDb"][0]["Link"]]
            print "Complete genome(s):", linked
            genome_record_id_list += linked


    for record_id in genome_record_id_list:
        print "ID:", record_id

        # get assembly
        handle_assembly = Entrez.elink(dbfrom="nuccore", db="assembly", id=record_id)
        record_assembly = Entrez.read(handle_assembly)
        try:
            assembly_link = record_assembly[0]["LinkSetDb"][0]["Link"][0]["Id"]
            print "assembly link", assembly_link
            # assembly 2 genome + plasmids
            handle_sequences = Entrez.elink(dbfrom="assembly", db="nuccore", id=assembly_link, term="srcdb+ddbj/embl/genbank[prop]")
            record_sequences =  Entrez.read(handle_sequences)
            wgs = False
            for i in range(0, len(record_sequences[0]["LinkSetDb"])):
                if record_sequences[0]["LinkSetDb"][i]['LinkName'] == 'assembly_nuccore_wgsmaster':
                    wgs = True
            if wgs:
                print "wgs link"
                continue
            
            sequences_links = [link["Id"] for link in record_sequences[0]["LinkSetDb"][0]["Link"]]            
        except:
            print "No assembly link, searching bioproject for %s..." % record_id

            # first: check if not WGS:
            handle = Entrez.efetch(db="nucleotide", id=record_id, rettype="gb", retmode="text")
            seq_records = list(SeqIO.parse(handle, "genbank"))


            print "annot", seq_records[0].annotations
            if "wgs" in seq_records[0].annotations or "WGS" in seq_records[0].annotations:
                print "Wgs, not downloading %s" % record_id, seq_records[0].description
                continue                
            elif isinstance(seq_records[0].annotations, dict) and len(seq_records[0].annotations.values()) > 1:
                tag = False
                for i in seq_records[0].annotations.values():
                    print "annot", i
                    if isinstance(i, list):
                        if "wgs" in i or "WGS" in i:
                            print "Wgs, not downloading %s" % record_id, seq_records[0].description
                            tag = True
                    elif isinstance(i, str):
                        if "wgs" in i or "WGS" in i:
                            print "Wgs, not downloading %s" % record_id, seq_records[0].description
                            tag = True
                    elif isinstance(i, int):
                        continue
                    else:
                        print "problem with:", i
                if tag:
                    continue

            handle_bioproject = Entrez.elink(dbfrom="nuccore", db="bioproject", id=record_id)
            record_bioproject = Entrez.read(handle_bioproject)

            # get bioproject link
            bioproject_link = record_bioproject[0]["LinkSetDb"][0]["Link"][0]["Id"]

            # as bioprojects are not linked with genbank, get link(s) to refseq (for both plasmids and chromosomes)
            handle_refseq = Entrez.elink(dbfrom="bioproject", db="nuccore", id=bioproject_link)
            record_refseq =  Entrez.read(handle_refseq)

            refseq_links = [link["Id"] for link in record_refseq[0]["LinkSetDb"][0]["Link"]]

            print "refseqlinks", refseq_links
            # for each refseq id, get link to genbank entry (which contain full record)
            sequences_links = []
            for link in refseq_links:
                handle_genbank = Entrez.elink(dbfrom="nuccore", db="nuccore", id=link, term="srcdb+ddbj/embl/genbank[prop]")
                record_genbank = Entrez.read(handle_genbank)
                sequences_links.append(record_genbank[0]["LinkSetDb"][1]["Link"][0]["Id"])
            



        print "Sequences links:", sequences_links

        if len(sequences_links) == 0:
            print "PROBLEM WITH ASSEMBLY"
            print record_sequences
        #import time
        #time.sleep(4)
        for one_id in sequences_links:
            eutils.get_genomic_data(one_id)





if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-t",'--taxon_id',type=str,help="get wgs link from taxonomic id")
    parser.add_argument("-w",'--wgs',type=str,help="download one wgs link")
    parser.add_argument("-c",'--complete_genomes',type=str,help="taxonomic id (get complete genomes only")

    
    args = parser.parse_args()

    if args.taxon_id:
        print "getting link to wgs records for taxon %s ..." % args.taxon_id
        all_links = multiple_wgs_links(args.taxon_id)
        for i in all_links:
            print i

    if args.wgs:
        print "downloading record %s..." % args.wgs
        download_one_wgs(args.wgs)
   
    if args.complete_genomes:
        print "downloading complete genomes from taxon %s..." % args.complete_genomes
        get_complete_genomes_data(args.complete_genomes)
