
```

update seqfeature_qualifier_value set value="RB11436b" where seqfeature_id=147658 and value="RB11436";

OK chlamdb-load-gbk.py -g *gbk -d 2019_06_PVC
OK chlamdb-load-orthofinder.py -m Orthogroups.txt -d 2019_06_PVC
OK chlamdb-setup-old_locus-table.py -d 2019_06_PVC
OK chlamdb-setup-genomes-statistics.py -d 2019_06_PVC
OK chlamdb-load-hash2locus.py -u nr_mapping.tab -d 2019_06_PVC


OK chlamdb-load-reference-phylogeny.py -r core_genome_phylogeny.nwk -d 2019_06_PVC -g ../../data/gbk_edited/*gbk

# load interproscan results
OK chlamdb-load-interproscan.py -i *tsv -d 2019_06_PVC -u ../../data/nr_mapping.tab

# add legacy table
OK chlamdb-load-interproscan.py -i *tsv -d 2019_06_PVC -u ../../data/nr_mapping.tab -l

# add add_SP_TM to orthology_ table
OK chlamdb-load-interproscan.py -a -d 2019_06_PVC

# load COG and legacy table
OK chlamdb-load-COG.py -i blast_COG.tab -d 2019_06_PVC -u ../../data/nr_mapping.tab -cc cog_corresp.tab -cl cog_length.tab -l

OK chlamdb-load-KO.py -k chunk*.tab -d 2019_06_PVC -c ../../data/nr_mapping.tab

OK chlamdb-load-uniprot-annotations.py -d 2019_06_PVC -um uniprot_mapping.tab -ud uniprot_data.tab -hm ../../data/nr_mapping.tab


OK chlamdb-load-alignments.py -a *faa -d 2019_06_PVC -c 100
OK chlamdb-load-phylogenies.py -t *nwk -d 2019_06_PVC

OK chlamdb-load-phylogenies-BBH.py -t *nwk -d 2019_06_PVC

OK chlamdb-load-swissprot-homology-search.py -i chunk_.*.tab -d 2019_06_PVC -t -p 2 -l -u ../../data/nr_mapping.tab

OK chlamdb-load-PRIAM.py -i sequenceECs.txt -d 2019_06_PVC -c ../../data/nr_mapping.tab

# comparative tables
OK chlamdb-setup-comparative-tables.py -d 2019_06_chlamydia -o # orthogroup
OK chlamdb-setup-comparative-tables.py -d 2019_06_chlamydia -c # COG
OK chlamdb-setup-comparative-tables.py -d 2019_06_chlamydia -p # pfam
OK chlamdb-setup-comparative-tables.py -d 2019_06_chlamydia -i # interpro
OK chlamdb-setup-comparative-tables.py -d 2019_06_chlamydia -k # ko
OK chlamdb-setup-comparative-tables.py -d 2019_06_chlamydia -e # EC PRIAM

OK chlamdb-setup-linear-taxonomy.py -d 2019_06_PVC -s linear_taxonomy.db
OK chlamdb-setup-gc-content-tables.py -d 2019_06_PVC

OK for i in {1..265}; do echo $i; chlamdb-load-TCDB.py -d 2019_06_PVC -b tcdb -f all.faa -x TCDB_RESULTS_chunk.$i/xml/ -t TCDB_RESULTS_chunk.$i/results.html -c ../../data/nr_mapping.tab; done

OK chlamdb-find-conserved-neighborhood.py -d 2019_06_PVC

```
