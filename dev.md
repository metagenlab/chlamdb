
# CHLAMDB setup

## setup priority 1

- [X] setup biosqldb schema
- [ ] setup kegg tables (ENZYME database)
  - [ ] separate script into multiple scripts and deal with incomplete tables
  - [X] move to separate script (not in load script anymore)
  - [X] ko2annotation table
    - [ ] redundancy with ko2module and ko2pathway, can be simplified and accelerated
  - [X] setup enzyme_dat and enzyme tables
  - [X] setup pathway table
  - [X] setup ko2pathway
  - [X] setup module table
  - [X] setup ko2module
- [X] setup COG tables (COG database)
- [X] setup linear taxonomy (with taxid for each rank see virulence db setup)
- [ ] setup TCDB tables
  - [X] download faa: setup_transporter table, tc_table and uniprot_table (TODO: rename table)
  - [ ] update table organization (dev), not all entries have uniprot accessions

## setup priority 2

- [ ] switch to new COG tables
- [ ] switch to new KEGG tables
- [ ] setup NOG_table_v5 and NOG_members_v5
- [ ] setup interpro master table
- [ ] setup pfam master table
- [ ] load DOOR2 data
    - [ ] accession table
    - [ ] operons data
- [ ] load RBBH data (for identity distribution plots)
- [ ] get effectiveT3 "eukaryote" domains
- [ ] setup alignnments from nucletide sequence with genoplotR (see locus2genoplotr.py)

### loading of annotation tables

- [X] load genomes
  - [X] setup features table
  - [X] setup seqfeature_id2locus_table
- [X] load orthofinder results  
- [X] load interproscan results
  - [X] add TM and SP columns to orthology_detail legacy table
- [X] load orthogroup alignments
- [X] get identity closest homolog table
- [X] get average identity table
- [X] setup genome statistics table
- [X] load COG hits
- [X] load KO hits
  - [X] load legacy tables
- [X] setup pairwise BBH identity tables
- [X] lead orthogroup alignments(identity matrices)
- [X] load orthogroup phylogenies
- [X] load orthogroup BBH phylogenies
- [X] load reference phylogeny
- [X] get genome table (homepage)
- [X] get conserved neighborhood
- [X] load uniprot annotations
- [X] load blast swissprot results
  - [X] download taxonomy-description information(s)
- [X] load blast refseq results
    - [X] load refseq taxonomy table
- [X] load TCDB annotations
- [X] get phylo profile
  - [ ] setup core_orthogroups_identity_msa_
- [X] legacy COG table
- [X] legacy locus2EC table
- [X] legacy PFAM table
- [X] setup blast databases
- [X] phylogenetic profiles
- [X] use celery for circos_main view
- [X] load PMID string mapping
- [X] load checkM
- [X] update PMID
- [X] load idmapping data
  - [X] from uniprot idmapping
  - [X] from RefSeq (locus tag and protein accession)
- [X] load T3SS effector predictions
  - [X] T3_MM
  - [X] deep_T3
  - [X] effective
  - [X] BPBAac
- [X] load cross-references
- [X] load psortdb results
- [X] load pdb best hits (with score)
- [ ] load synonymous table 
    - [X] from uniprot idmapping
    - [x] RefSeq locus_tag/protein id
    - [ ] uniparc cross references
- [ ] load GOA annotation
- [ ] similarity networks dervied from patristic distances (from RefSeq phylogeny?) ==> better than identity, score or e-value


# Web Interface

## TODO

### general priority 1

- [X] show confidence scores for PDB, KEGG, COG,... (identity, score, evalue,...)
- [X] keep comparative tables in memory
- [X] blast multiple proteins
- [X] blast link only when blasting locus databases
- [X] update browse genome view 
- [X] add formatdb ffn
- [X] database all for tblastn
- [X] add effector prediction to locus page
- [ ] switch to celery task
    - [X] plot region ==> celery task
    - [X] extract orthogroup ==> celery task
    - [X] plot phylogeny
    - [X] plot orthogroup/COG/KO,... heatmap ==> celery task
    - [X] extract interpro
    - [ ] BLAST ==> celery task
    - [ ] venn orthogroups
    - [ ] venn interpro
    - [ ] fam (crash with interpro entries with very large number of proteins) 
- [X] improve orthogroup table
    - [X] length distribution
    - [X] domain organization
    - [X] mapping to uniprot (n mapped, n reviewed,...)

- [X] 2019_06_PVC: add MSA faa to assets
- [X] update search using synonymous table single vs multiple matches
- [X] update to boolean search
- [X] deal with COG, PFAM,... missing from the database
- [X] check and update module, pathway, fam profile,... profile figure size
- [X] add download newick tree on pfam, TM tabs and phylogeny with BBH
- [ ] add uniprot proteome column (+ percent overlap) to genome table
- [ ] add explanations for hydropathy plot
- [ ] improve integration of transporters_family 
- [ ] orthogroup Venn: use consensus annotation
- [ ] add a vew "compare_cog" (similar to compare Pfam)
- [ ] add a vew "compare_interpro" (similar to compare Pfam)

### download page

- [ ] download page - bulk download all genome or for specific genomes 
  - [ ] KO annotation 
  - [ ] COG annotation 
  - [ ] TCDB annotation
  - [ ] interproscan
  - [ ] pfam only 
  - [ ] all phylogenies
  - [ ] all orthogroup fasta
  - [ ] orthology table 
  - [ ] all alignments
  - [ ] ...

### curated taxnonomy

- [X] setup a new view for manual curation of the taxonomy
- [ ] switch to curated taxnonomy (ref phylogeny, locus pages,...) 
- [ ] chose a reference strain for each species?
- [ ] search: show reference species data first 
 
### species tree display

- [ ] species tree based on curated taxnonomy
- [ ] phylo profile: display species tree by default 
    - [ ] show average or median identity as compared to reference locus
    - [ ] OR only show identity of the reference strain (rely on curated taxnonomy)

### effectors/euk like domains 

- [ ] plot distribution of euk proportion (how is the distribution?). See where to put the cutoff.
- [X] new view to compare data from the different T3SS effector predictiors: =effector_predictions=, 
- [ ] see also *venn_candidate_effectors* and *effector_pred*

### user interface

- [ ] gene sets: table allowing redordering (drag and drop)
- [ ] plot distribution normalized identity between species (rate of evolution)
- [ ] species tree (collapsed)
- [ ] paralogs tab?
- [ ] plot n species specific groups (collapse tree by species)
- [ ] plot n genus specific groups,...
- [ ] circos taxonomy
- [ ] overview number of different KO/group, COG/group,...
- [ ] overview taxonomy top 200 refseq hits
- [ ] blast page: highlight best non self phylum hit
- [ ] blast celery
- [ ] cross references
- [ ] go terms
- [ ] download page
- [ ] krona like plots for blast
- [ ] kronal like plots for 1 blast hit, 2 blast hit,... (choice) 

### ideas


- [ ] HAMAP rules as SPARQL A portable annotation pipeline for genomes and proteomes (https://www.biorxiv.org/content/10.1101/615294v1) 
- [ ] unifire annotation (UniRule and SAAS from interproscan results): https://gitlab.ebi.ac.uk/uniprot-public/unifire
  - [ ] https://www.ebi.ac.uk/training/online/sites/ebi.ac.uk.training.online/files/UniFIRE-URML.pdf
- [ ] kegg pathway => link show on circos
- [ ] kegg module => link show on circos
- [ ] integrate genome properties: table with property2steps avec IPP accession
- [ ] eficaz
- [ ] enzyme hierarchy
- [ ] interpro hierarchy 
- [ ] eggnog comparative table 
- [ ] cdd comparative table
- [ ] update to new eggnog version (and string)
- [ ] add uniprot kewords to locus page
- [ ] get species table 
  - [ ] display on homepage
  - [ ] use it for phylogenetic profiling rather than all strains (inclomplete genomes)?
- [ ] deal with search for KO, KEGG, IP absent from genomes included in the database
- [ ] integration of swissprot keywords (possibility to click on it and get complete list of prot, decsription,...)
- [ ] add tcdb classification to "fam" (annotation, phylogenetic profile,...). Include all classification levels (superfamilies,...)
- [ ] idem with EC classification system
- [ ] search bar: add option to search for TCDB accessions
- [ ] orthoinspector
