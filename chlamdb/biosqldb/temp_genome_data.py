#!/usr/bin/python

import biosql_own_sql_tables

biosql_own_sql_tables.collect_genome_statistics('corynedb5')
#biosql_own_sql_tables.collect_genome_statistics('2017_04_19_chlamydiae_taxo')
#biosql_own_sql_tables.collect_genome_statistics('2017_06_29b_motile_chlamydiae')
#biosql_own_sql_tables.collect_genome_statistics('rhabdochlamydia_07_16')
#biosql_own_sql_tables.collect_genome_statistics('2017_03_30_kcosson')

biosql_own_sql_tables.collect_genome_statistics('chlamydia_motile', '/home/trestan/work/dev/metagenlab/genome_pairwise_aa_identity/biosqldb.db')

#mat1, mat2 = biosql_own_sql_tables.get_comparative_subtable("chlamydia_04_16", "interpro", "id", ["46","48","49"],["52"],0.5)
#print dict((mat2.iloc[0:12000,0:10]> 0).sum(axis=1))
#print dict((mat2 > 0).sum(axis=1))