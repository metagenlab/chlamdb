#!/usr/bin/python

import biosql_own_sql_tables

biosql_own_sql_tables.collect_genome_statistics('chlamydia_04_16')
#biosql_own_sql_tables.collect_genome_statistics('2017_04_19_chlamydiae_taxo')
#biosql_own_sql_tables.collect_genome_statistics('2017_04_chlamydia_zurich')
#biosql_own_sql_tables.collect_genome_statistics('rhabdochlamydia_07_16')

#mat1, mat2 = biosql_own_sql_tables.get_comparative_subtable("chlamydia_04_16", "interpro", "id", ["46","48","49"],["52"],0.5)
#print dict((mat2.iloc[0:12000,0:10]> 0).sum(axis=1))
#print dict((mat2 > 0).sum(axis=1))