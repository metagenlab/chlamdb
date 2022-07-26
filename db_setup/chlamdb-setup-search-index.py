#!/usr/bin/env python


if __name__ == '__main__':
        import argparse
        import logging
        logging.basicConfig(filename="index.log", level=logging.DEBUG)
        from chlamdb.biosqldb import manipulate_biosqldb

        parser = argparse.ArgumentParser()
        parser.add_argument('-d', '--db_name', type=str, help='DB name', default='Biodb name')
        parser.add_argument('-v', '--interpro_version', type=str, help='Interpro version')

        args = parser.parse_args()

        server, db = manipulate_biosqldb.load_db(args.db_name)
        conn = server.adaptor.conn
        cursor = server.adaptor.cursor

        sql1= 'CREATE FULLTEXT INDEX GPF1 ON orthology_detail(gene);'
        sql2= 'CREATE FULLTEXT INDEX GPF2 ON orthology_detail(product);'
        sql3= 'CREATE FULLTEXT INDEX GPF3 ON orthology_detail(organism);'
        sql4= 'CREATE FULLTEXT INDEX GPF4 ON orthology_detail(gene,product,organism);'

        sql5= 'CREATE FULLTEXT INDEX ezf ON enzyme_enzymes_dat(value);'

        sql6= 'CREATE FULLTEXT INDEX koaf ON enzyme_ko_annotation_v1(definition);'
        # filter KO absent from DB
        sql7= 'CREATE FULLTEXT INDEX ko_full ON enzyme_ko_annotation_v1(ko_id,name,definition,EC);'

        sql8= 'CREATE FULLTEXT INDEX ipf ON interpro_signature(signature_description);'
        sql9= 'CREATE FULLTEXT INDEX ipf ON interpro_entry(description);'


        sql10= 'CREATE FULLTEXT INDEX modf ON enzyme_kegg_module_v1(description, module_sub_cat);'

        sql11= 'CREATE FULLTEXT INDEX kegf ON enzyme_kegg_pathway(description, pathway_category);'
        sql12= 'ALTER TABLE string_pmid2data_stringdb ADD FULLTEXT (title,journal,year,authors);'

        logging.info(sql1)
        try: 
            cursor.execute(sql1)
        except Exception as e:
            print(str(e))
            print('Duplicate' in str(e))
            if 'Duplicate' in str(e):
                pass 
            else:
                raise Exception(e)
        logging.info(sql2)
        try:      
            cursor.execute(sql2)
        except Exception as e:
            if 'Duplicate' in str(e):
                pass 
            else:
                raise Exception(e)
        logging.info(sql3)
        try:
                cursor.execute(sql3)
        except Exception as e:
            if 'Duplicate' in str(e):
                pass 
            else:
                raise Exception(e)
        logging.info(sql4)
        try:
            cursor.execute(sql4)
        except Exception as e:
            if 'Duplicate' in str(e):
                pass 
            else:
                raise Exception(e)
        logging.info(sql5)
        try:
            cursor.execute(sql5)
        except Exception as e:
            if 'Duplicate' in str(e):
                pass 
            else:
                raise Exception(e)
        logging.info(sql6)
        try:
            cursor.execute(sql6)
        except Exception as e:
            if 'Duplicate' in str(e):
                pass 
            else:
                raise Exception(e)
        logging.info(sql7)
        try:
            cursor.execute(sql7)
        except Exception as e:
            if 'Duplicate' in str(e):
                pass 
            else:
                raise Exception(e)
        logging.info(sql8)
        try:
            cursor.execute(sql8)
        except Exception as e:
            if 'Duplicate' in str(e):
                pass 
            else:
                raise Exception(e)

        logging.info(sql9)
        try:
            cursor.execute(sql9)
        except Exception as e:
            if 'Duplicate' in str(e):
                pass 
            else:
                raise Exception(e)
        
        logging.info(sql10)
        try:
            cursor.execute(sql10)
        except Exception as e:
            if 'Duplicate' in str(e):
                pass 
            else:
                raise Exception(e)
            
        logging.info(sql11)
        try:
            cursor.execute(sql1)
        except Exception as e:
            if 'Duplicate' in str(e):
                pass 
            else:
                raise Exception(e)
        logging.info(sql12)
        try:
            cursor.execute(sql12)
        except Exception as e:
            if 'Duplicate' in str(e):
                pass 
            else:
                raise Exception(e)
        conn.commit()
        