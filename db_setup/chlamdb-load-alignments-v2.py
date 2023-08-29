#!/usr/bin/env python

import pandas
import os 
from sqlalchemy import create_engine

def create_table(biodb):
    sqlpsw = os.environ['SQLPSW']
    engine = create_engine(f"mysql://root:{sqlpsw}@127.0.0.1/{biodb}")
    conn_mysql = engine.connect()
    sql = 'CREATE TABLE if not exists orthology_identity_v2 (' \
                ' seqfeature_id_a INTEGER,' \
                ' seqfeature_id_b INTEGER,' \
                ' pident FLOAT(5) ,' \
                ' length INTEGER)'
    conn_mysql.execute(sql)
    
    sql = 'CREATE TABLE if not exists orthology_average_identity (orthogroup_id INTEGER, identity float(5))'
    conn_mysql.execute(sql)
    
    conn_mysql.close()

def create_indexes(biodb):
    sqlpsw = os.environ['SQLPSW']
    engine = create_engine(f"mysql://root:{sqlpsw}@127.0.0.1/{biodb}")
    conn_mysql = engine.connect()
    
    sql1 = 'create index oidvsa on orthology_identity_v2(seqfeature_id_a)'
    sql2 = 'create index oidvsb on orthology_identity_v2(seqfeature_id_b)'
    sql3 = 'create index oaid on orthology_average_identity(orthogroup_id)'
    conn_mysql.execute(sql1)
    conn_mysql.execute(sql2)
    conn_mysql.execute(sql3)
    conn_mysql.close()


def get_parwise_id(identity_csv_list, 
                   biodb):
    
    sqlpsw = os.environ['SQLPSW']
    engine = create_engine(f"mysql://root:{sqlpsw}@127.0.0.1/{biodb}")
    conn_mysql = engine.connect()
    print("parse identity tables")
    create_table(biodb)
    print("Retrieve locus tag 2 seqfeature id")
    sql = 'select locus_tag,seqfeature_id from custom_tables_locus2seqfeature_id'
    locus2_seqfeature_id_df = pandas.read_sql(sql, conn_mysql).set_index("locus_tag")
    sql = 'select seqfeature_id,orthogroup_id from orthology_seqfeature_id2orthogroup'
    seqfeature_id2orthogroup = pandas.read_sql(sql, conn_mysql).set_index("seqfeature_id")
    sql_template = 'insert into orthology_average_identity values(%s, %s)'
    for n, one_id_file in enumerate(identity_csv_list):
        print(one_id_file)
        if n % 100 == 0:
            print(f"{n} / {len(identity_csv_list)}")
        # KJA58_RS02945,KJA62_RS03385,73.76760563380282,568
        df_identity = pandas.read_csv(one_id_file, names = ["locus_a", "locus_b", "pident", "length"])
        df_identity["pident"] = df_identity["pident"].round(2)

        #print("Merge with seafeature ids df")
        df_identity = df_identity.set_index("locus_a").join(locus2_seqfeature_id_df.rename(columns={"seqfeature_id": "seqfeature_id_a"})) #.reset_index(col_fill="locus_a")
        df_identity.index.name = "locus_a"
        df_identity = df_identity.reset_index().set_index("locus_b").join(locus2_seqfeature_id_df.rename(columns={"seqfeature_id": "seqfeature_id_b"})) #.reset_index(col_fill="locus_b")
        df_identity.index.name = "locus_b"
        
        #print(f"Inserting {n} \t {len(df_identity)} entries...")
        df_identity = df_identity.reset_index()[["seqfeature_id_a", "seqfeature_id_b", "pident", "length"]]

        df_identity.to_sql("orthology_identity_v2", conn_mysql, if_exists="append", index=False, chunksize=50000)
        
        orthogroup_id = seqfeature_id2orthogroup.loc[df_identity.loc[0].seqfeature_id_a].orthogroup_id
        average_identity = df_identity["pident"].mean()
        conn_mysql.execute(sql_template % (orthogroup_id, average_identity))


def get_average_id_table(db_name):
    
    sql = f'select orthogroup from comparative_tables_orthology'
                
    groups = [i[0] for i in self.server.adaptor.execute_and_fetchall(sql, )]
    #print len(groups)

    for group in groups:
        id_values = self.get_orthogroup_identity_table(group)
        if len(id_values) > 0:
            av_id = round(np.mean(id_values), 2)
            
            sql = 'insert into orthology_average_identity values ("%s", %s)' % (group, av_id)

            self.server.adaptor.execute(sql)
    self.server.commit()


    def get_orthogroup_identity_table(self, orthogroup):
        import os

        sql = 'SELECT identity FROM orthology_identity where orthogroup="%s"' % orthogroup

        values = [float(i[0]) for i in self.server.adaptor.execute_and_fetchall(sql,)]
        return values



if __name__ == '__main__':
    import argparse
    from chlamdb.biosqldb import manipulate_biosqldb
    from chlamdb.biosqldb import manipulate_biosqldb
    
    parser = argparse.ArgumentParser()

    parser.add_argument("-i",'--csv_files', type=str, help="input identitx csv files", nargs='+')
    parser.add_argument("-d",'--db_name', type=str, help="database_name")
    parser.add_argument("-a",'--average', action="store_true", help="Add average identity table")
    parser.add_argument("-x",'--index', action="store_true", help="Index table")   

    args = parser.parse_args()

    if args.index:
        create_indexes(args.db_name)
    
    else:
        get_parwise_id(args.csv_files, args.db_name)  

    if args.average:
        get_average_id_table(args.db_name)
    
    manipulate_biosqldb.update_config_table(args.db_name, "orthogroup_alignments")
    
