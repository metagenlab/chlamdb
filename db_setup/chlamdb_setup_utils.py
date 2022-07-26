
def get_hash2locus_list(hash2locus_tag, as_df=False):
    import pandas 
    df = pandas.read_csv(hash2locus_tag, sep="\t", names=["locus_tag","hash", "genome"])
    
    if as_df:
        return df
    else:
        hash2locus_list = {}
        for n,row in df.iterrows():
            if row.hash not in hash2locus_list:
                hash2locus_list[row.hash] = [row.locus_tag]
            else:
                hash2locus_list[row.hash].append(row.locus_tag)
        return hash2locus_list
