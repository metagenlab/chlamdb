
def get_hash2locus_list(hash2locus_tag):
    print("parsing hash2locus_tag")
    hash2locus_list = {}
    with open(hash2locus_tag, 'r') as f:
        for row in f:
            locus, hash, genome = row.rstrip().split('\t')
            if hash not in hash2locus_list:
                hash2locus_list[hash] = [locus]
            else:
                hash2locus_list[hash].append(locus)
    return hash2locus_list
