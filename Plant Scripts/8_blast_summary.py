## script to extract N match summary from the raw blast output
import re

max_matches = 1

# load blast data
handle = open('results/denovo_blast_plants.txt')
data = handle.readlines()
data = [x.rstrip() for x in data]
handle.close()

# PARSE the blast file 
# - find query lines to identify starts of blast blocks
query_index = [[i, x] for i, x in enumerate(data) if x.startswith('Query= ')]
query_index, query_name = zip(*query_index)
query_name = [x.replace('Query= ','') for x in query_name]
blocks = zip(query_index, query_index[1:] + (len(data),))

handle = open('results/denovo_blast_summary.txt', 'w')
handle.writelines('results\tdenovo_nbp\tblast_bits\tblast_E\tblast_ident_p\tblast_ident_n\tblast_gaps_p\tblast_gaps_n\tblast_taxon\n')


bits_re = re.compile(' [0-9]+ ?')
expect_re = re.compile(' [0-9e-]+ ?')
id_perc_re = re.compile('[0-9]+%')
id_n_re = re.compile('[0-9/]+')


# look at blocks
for denovo_ind, block in enumerate(blocks):
    
    # start looking at the start of the query
    block_data = data[block[0]:block[1]]
    seq_ln = block_data[2].replace('Length=','') # positional - bit fragile
    
    # extract matches - find the first sequence name
    names = [i for i, x in enumerate(block_data) if x.startswith('>')]
    if not names:
        handle.writelines(query_name[denovo_ind] + '\tNo hits found\t\t\t\t\t\t\n')
    else:
        match_data = block_data[names[0]:]
        
        # now find indices of bits of interest
        match_index = [i for i, x in enumerate(match_data) if x.startswith('>')]
        align_len_ind, align_len = zip(*[(i,x) for i, x in enumerate(match_data) if x.startswith('Length')])
        align_len = [x.replace('Length=','') for x in align_len]
        
        score = [x for x in match_data if x.startswith(' Score')]
        score = [x.split(',') for x in score]
        bits = [x[0] for x in score]
        bits = [bits_re.search(x).group() for x in bits]
        expect = [x[1] for x in score]
        expect = [expect_re.search(x).group() for x in expect]
        
        identity = [x for x in match_data if x.startswith(' Identities')]
        identity = [x.split(',') for x in identity]
        identities = [x[0] for x in identity]
        id_perc = [id_perc_re.search(x).group() for x in identities]
        id_n = [id_n_re.search(x).group() for x in identities]
        gaps = [x[1] for x in identity]
        gap_perc = [id_perc_re.search(x).group() for x in gaps]
        gap_n = [id_n_re.search(x).group() for x in gaps]
        
        # get the name
        names = [match_data[s:e] for s, e in zip(match_index, align_len_ind)]
        names = [' '.join(x) for x in names]
        
        # format the outputs
        n_matches = min(len(names), max_matches)
        
        for n in range(n_matches):
            handle.writelines('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(query_name[denovo_ind], 
                              seq_ln, bits[n], expect[n], id_perc[n], id_n[n], gap_perc[n], gap_n[n], names[n],))

handle.close()
