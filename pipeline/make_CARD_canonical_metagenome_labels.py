import sys

if __name__ == '__main__':
    
    labels_fh = open(sys.argv[2], 'w')
    with open(sys.argv[1]) as fh:
        for line in fh:
            line = line.strip().replace('@', '')
            read_name = line
            contig_name = '-'.join(line.split('-')[:-1])

            split_line = line.split('|')
            aro = split_line[4].replace('ARO:', '')
            amr_name = '-'.join(split_line[5].split('-')[:-1])
            
            # as they are all straight from the canonical sequence
            # and can't partially overlap 
            amr_cutoff = 'Perfect'
            amr_overlap = '250'
            
            labels_fh.write('\t'.join([read_name, contig_name, aro, amr_name,
                                       amr_cutoff, amr_overlap]) + '\n')
    labels_fh.close()
    print('Labels generated')
