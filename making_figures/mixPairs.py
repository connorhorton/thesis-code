import os, sys, gzip
import numpy as np
from itertools import combinations

def main():
    keys = []
    pairs_in = {}
    pairs_out = {}
    lines = {}
    lengths = {}
    for i in range(1, len(sys.argv)-2):
        keys.append(sys.argv[i])
        pairs_in[sys.argv[i]] = gzip.open(sys.argv[i]) #pairs file containing scHiC contacts from one cell
        lines[sys.argv[i]] = pairs_in[sys.argv[i]].readlines()
        lengths[sys.argv[i]] = len(lines[sys.argv[i]])
    num_reads = int(sys.argv[len(sys.argv)-2]) #number of reads to put in the output file
    header_len = int(sys.argv[len(sys.argv)-1]) #how many lines to skip in pairs header (varies by species)

    np.random.seed(67)
    total_length = 0

    for i in range(1, len(keys)+1):
        count = 0
        for subset in combinations(keys, i): # each combination of files
            outfile_name = "%dchoose%d_%0.3d.pairs" % (len(keys),i,count)
            print("preparing %s" % outfile_name)
            pairs_out = open(outfile_name, 'w')
            for file in subset:
                total_length += lengths[file]
            try:
                selected_lines = np.random.choice(total_length, num_reads, replace=False) # sample without replacement specified number of lines
                for rand in selected_lines: # loop through random numbers
                    k = 0
                    while rand >= lengths[subset[k]]:
                        rand -= lengths[subset[k]]
                        k += 1
                    if rand > header_len:
                        pairs_out.write(lines[subset[k]][rand].decode('utf-8'))
                pairs_out.close()             
            except ValueError:
                print("subset %s does not have 50,000 reads. It only has %d reads." % (subset, total_length))
                os.unlink(outfile_name)        

            # prepare for the next iteration
            total_length = 0
            count += 1

    for f in pairs_in:
        pairs_in[f].close()


if __name__ == "__main__":
    main()