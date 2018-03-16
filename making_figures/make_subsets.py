import sys, gzip
import numpy as np

def main():
    prefix = sys.argv[1].split('.')[0]
    pairs_in = gzip.open(sys.argv[1])
    lines = pairs_in.readlines()
    length = len(lines)
    num_reads = np.arange(10000,length,10000) #number of reads to put in the output file
    header_len = int(sys.argv[2]) #how many lines to skip in pairs header (varies by species)

    np.random.seed(67)

    for i in num_reads:
        outfile_name = "%s_%dreads.pairs" % (prefix,i)
        pairs_out = open(outfile_name, 'w')
        selected_lines = np.random.choice(np.arange(header_len, length), i, replace=False)
        for rand in selected_lines:
            pairs_out.write(lines[rand].decode('utf-8'))
        pairs_out.close()

    pairs_in.close()

if __name__ == '__main__':
    main()