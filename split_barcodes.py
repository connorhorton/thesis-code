import os, sys, re
from operator import itemgetter
from itertools import izip, islice

def getExpectedBarcodes(ref_file, archive):
    barcodes = {}
    for line in ref_file:
        split = line.split('\t')
        if archive == split[0]:
            split = line.split('\t')
            barcodes[(split[2], split[3])] = split[4].strip('\n')
    if len(barcodes) == 0:
        raise ValueError("Invalid archive name. Must be in a form like Haploid_10")
    return barcodes

def statistics(statfile, outfiles, pfix, bc_dict, archive_name):
    print >> statfile, "%s\t%s\t%s\t%s\t%s" % ("Cell name", "5` barcode", "3` barcode", "# of reads", "Ploidy")
    stats = {}
    # count the number of reads in each file
    for filename in outfiles:
        # for now just count one out of the pair
        fq = open(pfix+filename+".1.fastq")
        num_reads = 0
        for line in fq:
            num_reads += 1
        for barcodes, name in bc_dict.iteritems():
            if name == filename:
                bc1, bc2 = barcodes
        if "Diploid" in archive_name:
            ploidy = "Diploid"
        else:
            ploidy = "Haploid"
        stats[(filename, bc1, bc2, ploidy)] = num_reads/4 # 4 lines per read
        fq.close()
    # sort files by number of reads
    sortedByReads = sorted(stats.items(), key=itemgetter(1), reverse=True)
    # print barcodes and number of reads in sorted order to stat file
    for item in sortedByReads:
        print >> statfile, "%s\t%s\t%s\t%d\t%s" % (item[0][0], item[0][1], item[0][2], item[1], item[0][3])
    statfile.close()

def demultiplex(fq1, fq2, bcfile1, bcfile2, outfiles, pfix, expBC, statfile, name):
    filename = ""
    forFiles = []
    revFiles = []

    fq1lines = list(islice(fq1, 4))
    fq2lines = list(islice(fq2, 4))
    bc1lines = list(islice(bcfile1, 4))
    bc2lines = list(islice(bcfile2, 4))

    while fq1lines and fq2lines and bc1lines and bc2lines:
        bc1 = bc1lines[1].strip('\n')
        bc2 = bc2lines[1].strip('\n')

        if (bc1, bc2) in expBC:
            filename = expBC[(bc1, bc2)]
            if filename in outfiles:
                idx = outfiles.index(filename)
                fout1 = forFiles[idx]
                fout2 = revFiles[idx]
            else:
                outfiles.append(filename)
                forFiles.append(open(pfix+filename+".1.fastq", 'w'))
                revFiles.append(open(pfix+filename+".2.fastq", 'w'))
                fout1 = forFiles[-1]
                fout2 = revFiles[-1]
                print "Making file for cell: %s" % (filename)
            fout1.write(''.join(fq1lines))
            fout2.write(''.join(fq2lines))
        #else do nothing
        fq1lines = list(islice(fq1, 4))
        fq2lines = list(islice(fq2, 4))
        bc1lines = list(islice(bcfile1, 4))
        bc2lines = list(islice(bcfile2, 4))
    else:
        # close files
        for file1, file2 in izip(forFiles, revFiles):
            file1.close()
            file2.close()
        # calculate stats
        statistics(statfile, outfiles, pfix, expBC, name)

def main():
    fq1 = open(sys.argv[1]) #fastq file 1 containing multiplexed sequence
    fq2 = open(sys.argv[2]) #fastq file 2 containing multiplexed sequence
    bcfile1 = open(sys.argv[3]) #fastq file 3 containing barcodes corresponding to fq1
    bcfile2 = open(sys.argv[4]) #fastq file 4 containing barcodes corresponding to fq2
    reference = open(sys.argv[5]) #text file containing expected barcodes
    stats = open(sys.argv[6], 'w') #output text file that will contain stats (num of reads)
    outfile_prefix = sys.argv[7] #gives the directory etc. where this is all located
    archive_name = sys.argv[8] #in the format "Haploid_10"
    outfiles = []

    exp_barcodes = getExpectedBarcodes(reference, archive_name)
    demultiplex(fq1, fq2, bcfile1, bcfile2, outfiles, outfile_prefix, exp_barcodes, stats, archive_name)

    fq1.close()
    fq2.close()
    bcfile1.close()
    bcfile2.close()
    reference.close()
    stats.close()

if __name__ == "__main__":
    main()