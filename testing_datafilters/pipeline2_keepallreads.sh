#!/bin/bash

:<<'end_long_comment'
---- COPYRIGHT ----------------------------------------------------------------
Copyright (C) 2017-2018
Connor Horton (Harvard University)

---- LICENSE ------------------------------------------------------------------
This file is part of HiClass.

HiClass is free software: you can redistribute it and/or modify it under the
terms of the GNU Lesser General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option) any
later version.

HiClass is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
details. 

You should have received a copy of the GNU Lesser General Public License along
with this software.  If not, see <http://www.gnu.org/licenses/>.
end_long_comment
#add pairix to path
export PATH=~/pairix/bin/:$PATH
#parse command line arguments
usage() {
	cat <<EOM
Usage: bash HiClass.sh [-h] [-d] -i DIRECTORY -fq FASTQ1 FASTQ2 -e RE1 [RE2 RE3] -s FRAG_SIZE -g GENOME -b BWA_PATH -r GENOME_INDEX -q PAIRSQC_PATH
	-h				prints this message
	-d 				indicates a diploid cell
	-i input			prefix for fastq files (e.g. if fastq files are cell.1.fastq and cell.2.fastq, then the prefix will be cell)
	-o output			desired name for output files, including pairs file
	-e restriction_enzyme		use -e for each restriction enzyme (e.g. -e MboI -e AluI). See enzymes.conf for allowed enzymes
	-s fragment_size		allowed range of sequenced molecule sizes in bp (e.g. 150-3000)
	-g genome			genome build (either mm10 or hg38)
	-b bwa_path			directory where Burrows-Wheeler Aligner is installed
	-r ref_genomes_path		directory where reference genomes lie
	-q pairsQC_path			directory where PairsQC is installed
EOM
	exit 1;
}
if [ $# == 0 ] ; then
    usage
    exit 1;
fi
#default options
DIPLOID=false
DIR=$(pwd)
NUM_CPU=1
RE1="MboI"
GENOME="mm10"
POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"
case $key in
    -h|--help)
	    usage
	    shift # past argument
	    shift # past value
	    ;;
    -i|--directory)
	    DIR="$2"
	    shift # past argument
	    shift # past value
	    ;;
    -f)
	    FASTQ1="$2"
	    FASTQ2="$3"
	    shift # past argument
	    shift # past value
	    ;;
    -e)
		RE1="$2"
        if [ -n "$3" ]; then
            RE2="$3"
        fi
        if [ -n "$4" ]; then
        	RE3="$4"
        fi
        shift
        shift
        ;;
	-a|--assembly)
		if [ "$2" == "hg38" ] || [ "$2" == "mm10" ]; then
	    	GENOME="$2"
	    else
	    	echo "Must specify mm10 or hg38 genome assembly"
	    	usage
	    fi
	    shift # past argument
	    shift # past value
    	;;
	-g|--genomeindex)
		GENOME_INDEX="$2"
	    shift # past argument
	    shift # past value
    	;;
    -b|--bwa)
	    BWA_PATH="$2"
	    shift # past argument
	    shift # past value
    	;;
    -q|--pairsqc)
	    PAIRSQC_PATH="$2"
	    shift # past argument
	    shift # past value
	    ;;
    -n)
	    NUM_CPU="$2"
	    shift # past argument
	    shift # past value
	    ;;
    -d)
	    DIPLOID=true
	    shift # past argument
	    shift
	    ;;
    *)    # unknown option
    	POSITIONAL+=("$1") # save it in an array for later
    	shift # past argument
    	;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters
#check to make sure parameters are here
if [ -z ${FASTQ1+x} ] && [ -z ${FASTQ2+x} ] ; then echo "Error: Must specify two FASTQ files" ; echo ; usage ; fi
if [ -z ${BWA_PATH+x} ] ; then echo "Error: Must specify directory where BWA is located" ; echo ; usage ; fi
if [ -z ${PAIRSQC_PATH+x} ] ; then echo "Error: Must specify directory where PairsQC is located" ; echo ; usage ; fi
if [ -z ${GENOME_INDEX+x} ] ; then echo "Error: Must specify where genome index is located" ; echo ; usage ; fi
FASTQ_ROOT=$(echo $FASTQ1 | cut -f 1 -d '.')
WD=$(pwd)
# clip fastq -- make sure to do file prefix
# make sure to also get re stuff
# maybe make processing folder for all these rando files we generate?

echo "clipping reads and discarding short reads..."
python clip_reads.py $FASTQ1 $DIR/$FASTQ_ROOT $RE1 False
python clip_reads.py $FASTQ2 $DIR/$FASTQ_ROOT $RE1 True
#edit fastq for bwa-mem
echo "editing fastq files in preparation for alignment..."
python edit_fastq.py $DIR/${FASTQ_ROOT}_reads1_clipped.fastq $DIR/${FASTQ_ROOT}_reads2_clipped.fastq $DIR/$FASTQ_ROOT
#map (bwa-mem)
echo "aligning reads..."
cd $BWA_PATH
./bwa mem -t $NUM_CPU -SP5M $GENOME_INDEX $DIR/${FASTQ_ROOT}_reads1_edited.fastq $DIR/${FASTQ_ROOT}_reads2_edited.fastq > $DIR/${FASTQ_ROOT}_raw.sam 2> $DIR/align.err
# #pairsamtools
cd $WD
if [ ! -f helper_files/$GENOME.$RE1.fragments.bed ]; then
    echo "digesting genome with $RE1"
    echo cooler digest helper_files/$GENOME.chrom.sizes.mainonly $(dirname "${GENOME_INDEX}") $RE1 -o helper_files/$GENOME.$RE1.fragments.bed
fi
echo "pairsamtools!" 
pairsamtools parse -c helper_files/$GENOME.chrom.sizes.mainonly --assembly $GENOME --output-stats $DIR/${FASTQ_ROOT}_stats.txt $DIR/${FASTQ_ROOT}_raw.sam | {
pairsamtools sort --nproc 1
} | {
pairsamtools select '(pair_type=="UU") or (pair_type=="UR") or (pair_type=="RU")' --output-rest >( pairsamtools split \
            --output-pairs $DIR/${FASTQ_ROOT}.unmapped.pairs \
            --output-sam $DIR/${FASTQ_ROOT}.unmapped.sam ) 
} | {
pairsamtools restrict -f helper_files/$GENOME.$RE1.fragments.bed
} | {
pairsamtools split \
            --output-pairs $DIR/${FASTQ_ROOT}.pairs \
            --output-sam $DIR/${FASTQ_ROOT}.sam
}
#filter (each separate script?)
# generate ncc for modeling
# python pairsam2ncc.py $DIR/${FASTQ_ROOT}.pairs $DIR/${FASTQ_ROOT}.sam $DIR/${FASTQ_ROOT}.ncc
#remove promiscuous
# if [ $DIPLOID ] ; then
# 	python remove_promiscuous.py $DIR/$FASTQ_ROOT
# fi
echo "zipping and indexing pairs file"
bgzip -f $DIR/$FASTQ_ROOT.pairs
pairix -f $DIR/$FASTQ_ROOT.pairs.gz
#run pairsqc!
REPORT_DIR=$(echo $DIR | rev | cut -f 1 -d '/' | rev) # get report directory for PairsQC
echo "running PairsQC"
if [[ $GENOME == "mm10" ]]; then
		python $PAIRSQC_PATH/pairsqc.py -p $DIR/$FASTQ_ROOT.pairs.gz -c helper_files/$GENOME.chrom.sizes.mainonly -t P -s $REPORT_DIR/$FASTQ_ROOT -M 8.2
elif [[ $GENOME == "hg38" ]]; then
		python $PAIRSQC_PATH/pairsqc.py -p $DIR/$FASTQ_ROOT.pairs.gz -c helper_files/$GENOME.chrom.sizes.mainonly -t P -s $REPORT_DIR/$FASTQ_ROOT
fi
echo "calculating gini index"
cooler cload pairix helper_files/$GENOME.chrom.sizes.mainonly:1000000 $DIR/$FASTQ_ROOT.pairs.gz $DIR/$FASTQ_ROOT.gini.cool 2>> $DIR/${FASTQ_ROOT}_std.err
python gini.py $DIR/$FASTQ_ROOT.gini.cool report/$REPORT_DIR/$FASTQ_ROOT.cis_to_trans.out
echo "generating cool files"
cooler cload pairix helper_files/$GENOME.chrom.sizes.mainonly:10000 $DIR/$FASTQ_ROOT.pairs.gz $DIR/$FASTQ_ROOT.cool 2>> $DIR/${FASTQ_ROOT}_std.err
cooler zoomify $DIR/$FASTQ_ROOT.cool -r 10000,50000,100000,250000,500000,1000000,2500000,5000000,10000000 2>> $DIR/${FASTQ_ROOT}_std.err
#remove unnecessary files
rm $DIR/${FASTQ_ROOT}_reads1_clipped.fastq $DIR/${FASTQ_ROOT}_reads2_clipped.fastq
rm $DIR/${FASTQ_ROOT}_reads1_edited.fastq $DIR/${FASTQ_ROOT}_reads2_edited.fastq
rm $DIR/${FASTQ_ROOT}.unmapped.pairs $DIR/${FASTQ_ROOT}.unmapped.sam
rm $DIR/${FASTQ_ROOT}_raw.sam