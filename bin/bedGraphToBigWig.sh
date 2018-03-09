#!/bin/sh

if [ $args_return -gt 0 ]; then
	print_stderr "
	BigWig file for chipseq from bedgraph file
	
	++ Mandatory args: procs mem organism inputFile
	++ recommended: mem=20000   procs=1

	==== bedGraphFile.txt
	10_MUT_NPM1_macs_0.05.spmr_treat_pileup.bdg
	10_TY1_macs_0.05.spmr_treat_pileup.bdg
	10_WT_NPM1_macs_0.05.spmr_treat_pileup.bdg
	30_MUT_NPM1_macs_0.05.spmr_treat_pileup.bdg
	30_TY1_macs_0.05.spmr_treat_pileup.bdg
	30_WT_NPM1_macs_0.05.spmr_treat_pileup.bdg


	**** You should have chr.sizes as ~/GENOMES/genomeSizes/genome.chr.sizes
	"
	exit 2
fi


genome=$organism
gtype=$chrType
bedGraphList=$inputFile
genomeSizes="~/GENOMES/genomeSizes/${genome}.chr.sizes"
mem=$mem
nProcs=$procs
#genome=$1; shift
#gtype=$1; shift
#bedGraphList=$1;shift
#genomeSizes="~/GENOMES/genomeSizes/${genome}.chr.sizes"
#mem=20000
#nProcs=1
#echo "bedGraphList = $bedGraphList"


while read line
do
	l=($line)
	file=${l[0]}
	echo "File = $file"
	f=`basename $file .bdg`
	tfile=`tempfile -d ./`
	command=""
	if [ "$gtype" == "ensembl" ]; then
		echo "awk ' OFS=\"\t\" { if(\$1 ~/^[0-9]|^X|^Y/) print \"chr\"\$0 }  ' $file | sortBed -i - > $tfile" > $f.bdg2bw.sh
	else 
		echo "sortBed -i $file > $tfile" > $f.bdg2bw.sh
	fi
	echo "bedGraphToBigWig $tfile $genomeSizes $f.bw" >> $f.bdg2bw.sh
	echo "rm -f $tfile" >> $f.bdg2bw.sh
	farm..sub -q $queue -n bw -m $mem -p $procs -s $f -c "sh $f.bdg2bw.sh"
done < $bedGraphList

