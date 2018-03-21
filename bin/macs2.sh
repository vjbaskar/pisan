
if [ $args_return  -gt 0 ]; then
print_stderr "
	++ Mandatory args: title, organism, paired, qval, inputFile, mem, procs
	
	organism options: hs,mm
	paired: 0 or 1
	qval: 0.05
	mem: 10000
	procs: 1
	***** inputFile format
	bamfile1.bam control1.bam
	bamfile2.bam control1.bam

"
	exit 2
fi

# Conda env installation of macs2
source activate p27

cmds="$title.macs.cmds"
rm -f $cmds
while read line
	do
	l=($line)
	bamfile=${l[0]}
	controlfile=${l[1]}
	if [ $paired -eq 1 ]; then
		fileFormat="BAMPE"
	else
		fileFormat="BAM"
	fi
	b=`basename $bamfile .bam` 
	b=`echo "${b}_macs_$qval"`
	c=`basename $controlfile .bam`
	odir="${b}_${c}.narrowpeak"
	macs_np="macs2 callpeak -t $bamfile -c $controlfile -f $fileFormat -g $organism -n $b -q $qval -B --SPMR --outdir $odir"
	odir="${b}_${c}.broadpeak"
	macs_bp="macs2 callpeak -t $bamfile -c $controlfile -f $fileFormat -g $organism -n $b -q $qval -B --SPMR --outdir $odir --broad -g hs --broad-cutoff 0.1"
	macs_vbp="macs2 callpeak -t $bamfile -c $controlfile -f $fileFormat -g $organism -n $b -q $qval -B --SPMR --outdir $odir --broad -g hs --broad-cutoff 0.2"
	echo $macs_np >> $cmds
	echo $macs_bp >> $cmds

done < $inputFile

farm..sub -q $queue -n macs2 -s $title -m $mem -p $procs -t 20 -a $cmds
