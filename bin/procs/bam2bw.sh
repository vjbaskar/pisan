if [ $args_return  -gt 0 ]; then
print_stderr "
	++ Mandatory args: 
	-t|--title
	-c|--comment
	-f|--inputFile
	-g|--genome = hs,mm
	-e|--expType = rnaseq, chipseq
	
	++ recommended: 
	-F|--farm = flag
	-m|--mem=20000   # becomes bigger with bam file size
	-p|--procs=1
	
	
	***** inputFile format
	bamfile1.bam 
	bamfile2.bam 

"
	exit 2
fi


if [ "$genome" == "hs" ]; then
    gsize=3200000000
fi
if [ "$genome" == "mm" ]; then
    gsize=2159570000
fi

commandFile="$title.$comments.cmds"
echo -n "" > $commandFile
while read line
do
	l=($line)
	bam=${l[0]}
	b=`basename $bam .bam`
	cmd=""
	if [ ! -e $bam.bai ]; then
			cmd="samtools index $bam ;"
	fi
	if [ "$expType" == "rnaseq" ]; then
		cmd="$cmd  $mconf_deeptools/bamCoverage -b $bam -o $b.fw.bw  --minMappingQuality 1 -bs 2 --normalizeUsingRPKM -p $procs --filterRNAstrand forward" 
		echo "$cmd" >> $commandFile
		cmd=""
		cmd="$cmd  $mconf_deeptools/bamCoverage -b $bam -o $b.rev.bw  --minMappingQuality 1 -bs 2 --normalizeUsingRPKM -p $procs --filterRNAstrand reverse" 
        echo "$cmd" >> $commandFile
    else
        cmd="$cmd $mconf_deeptools/bamCoverage -b $bam -o $b.bw  --minMappingQuality 1 -bs 2 --normalizeUsingRPKM -p $procs" 
        echo "$cmd" >> $commandFile
    fi
done < $inputFile

#farm submission
if [ ! -z "$farm" ]; then
	export fileOfCommands=$commandFile
	export concurrentJobs=20
	$mconf_installdir/bin/farmsub.sh
else 
	warnsms "Not running as farm: Not recommended"
	warnsms "$mconf_bashexec $commandFile"
fi
