# All args come from sanpi command

Help() { 
echo -en "\033[31m"
cat << EOM

DESC:
Uses BWA aligner to align reads. 
Gives you sorted bam files.

++ Mandatory args: 
-t|--title [project_Xmen]
-c|--comments [crisprcount]
-f|--inputFile [input file. see below]
-F|--farm [no input]
-m|--mem [20000]   # becomes bigger with bam file size
-o|--organism ~/GENOMES/bwa/{genome_name}_bwa


++ inputFile format
file1.fq.gz
file2.fq.gz

++ Packages used
samtools, bwa

EOM
echo -en "\033[30m"
}


if [ $args_return  -gt 0 ]; then
	Help
	exit 2
fi

mandatory_fails=`mandatory title comments inputFile farm mem organism`
if [ `echo "$mandatory_fails" | wc -w` -gt 0 ]; then
	warnsms "mandatory check failed"
	errorsms "Set: $mandatory_fails"
fi

commandFile="$title.$comments.cmds"
echo -n "" > $commandFile
while read line
do
	l=($line)
	sampleFile=${l[0]}
	echo "export sampleFile=${sampleFile} ; export title=${title} ; export cmd0=$organism ; $mconf_Rscript $mconf_installdir/bin/crispr/crisprCount.R" >> $commandFile
done < $inputFile

# farm submission: commandfile var = commandFile
if [ ! -z "$farm" ]; then
 	export fileOfCommands=$commandFile
 	$mconf_installdir/bin/farmsub.sh
else 
 	warnsms "Not running as farm: Not recommended"
 	warnsms "$mconf_bashexec $commandFile"
fi
