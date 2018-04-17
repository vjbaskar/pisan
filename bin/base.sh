# All args come from pisan command

if [ $args_return  -gt 0 ]; then
print_stderr "
	++ Mandatory args: 
	-t|--title
	-c|--comment
	-f|--inputFile
	-g|--genome = hs,mm
	-e|--expType
	
	++ recommended: 
	-F|--farm = flag
	-m|--mem=20000   # becomes bigger with bam file size
	-p|--procs=1
	
	organism options: hs,mm

	***** inputFile format
	bamfile1.bam control1.bam
	bamfile2.bam control1.bam

"
	exit 2
fi



# farm submission
# if [ ! -z "$farm" ]; then
# 	export fileOfCommands=$commandFile
# 	export concurrentJobs=20
# 	$mconf_installdir/bin/farmsub.sh
# else 
# 	warnsms "Not running as farm: Not recommended"
# 	warnsms "$mconf_bashexec $commandFile"
# fi
