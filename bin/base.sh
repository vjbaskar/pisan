# Script: MACS2 data

# Main config file
#source ~/.$USER/pisan.conf
#! $mconf_bashexec

# get the exec dir
# basedir=`readlink -f $0 | xargs dirname`
#source $mconf_installdir/src/general_args.sh
#source $mconf_installdir/src/basic_functions.sh

#if [ ! -z "$conf" ]; then
#	source $conf # config file in the local directory can override the default ones
#fi


if [ $args_return  -gt 0 ]; then
print_stderr "
	++ Mandatory args: organism, paired, qval, inputFile, mem, procs
	++ recommended: mem=10000   procs=1
	
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