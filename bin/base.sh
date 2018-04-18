# All args come from pisan command

Help() { 
echo -en "\033[31m"
cat << EOM
++ Mandatory args: 
-t|--title
-c|--comments
-f|--inputFile

++ recommended: 
-F|--farm = flag
-m|--mem=20000   # becomes bigger with bam file size
-p|--procs=1

organism options: hs,mm

***** inputFile format
4010STDY7036819	C_3143NTA
4010STDY7036820	C_3143NTB
4010STDY7036821	C_3143NTC
4010STDY7036822	C_3143CISWTA
4010STDY7036823	C_3143CISWTB
4010STDY7036824	C_3143CISWTC
4010STDY7036825	C_3143CISMUTA
4010STDY7036826	C_3143CISMUTB
4010STDY7036827	C_3143CISMUTC
	
EOM
echo -en "\033[30m"
}


if [ $args_return  -gt 0 ]; then
	Help
	exit 2
fi

mandatory_fails=`mandatory title comments inputFile`
if [ `echo "$mandatory_fails" | wc -w` -gt 0 ]; then
	warnsms "mandatory check failed"
	errorsms "Set: $mandatory_fails"
fi



# farm submission: commandfile var = commandFile
# if [ ! -z "$farm" ]; then
# 	export fileOfCommands=$commandFile
# 	$mconf_installdir/bin/farmsub.sh
# else 
# 	warnsms "Not running as farm: Not recommended"
# 	warnsms "$mconf_bashexec $commandFile"
# fi
