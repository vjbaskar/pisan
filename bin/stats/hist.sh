# All args come from legopipe command

export nohistory=1 ; 
Help() { 
echo -en "\033[31m"
cat << EOM
++ Mandatory args: 
-P|--progArgs: cat, more, less, edit, grep


++ optional: 
-c|--comments: simple search. For eg. use date, command names etc. Uses linux grep.
-i|--id: local

++ example:

legopipe hist -P cat
legopipe hist -P cat -c "Wed Apr 2"

EOM
echo -en "\033[30m"
}

if [ "$id" == "local" ]; then
	legopipe_cmds="./.legopipe/hist.cmds"
else
	legopipe_cmds="$HOME/.legopipe/hist.cmds"
fi



if [ $args_return  -gt 0 ]; then
	Help
	exit 2
fi

mandatory_fails=`mandatory progArgs`
if [ `echo "$mandatory_fails" | wc -w` -gt 0 ]; then
	warnsms "mandatory check failed"
	errorsms "Set: $mandatory_fails"
fi

if [ ! -z "$comments" ]; then
	grep "$comments" $legopipe_cmds
else 
	if [ "$progArgs" == "edit" ]; then
		progArgs="vim"
	fi
	if [ "$progArgs" == "grep" ]; then
		progArgs="grep -i "
	fi
	$progArgs $legopipe_cmds
fi
