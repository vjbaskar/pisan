# All args come from pisan command

Help() { 
echo -en "\033[31m"
cat << EOM
++ Mandatory args: 
-P|--progArgs: cat, more, less, edit

++ optional: 
-c|--comments: simple search. For eg. use date, command names etc. Uses linux grep.

++ example:

pisan hist -P cat
pisan hist -P cat -c "Wed Apr 2"

EOM
echo -en "\033[30m"
}

pisan_cmds="$HOME/.${USER}/pisan/pisan.cmds"

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
	grep "$comments" $pisan_cmds
else 
	if [ "$progArgs" == "edit" ]; then
		progArgs="vim"
	fi
	$progArgs $pisan_cmds
fi