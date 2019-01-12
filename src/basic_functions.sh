d=`date +"%d/%m/%y %T"`
legopipeID=`date +%Y-%m-%d-%H-%M-%S`
legopipeID="legopipe.${legopipeID}.${RANDOM}"
## Functions for messaging to user
sms(){
	messaging=$1
	echo "[ $d ] $messaging" > /dev/stderr
}

sms_nonewline(){
	messaging=$1
	echo -n "[ $d ] $messaging"
}

infosms(){
	messaging=$1
	echo -e "\033[0;33m[ $d ] $messaging\033[0m" > /dev/stderr
}

warnsms(){
	messaging=$1
	echo -e "\033[0;33m[ Warning ] $messaging\033[0m" > /dev/stderr
}
errorsms(){
	messaging=$1
	echo -e "\033[0;31m[ Error ] $messaging\033[0m" > /dev/stderr
	exit 1
}
print_stderr(){
	RED='\033[0;31m'
	NC='\033[0m' # No Color
	messaging=$1
	echo -e "${RED} $messaging ${NC}" > /dev/stderr
	exit 1
}

mandatory(){
	mandatory_check=0
	while [ $# -gt 0 ]; 
	do
		x=$1
		if [ -z ${!x} ]; then
			printf "$x "			
		fi
		shift
	done
}

commandsHelp(){
    mainCommand=$commandName
    grep "#@" $mainCommand | awk -F "#@" ' { command=gensub("\"|)","","g", $1) }  { print command, $2 } '   >&2
    exit 0
}

copy_command(){
	cmd=$1
	d=`date +"%d_%m_%y_%T"`
	c=`echo $cmd | awk ' { print $NF } '`
	if [ ! -z $cwd/.legopipe ]; then
		mkdir -p $cwd/.legopipe	
	fi
	temp=`basename $c`
	cp $c $cwd/.legopipe/$legopipeID.$temp

}

toLower() {
	echo "$1" | tr '[A-Z]' '[a-z]'
}

toUpper() {
	echo "$1" | tr '[a-z]' '[A-Z]' 
}


export -f infosms
export -f sms
export -f warnsms
export -f errorsms
export -f sms_nonewline
export legopipeID=$legopipeID
