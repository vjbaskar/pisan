d=`date +"%d/%m/%y %T"`
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
	echo "[ Info ] $messaging" > /dev/stderr
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


export -f infosms
export -f sms
export -f warnsms
export -f errorsms
export -f sms_nonewline
