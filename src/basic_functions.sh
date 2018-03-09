
## Functions for messaging to user
sms(){
	messaging=$1
	d=`date`
	echo "[ $d ] $messaging"
}

infosms(){
	messaging=$1
	echo "[ Info ] $messaging"
}

warnsms(){
	messaging=$1
	echo "[ Warning ] $messaging" > /dev/stderr
}
errorsms(){
	messaging=$1
	echo "[ Error ] $messaging" > /dev/stderr
	exit 1
}
print_stderr(){
	RED='\033[0;31m'
	NC='\033[0m' # No Color
	messaging=$1
	echo -e "${RED} $messaging ${NC}" > /dev/stderr
	exit 1
}

export -f infosms
export -f sms
export -f warnsms
export -f errorsms
