# Script: Get data from irods

# Main config file
source ~/.$USER/pipeline/conf
#! $mconf_bashexec

# get the exec dir
basedir=`readlink -f $0 | xargs dirname`
source $mconf_installdir/src/general_args.sh
source $mconf_installdir/src/basic_functions.sh
echo $args_return
if [ $args_return  -gt 0 ]; then
	echo "
	++ Mandatory args: sampleFile, procs
	"
	exit 2
fi

source $conf # config file in the local directory can override the default ones
