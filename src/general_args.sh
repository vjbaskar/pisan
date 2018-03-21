#!/usr/bin/env bash
# s = sampleFile
# i = id
# e = expType
# d = dataType
# c = comments
# p = procs
# m = mem
# t = title
# q = queue
# f = conf
# o = organism
# r = paired, 
# v = qval
# u = inputFile
# y = chrType


source ~/.$USER/pisan.conf
args_return=0
print_help() {
	cat $mconf_installdir/cfgs/args
	args_return=2
}

if [ $# -le 1 ]; then
	print_help
fi


#set -o errexit -o noclobber -o nounset -o pipefail
set -o errexit -o pipefail

params="$(getopt -o s:i:e:d:c:p:m:f:t:q:o:r:v:u:y:g:0:1:2:3:4:5:h -l sampleFile:,id:,expType:,dataType:,comments:,procs:,mem:,conf:,title:,queue:,organism:,paired:,qval:,inputFile:,chrType:,gtf:,file0:,file1:,file2:,file3:,file4:,file5:,help --name "$0" -- "$@")"
eval set -- "$params"

# --sampleFile samplefile.csv --id 3900STDY1234 --expType rnaseq --dataType genecounts --comments "Run by Vijay with DESeq2"

while true
do
	case "$1" in
		-t|--title)
			title=$2
			shift 2
		;;
		-s|--sampleFile)
			sampleFile=$2
			shift 2
		;;
		-i|--id)
			id=$2
			shift 2
		;;
		-e|--expType)
			expType=$2
			shift 2
		;;
		-d|--dataType)
			dataType=$2
			shift 2
		;;
		-c|--comments)
			comments=$2
			shift 2
		;;
		-p|--procs)
			procs=$2
			shift 2
		;;
		-m|--mem)
			mem=$2
			shift 2
		;;
		-q|--queue)
			queue=$2
			shift 2
		;;
		-f|--conf)
			conf=$2
			shift 2
		;;
		-o|--organism)
			organism=$2
			shift 2
		;;
		-r|--paired)
			paired=$2
			shift 2
		;;
		-v|--qval)
			qval=$2
			shift 2
		;;
		-u|--inputFile)
			inputFile=$2
			shift 2
		;;
		-y|--chrType)
			chrType=$2
			shift 2
		;;
		-g|--gtf)
			gtf=$2
			shift 2
		;;
		-0|--file0)
			file0=$2
			shift 2
		;;
		-1|--file1)
			file1=$2
			shift 2
		;;
		-2|--file2)
			file2=$2
			shift 2
		;;
		-3|--file3)
			file3=$2
			shift 2
		;;
		-4|--file4)
			file4=$2
			shift 2
		;;
		-5|--file5)
			file5=$2
			shift 2
		;;
		--)
			shift
			break
		;;
		-h|--help)
			print_help; exit 2
		;;
		?)
			print_help; exit 2
		;;
		*)
			echo "Not implemented: $1" >&2
			args_return=1
		;;
	esac
done


if [ -z "$queue" ]; then
	queue="normal"
fi
if [  -z "$mem" ]; then
	mem=1000
fi
if [  -z "$procs" ]; then
	procs=1
fi

