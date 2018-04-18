#!/usr/bin/env bash
# c = comments
# d = dataType
# e = expType
# f = conf
# i = id
# m = mem
# o = organism
# p = procs
# q = queue
# r = paired, 
# s = sampleFile
# t = title
# u = inputFile
# v = qval
# x = command
# y = chrType
# a = array
# W = walltime
# T = concurrent
# j = waitForJob
# w = wait
# G  = farm group
# F = submit in farm or HPC



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

if [ $# -gt 4 ]; then
	d=`date`
	echo "[ $d ] $@" > .pisan.cmds
fi


SHORT="s:i:e:d:c:p:m:f:t:q:o:r:v:u:y:x:g:a:W:T:j:G:Fw0:1:2:3:4:5:h"
LONG="sampleFile:,id:,expType:,dataType:,comments:,procs:,mem:,inputFile:,title:,queue:,organism:,paired:,qval:,conf:,chrType:,command:,gtf:,array:,walltime:,concurrent:,waitForJob:,farmgroup:,farm,wait,file0:,file1:,file2:,file3:,file4:,file5:,help"
#params="$(getopt -o s:i:e:d:c:p:m:f:t:q:o:r:v:u:y:x:g:a:W:T:j:w1:2:3:4:5:h -l sampleFile:,id:,expType:,dataType:,comments:,procs:,mem:,conf:,title:,queue:,organism:,paired:,qval:,inputFile:,chrType:,command:,gtf:,array:,walltime:, concurrent:,waitForJob:,wait,file0:,file1:,file2:,file3:,file4:,file5:,help "$0" -- "$@")"
params="$(getopt -o $SHORT -l $LONG --name "$0" -- "$@")"
eval set -- "$params"

# --sampleFile samplefile.csv --id 3900STDY1234 --expType rnaseq --dataType genecounts --comments "Run by Vijay with DESeq2"


while true
do
	case "$1" in
		-t|--title) title=$2; shift 2;; # general
		
		-s|--sampleFile) sampleFile=$2; shift 2;; # general
		
		-i|--id) id=$2; shift 2;; # sanger id
		
		-e|--expType) expType=$2; shift 2;; # rna-seq, chip-seq
		
		-d|--dataType) dataType=$2; shift 2;; # genecount,exoncount,tpms,fpkms
		
		-c|--comments) comments=$2; shift 2;; # Use this for creating sample ids in farm submissions
		
		-p|--procs) procs=$2; shift 2;;
		
		-m|--mem) mem=$2; shift 2;;
			
		-q|--queue)	queue=$2; shift 2;;

		-f|--inputFile) inputFile=$2; shift 2;;
			
		-o|--organism) organism=$2; shift 2;; 
			
		-r|--paired) paired=$2; shift 2;;
		
		-v|--qval) qval=$2; shift 2;;
		
		-u|--conf) conf=$2; shift 2;; # command specific
			
		-y|--chrType) chrType=$2; shift 2;; # ensembl or ucsc
			
		-x|--command) command=$2; shift 2;; # farm sub option
			
		-g|--gtf) gtf=$2; shift 2;;
		
		-a|--array) fileOfCommands=$2; shift 2;; # farm sub option
		
		-W|--walltime) walltime=$2; shift 2;; # farm sub option
			
		-T|--concurrent) concurrent=$2; shift 2;; # farm sub option
		
		-j|--waitForJob) jobIDtoWait=$2; shift 2;; # farm sub option
		
		-G|--farmgroup) farmgroup=$2; shift 2;; # farm sub option
		
		-w|--wait) wait_to_finish=1; shift;; # farm sub option
		
		-F|--farm) farm=1; shift;; # farm sub option
		
		-0|--file0) file0=$2; shift 2;;
		
		-1|--file1) file1=$2; shift 2;;
			
		-2|--file2) file2=$2; shift 2;;

		-3|--file3) file3=$2; shift 2;;
			
		-4|--file4) file4=$2; shift 2;;
			
		-5|--file5) file5=$2; shift 2;;
			
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
if [  -z "$farmgroup" ]; then
	farmgroup="team163-grp"
fi

