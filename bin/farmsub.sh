if [ $args_return  -gt 0 ]; then
print_stderr "
	++ Mandatory args: mem,procs,queue,command,title,comments

"
	exit 2
fi



## FUnctions for setting the variables if they are optional
optional(){
	if [  -z "$queue" ]; then
		warnsms "Queue is undefined. Setting it to normal"
		export queue="normal"
	fi
	if [  -z "$mem" ]; then
		warnsms "Memory not set. Setting it to default = 1000 mb"
		export mem=1000
	fi
	if [  -z "$procs" ]; then
		warnsms "procs not set. Setting it to default = 1"
		export procs=1
	fi
}

mandatory(){
	echo "Mandatory Args"
}

## Functions to generate bsub script

farmGroup(){
	c=`echo $HOSTNAME | grep cgp -c`
	if [ "$c" -eq "0" ]; then
		echo "#BSUB -G ${farmgroup}"
	fi
}

coreBSUB() {
# - core bsub function
	#command=$1
	d=`date +%_F.%Hhr`
	ofile="${title}.${comments}"
echo "
#!/bin/sh
#source ~/.bashrc
#BSUB -M ${mem} 
#BSUB -R \"select[mem>${mem}] rusage[mem=${mem}] span[hosts=1] \" 
#BSUB -n $procs 
#BSUB -q ${queue}"
farmGroup
}

jobName(){
# - Job name
	echo "#BSUB -J ${ofile}"
}

bsubFiles(){
# - Out and err files
	echo "
#BSUB -o .bsub/${ofile}.farm
#BSUB -e .bsub/${ofile}.farm"
}
waitBSUB(){
# - Wait till finish
	echo "#BSUB -K"
}

waitForJob(){
# - Wait till finish of a specific job
echo "#BSUB -w 'done(${jobIDtoWait})'"
}

arrayBSUB(){
# - Run as array
	totalJobs=`cat $fileOfCommands | wc -l`
	if [ -z ${concurrentJobs} ]; then
		concurrentJobs=${totalJobs}
	fi
echo "
#BSUB -o .bsub/${ofile}.%I.farm
#BSUB -e .bsub/${ofile}.%I.farm
#BSUB -J \"${ofile}[1-$totalJobs]%${concurrentJobs}\"
command=\`sed -e \"\${LSB_JOBINDEX}q;d\" $fileOfCommands\`
echo \$command
\$command
"
}

timeBSUB(){
	t=$1
echo "#BSUB -W $t"
}


commandBSUB(){
echo "
echo \"Job id = \$LSB_JOBID \"
d=\`date\`
echo \"[ \$d ] ${ofile} starts\"
$command
d=\`date\`
echo \"[ \$d ] ${ofile} ends\"
" 
}

writeBSUB() {
	# -- A function to create submission scripts for LSF

	if [ ! -d .bsub ]; then
			mkdir .bsub
	fi

	coreBSUB
	if [ ! -z "$wait_to_finish" ]; then
		waitBSUB
	fi
	if [ ! -z "$jobIDtoWait" ]; then
		waitForJob
	fi
	if [ ! -z "$walltime" ]; then
		timeBSUB $walltime
	fi
	if [ ! -z "$fileOfCommands" ]; then
		command=`arrayBSUB`
	else
		bsubFiles
		jobName
	fi
commandBSUB

}



#### End of functions

# If the values are undefined, default them
optional

# Generate the bsub script
writeBSUB "$command" > ${title}.${comments}.bsub

# Submit job
temp=`tempfile`
bsub < ${title}.${comments}.bsub > $temp


# Get job id
jid=`head -1 $temp | awk ' { print $2 } ' | sed -e "s/>//g" | sed -e "s/<//g"`
print_stderr ">> Job title: ${title}.${comments}.bsub >> Job ID = $jid"
echo "[date = `date`] [operation = $title] [comments = $comments] [bsub_file = $ofile] [job_id = $jid] [log = ${ofile}] [command = ${fullcommand}]" >> joblists.lsf
print_stderr "" >> joblists.lsf
echo "$jid"
