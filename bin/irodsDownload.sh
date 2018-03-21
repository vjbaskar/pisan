# Get data from irods


if [ $args_return  -gt 0 ]; then
print_stderr "
++ Mandatory args: title, sampleFile
"
exit 2
fi


#Set Script Name variable
THREADS=$procs



function md5sumCheck(){
	irods_file=$1
	local_file=$2
	returnVal=0
	md5=`imeta ls -d $irods_file md5 | awk ' $1 ~ /value/ { print $2 } '`
	lmd5=`md5sum $local_file ! awk ' { print $1 } '`
	echo $md5
	echo $lmd5
	if [ $md5 == $lmd5 ]; then
		returnVal=1
	fi
	echo $returnVal
}

# Get data from irods
function getDataFromIrods(){

function md5sumCheck(){
	irods_file=$1
	local_file=$2
	returnVal=0
	md5=`imeta ls -d $irods_file md5 | awk ' $1 ~ /value/ { print $2 } '`
	lmd5=`md5sum $local_file | awk ' { print $1 } '`
	if [ $md5 == $lmd5 ]; then
		returnVal=1
	fi
	echo $returnVal
}

	filetype="cram"
	runId=$1; shift
	sample=$1; shift # MOLM_150_C
	lanes=$1; shift # 3,4,5,6
	lanes=`echo "$lanes" | sed -e "s/,/ /g"`

	echo "
	*******
	run id = $runId
	sample name = $sample
	file type = $filetype
	lanes = $lanes
	*******
	" 

	ofile=${sample}
	targetFolder=$sample
	
	rm -rf $targetFolder
	if [ ! -z $targetFolder ]; then
		infosms "Dumping data in $targetFolder"
		mkdir $targetFolder
	fi
#imeta qu -z seq -d id_run = $id_run and target = 1 and manual_qc = 1 and sample = $sample > $imetafile
	rm -f $imetafile

	mkdir ${targetFolder}/finalfiles/
	for i in $lanes
		do
			md5checkVal=0
			imetafile=$targetFolder/${runId}.${sample}.$i.imeta
			echo "Getting data for lane = $i"
			imeta qu -z seq -d type = $filetype and target = 1 and lane = $i and id_run = $runId and sample = $sample > $imetafile
			wc=`cat $imetafile | grep cram | wc -l`
			p=`pwd`
			if [ $wc -gt 0 ]; then
				seqfolder=(`awk ' $1 ~ "collection:" { print $2 } ' $imetafile`)
				files=(`awk ' $1 ~ "dataObj:" { print $2 } ' $imetafile`)
				fileName=${seqfolder}/${files}
				mkdir -p $targetFolder/${runId}
				iget $fileName $targetFolder/${runId}/${files}
				md5checkVal=$(md5sumCheck $fileName $targetFolder/${runId}/${files})
				if [ "$md5checkVal" -eq 1 ]; then
					mv $targetFolder/${runId}/${files} $targetFolder/finalfiles/
				else
					echo "md5checksum failed"
				fi
				echo "$runId	$sample	$i	$seqfolder/$files	mdcheck=$md5checkVal" >> samples.irods.log
					
			fi
		done
	return 0
}


#### codes start
#### Get data
pipe..createProject --sampleSheet $sampleFile -c RUNID,SAMPLEID,LANES | sed -e '1d' >| samples.irods.txt
rm -f samples.irods.log
export -f getDataFromIrods
#getDataFromIrods 24076 3960STDY7100225     1
rm -f samples.irods.log
cat samples.irods.txt | xargs -n3 -I@ -P $THREADS bash -c "$(declare -f getDataFromIrods); getDataFromIrods @"

