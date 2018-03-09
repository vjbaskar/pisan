#!/bin/sh


if [ $args_return -gt 0 ]; then
print_stderr "
	++ Mandatory args: procs mem paired organism inputFile 
	++ recommended: mem=32000   procs=1
	
	organism options: Your genome should be present in /nfs/users/nfs_v/vm11/GENOMES/STAR_INDEX/
"
	exit 2
fi

#threads=$1; shift
#paired=$1; shift
#genome=$1; shift
#namefile=$1; shift
#index=/nfs/users/nfs_v/vm11/GENOMES/STAR_INDEX/${genome}
#mem=32000


threads=$procs
genome=$organism
namefile=$inputFile
index=/nfs/users/nfs_v/vm11/GENOMES/STAR_INDEX/${organism}



seAlign() {
	id=$1
	name=$2
	STAR_INPUT="$name.fq.gz"
	OUTPUT_PREFIX="$name.star."
	echo "
		STAR --runThreadN  ${threads} \
		--genomeDir ${index} \
		--readFilesIn ${STAR_INPUT} \
		--outFileNamePrefix ${OUTPUT_PREFIX} \
		--readFilesCommand zcat \
		--outSAMtype BAM Unsorted \
		--outSAMunmapped Within \
		--quantMode GeneCounts
		bamsormadup SO=coordinate M=temp.metrics indexfilename=${OUTPUT_PREFIX}sormadup.bam.bai < ${OUTPUT_PREFIX}Aligned.out.bam > ${OUTPUT_PREFIX}sormadup.bam
	"
}


peAlign() {
	id=$1
	name=$2
	STAR_INPUT="$name.pair1.fq.gz $name.pair2.fq.gz"
	OUTPUT_PREFIX="$name.star."
	echo "
		STAR --runThreadN  ${threads} \
		--genomeDir ${index} \
		--readFilesIn ${STAR_INPUT} \
		--outFileNamePrefix ${OUTPUT_PREFIX} \
		--readFilesCommand zcat \
		--outSAMtype BAM Unsorted \
		--outSAMunmapped Within \
		--quantMode GeneCounts 
		bamsormadup SO=coordinate M=temp.metrics indexfilename=${OUTPUT_PREFIX}sormadup.bam.bai < ${OUTPUT_PREFIX}Aligned.out.bam > ${OUTPUT_PREFIX}sormadup.bam
		#samtools sort -@ $threads ${OUTPUT_PREFIX}Aligned.out.bam -o ${OUTPUT_PREFIX}bam
	"
}

echo -n "" > $namefile.star.commands
while read line
do
		l=($line)
		id=${l[0]}
		name=$id

		if [ $paired -eq 1 ]; then
			echo "Running in paired end mode"
			peAlign $id $name > $id.$name.star.sh
		else
			echo "Running in single end mode"
			seAlign $id $name > $id.$name.star.sh
		fi
		#sub_ncores normal star $mem $threads $name "sh $id.$name.star.sh"
		echo "sh $id.$name.star.sh" >> $namefile.star.commands

done < $namefile
#sub_ncores normal star $mem $threads $name "sh $id.$name.star.sh"
echo "farm..sub -q normal -n star -m $mem -p $threads -s batch -t 20 -a $namefile.star.commands"
