
if [ $args_return  -gt 0 ]; then
print_stderr "
	++ Mandatory args: title, organism
	
	options: 
	-t|--title = short name of the hub
	-c|--comments = long name of the hub
	-o|--organism = ucsc organism names such as hg19,hg38,mm10 etc
	-f|--inputFile
	
	inputFile format:
	-----------------
	file1.bw	10,12,255
	file2.bw	255,40,1
	
	Get rgb colours from https://www.rapidtables.com/web/color/RGB_Color.html
"
	exit 2
fi

p=`pwd`
BASE_DIR=$mconf_installdir
name=$title
genome=$organism
inputFile=`readlink -f $inputFile`
#### generate bigwig track
getbw() {
	bw=$1
	c=$2
	track_name=`basename $bw .bw`
	#c=`colour`
	echo "
track $track_name
bigDataUrl $bw
shortLabel $track_name
longLabel $track_name
maxHeightPixel 40:40:11 
color $c
visibility full
smoothingWindow 2
windowingFunction mean
type bigWig 0 200
autoScale on
viewLimits 0:20

	"

}


getbb() {
	bb=$1
	c=$2
	track_name=`basename $bb .bb`
	#c=`colour`
	echo "
track $track_name
bigDataUrl $bb
shortLabel $track_name
longLabel $track_name
maxHeightPixel 40:40:11 
color $c
type bigBed
visibility dense
viewLimits 0:20
	"

}

#### generate hubs
#longLabel=`pwd  | xargs basename`
#cat > hub.txt <<'EOF' 
echo "hub $title
shortLabel $title
longLabel $comments
genomesFile genomes.txt
email vm11@sanger.ac.uk
descriptionUrl desc.html" > hub.txt

#### trackDb file

cat <<EOF > genomes.txt
genome $genome
trackDb $genome/trackDb.txt
EOF

#### Copy files to genome

mkdir $genome
#cp *.bw $genome

#### Populate the trackDb


cd $genome
echo -n "" > trackDb.txt
while read line
do
	l=($line)
	file=${l[0]}
	colour=${l[1]}
	
	suffix=`echo $file | awk -F "." ' { print $NF } '`
	if [ "$suffix" == "bw" ]; then
		sms "Processing as bigWig >> $file"
		getbw $file $colour >> trackDb.txt
		cp $cwd/$file .
	fi		
	if [ "$suffix" == "bb" ]; then
		sms "Processing as bigBed >> $file"
		getbb $file $colour >> trackDb.txt
		cp $cwd/$file .
	fi		
	
	
done < $inputFile 

