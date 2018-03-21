
if [ $args_return  -gt 0 ]; then
print_stderr "
	++ Mandatory args: title, organism
	
	options: 
	organism: hg19,hg38 etc.
	title: alpha numeric 


"
	exit 2
fi

p=`pwd`
BASE_DIR=$mconf_installdir
name=$title
genome=$organism

#### Generate random colours
colour() {
	r=`shuf -i 1-255 -n 1`
	g=`shuf -i 1-255 -n 1`
	b=`shuf -i 1-255 -n 1`
	echo "$r,$g,$b"
}

colour() {
	c=`randomLines $BASE_DIR/../data/palette1.txt 1 stdout`
	echo $c
}

#### generate bigwig track
getbw() {
	bw=$1
	track_name=`basename $bw .bw`
	c=`colour`
	echo "
track $track_name
bigDataUrl $bw
shortLabel $track_name
longLabel $track_name
maxHeightPixel 40:40:11 
color $c
visibility full
smoothingWindow 5
windowingFunction mean
type bigWig 0 200
autoScale off
viewLimits 0:20

	"

}

#### generate hubs
longLabel=`pwd  | xargs basename`
#cat > hub.txt <<'EOF' 
echo "hub $name
shortLabel $name
longLabel $longLabel
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
cp *.bw $genome

#### Populate the trackDb

cd $genome
for file in `ls *.bw`
do
	getbw $file
done > trackDb.txt

