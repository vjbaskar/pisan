# Given an tsv/csv file get the required columns

library(getopt)
mconf_installdir=Sys.getenv("mconf_installdir")
source(paste0(mconf_installdir,"/src/basic_functions.R"))

args_return=Sys.getenv("args_return")
if(args_return > 0){
    message("
* Required Args:

-t|--title
-s|--sampleFile: file to be filtered
-c|--comments: colnames to be filtered SAMPLEID,SAMPLENAME
-0|--cmd0: Input file type [tsv,csv]
-1|--cmd1: Output file type [tsv,csv]
-O|--outpt: Output file
		")
    quit(status = 2)
}

title = Sys.getenv("title")
sampleFile = Sys.getenv("sampleFile")
cols = Sys.getenv("comments")
inputType=Sys.getenv("cmd0")
outputType=Sys.getenv("cmd1")
outpt=Sys.getenv("outpt")

if(inputType=="tsv"){
	info("Reading file as tsv:")
	x = read.table(sampleFile, sep="\t", header=T)
}
if(inputType=="csv"){
	info("Reading file as csv:")
	x = read.csv(sampleFile, sep="\t", header=T)
}

info("Filtering required columns:")
reqCols = strsplit(cols, ",")[[1]]
x = x[,reqCols]




if(outputType=="tsv"){
info("Writing file as tsv:")
	x = write.table(x, file=outpt, sep="\t", quote=F ,row.names=F)
}
if(inputType=="csv"){
info("Writing file as csv:")
	x = read.csv(x, file=outpt,row.names=F	)
}