# basic commands for an R program in pisan

library(getopt)
mconf_installdir=Sys.getenv("mconf_installdir")
source(paste0(mconf_installdir,"/src/basic_functions.R"))

args_return=Sys.getenv("args_return")
if(args_return > 0){
    message("
* Required Args:

-t|--title
-s|--sampleFile: master sample file csv
-c|--comments: colnames eg. SAMPLEID,SAMPLENAME

		")
    quit(status = 2)
}

title = Sys.getenv("title")
sampleFile = Sys.getenv("sampleFile")
cols = Sys.getenv("comments")


