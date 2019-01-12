library(methylKit)
library(genomation)
outpt_rdata = "methylkit.RData"
sampleFile = "bs.test.txt"
qcut = 0.05
percCutOff = 10
organism = "mm10"
ncores = 1




# Read sample file
x = read.table(sampleFile)
colnames(x) <- c("files", "names", "treatment")
file.list=as.list(as.character(x$files))
sampleids = as.list(as.character(x$names))
t = as.list(x$treatment)

# Create methylkit objects
methRawList_data=methRead(file.list,
		sample.id=sampleids,
		assembly=organism,
		treatment=as.numeric(x$treatment),
		context="CpG"
		)


methTabixList_data=methRead(file.list,
              sample.id=sampleids,
              assembly=organism,
              treatment=t,
              context="CpG",
			  dbtype = "tabix",
              dbdir = "methylDB"
              )

save.image(outpt_rdata)

# Filter by sample cov

meth_filt=filterByCoverage(methRawList_data,lo.count=10,lo.perc=NULL, hi.count=NULL,hi.perc=99.9)

# Unite data

meth_unite=unite(meth_filt, destrand=FALSE)

# Unite by groups
meth_condn = unite(meth_filt, min.per.group=1L)

## Unite
# Tiling

tiles=tileMethylCounts(meth_filt,win.size=1000,step.size=1000)
head(tiles[[1]],3)

tiles_filt = filterByCoverage(tiles,lo.count=10,lo.perc=NULL, hi.count=NULL,hi.perc=99.9)

tiles_unite = unite(tiles_filt, destrand = FALSE)


# Correlation, plotting and pca

pdf("corr.pdf")
getCorrelation(meth_unite, plot=T)
clusterSamples(meth_unite, dist="correlation", method="ward", plot=TRUE)
dev.off()

hc = clusterSamples(meth_unite, dist="correlation", method="ward", plot=FALSE)

PCASamples(meth_unite, screeplot=TRUE)

PCASamples(meth_unite)

# diff meth
save.image(outpt_rdata)
diffmeth = calculateDiffMeth(meth_unite, overdispersion="MN",test="F", adjust = "BH", mc.cores=1)
save.image(outpt_rdata)
diffmeth_tiles = calculateDiffMeth(tiles_unite, overdispersion="MN", test="F", adjust = "BH", mc.cores=ncores)
save.image(outpt_rdata)

# Get diffmeth data

diffmeth_res = getMethylDiff(diffmeth,difference=percCutOff,qvalue=qcut) # percCutOff= 0-100 ; qcut = 0 - 1; type = hyper|hypo|all 
diffmeth_tiles = getMethylDiff(diffmeth_tiles,difference=percCutOff,qvalue=qcut)



# Annotate regions
library(genomation)
# read the gene BED file
gene.obj=readTranscriptFeatures(system.file("extdata", "refseq.hg18.bed.txt", 
                                            package = "methylKit"))


temp = readTranscriptFeatures("../../DATA/temp.bed")