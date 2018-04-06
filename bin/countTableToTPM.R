#!/usr/bin/env Rscript


args_return=Sys.getenv("args_return")
if(args_return > 0){
    message("
* Required Args:

title: 
gtfFile: gtf file using which the count table was generated
inputFile: count table
sampleFile: file containing <name> <condn> <library type>

SampleFile:  <name> <condn> <library type>
BCOR_shRNA1.star.ReadsPerGene.out.tab	BCOR_shRNA1	BCOR_shRNA	paired-end
BCOR_shRNA2.star.ReadsPerGene.out.tab	BCOR_shRNA2	BCOR_shRNA	paired-end
Control_shRNA1.star.ReadsPerGene.out.tab	Control_shRNA1	Control	paired-end
Control_shRNA2.star.ReadsPerGene.out.tab	Control_shRNA2	Control	paired-end

Note: <name> will not be used in this case.

file0: comparisons.txt
Control BCOR_shRNA

		")
    quit(status = 2)
}

title = Sys.getenv("title")
gtfFile=Sys.getenv("gtf")
countTable=Sys.getenv("inputFile")
sampleMatrix = Sys.getenv("sampleFile")

message("Loading libraries")
suppressMessages({
    library(DESeq2)
    library(rtracklayer)
    library(GenomicFeatures)
    library(WriteXLS)
})
# title = "trial1"
# gtfFile = "ensembl_75_transcriptome-1000Genomes_hs37d5.gtf"
# countTable = "featurecount_table.txt"
# sampleMatrix = "samples.txt"



" Read sample matrix "
" -- Note: Create a contrast matrix on your own and parse it in"
sample.matrix = read.table(sampleMatrix)
colnames(sample.matrix) = c("bam","names", "condition", "libtype")
rownames(sample.matrix) = sample.matrix$names
sample.matrix = sample.matrix [, c("condition","libtype") ]
#head(sample.matrix)
" conditions "
conditions = as.character(unique((sample.matrix$condition)))
message("[ Info ] Conditions =  ", paste0(conditions, sep=" "))
#sample.matrix

" Read in GTF file "
gtfData <- import.gff2(gtfFile)
annot_data <- as.data.frame(gtfData)
annot_data <- annot_data[ annot_data$type == "gene", ]
annot_data <- data.frame("gene_name" = annot_data$gene_name,
                         "gene_id" = annot_data$gene_id,
                         "coords" = paste0(annot_data$seqnames, ":", annot_data$start, "-", annot_data$end), 
                         "strand" = annot_data$strand,
                         "source" = annot_data$source
)
annot_data = unique(annot_data)

" -- Convert to txdb"
txdb <- makeTxDbFromGRanges(gtfData)
" -- Convert to genes GRList"
genes_txdb <- exonsBy(txdb, by = "gene")
" -- Merge exonic regions per gene "
genes_fusetxdb <- sapply(genes_txdb, reduce)
" -- Computing length for the genes" 
geneLengths <- sapply(genes_fusetxdb, function(x) sum(width(x)))
" -- gene lengths available for: "
genesToConsider = as.character(names(geneLengths))


" Read count table "
x = read.table(countTable, head=T, row.names=1)
#tail(counttable)
#dim(counttable)
genes.table <- x [ ! grepl("__", rownames(x)), ]
nonGene.table <- x [ grepl("__", rownames(x)), ]
counttable = genes.table
"--> Filter out genes that you can not obtain exon size information"
cat(" Total genes for which length available: ", length(genesToConsider))
cat(" Total genes in count table : ", nrow(counttable))
commonGenes <- genesToConsider [ genesToConsider %in% rownames(counttable)]
cat(" Total genes in common : ", length(commonGenes))
counttable <- counttable [ commonGenes, ]
geneLengths <- geneLengths [ commonGenes ]


colData <- data.frame(
    "condition" = colnames(counttable),
    "libtype" = "paired-end"
)
rownames(colData) <- colnames(counttable)
dds <- DESeqDataSetFromMatrix(countData = counttable, colData = colData, design = ~ 1)
dds <- estimateSizeFactors(dds)
geneCounts <- counts(dds, normalized = FALSE)
normCounts <- counts(dds, normalized = TRUE)

"Calculate TPMS"
geneLengths.kb <- geneLengths [ rownames(normCounts)] / 10^3 # geneLengths in kb
tpm <- normCounts/geneLengths.kb  # Normalise for gene lengths
dataNorm <- colSums(tpm)/10^6 # The column/expt wise normalised is on the geneLength-normalised-matrix
#tpm <- tpm/dataNorm
for(i in colnames(tpm)){
    message("Normalising for ", i)
    tpm[,i] = tpm[,i]/dataNorm[i]
}

" Write TPM to csv"
write.table(tpm, file=paste(title,".tpm.tsv",sep=""), sep="\t", quote=F)
write.csv(tpm, file=paste(title,".tpm.csv",sep=""))

"Calculating FPKMs"
dataNorm <- colSums(normCounts)/10^6
rpkm <- normCounts/geneLengths.kb
for(i in colnames(rpkm)){
    message("Normalising for ", i)
    rpkm[,i] = rpkm[,i]/dataNorm[i]
}
" Write RPPM to csv"
write.table(rpkm, file=paste(title,".rpkm.tsv",sep=""), sep="\t", quote=F)
write.csv(rpkm, file=paste(title,".rpkm.csv",sep=""))




#gene = "ENSG00000000938"
# > head(counttable)
# BCOR_shRNA1 BCOR_shRNA2 Control_shRNA1 Control_shRNA2 Control2_shRNA1 Control2_shRNA2 Ring1_RNF2_shRNA1
# ENSG00000223972.5           1           0              0              0               0               0                 0
# ENSG00000227232.5           1           0              0              0               1               0                 1
# ENSG00000278267.1           0           0              0              0               0               0                 0
# ENSG00000243485.5           0           0              0              0               2               0                 2
# ENSG00000284332.1           0           0              0              0               0               0                 0
# ENSG00000237613.2           0           0              0              0               0               0                 0
# Ring1_RNF2_shRNA2
# ENSG00000223972.5                 0
# ENSG00000227232.5                 0
# ENSG00000278267.1                 0
# ENSG00000243485.5                 1
# 
# > head(sample.matrix)
# condition    libtype
# BCOR_shRNA1     BCOR_shRNA paired-end
# BCOR_shRNA2     BCOR_shRNA paired-end
# Control_shRNA1     Control paired-end
# Control_shRNA2     Control paired-end
# Control2_shRNA1   Control2 paired-end
# Control2_shRNA2   Control2 paired-end
# 
# > head(annot_data)
# gene_name           gene_id           coords strand  source
# 1     DDX11L1 ENSG00000223972.5 chr1:11869-14409      +  HAVANA
# 2      WASH7P ENSG00000227232.5 chr1:14404-29570      -  HAVANA
# 3   MIR6859-1 ENSG00000278267.1 chr1:17369-17436      - ENSEMBL
# 4 MIR1302-2HG ENSG00000243485.5 chr1:29554-31109      +  HAVANA
# 5   MIR1302-2 ENSG00000284332.1 chr1:30366-30503      + ENSEMBL
# 6     FAM138A ENSG00000237613.2 chr1:34554-36081      -  HAVANA
# 
#res.results = getDESeq2(counttable, sample.matrix, annot_data, "Control", "BCOR_shRNA")


