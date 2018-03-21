#!/usr/bin/env Rscript

library(getopt)

args_return=Sys.getenv("args_return")
if(args_return > 0){
    message("Help: 
		title: 
		gtfFile: gtf file using which the count table was generated
		inputFile: count table
		sampleFile: file containing <name> <condn> <library type>
		
		")
    quit(status = 2)
}

title = Sys.getenv("title")
message(title)
gtfFile=Sys.getenv("gtf")
countTable=Sys.getenv("inputFile")
sampleMatrix = Sys.getenv("sampleFile")

message("Loading libraries")
suppressMessages({
    library(DESeq2)
    library(rtracklayer)
    library(GenomicFeatures)
    # library("preprocessCore")
    # library(ggplot2)
    # library(ggthemes)
    # library(reshape2)
})

#gtfFile = "~/nfs_vijay/GENOMES/SANGER/GRCm38/ensembl_84_transcriptome-GRCm38.gtf"
#countTable = "counts.txt"
#sampleMatrix = "samples.GB_Imran.deseq.txt"
#c1="brainControl"
#c2="brainTumour"
#
" Read contrast matrix "
" -- Note: Create a contrast matrix on your own and parse it in"
sample.matrix = read.table(sampleMatrix)
colnames(sample.matrix) = c("bam","names", "condition", "libtype")
rownames(sample.matrix) = sample.matrix$names
sample.matrix = sample.matrix [, c("condition","libtype") ]
#head(sample.matrix)
" conditions "
conditions = as.character(unique((sample.matrix$condition)))
message("[ Info ] Conditions =  ", paste0(conditions, sep=" "))
sample.matrix

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
genesToConsider = as.character(annot_data$gene_id)

" Read count table "
x = read.table(countTable, head=T, row.names=1)
#tail(counttable)
#dim(counttable)
genes.table <- x [ ! grepl("__", rownames(x)), ]
nonGene.table <- x [ grepl("__", rownames(x)), ]
counttable = genes.table
"--> Filter out genes that you can not obtain exon size information"
counttable <- counttable [ genesToConsider, ]
dim(counttable)
" Read in GTF file "
txdb <- makeTxDbFromGFF(gtfFile, format = "gtf", circ_seqs = character())
ebg <- exonsBy(txdb, by="gene")
genesToConsider <- names(ebg)





res.filt <- function(res.data, fc, p.adj, annot_data, fpkmMatrix) {
    temp = as.data.frame(subset(res.data, abs(log2FoldChange) >= fc & padj <= p.adj ))
    x = rownames(temp)
    x = cbind(x, temp)
    colnames(x) = c("gene", colnames(temp))
    x = as.data.frame(x)
    x = merge(annot_data, x, by.x = "gene_id",by.y="gene", all.y = T)
    fpkmMatrix <- fpkmMatrix[, ! grepl("gene_name", colnames(fpkmMatrix)) ]
    x = merge(x, fpkmMatrix, by="gene_id", all.x=T)
    x = x [ order(x$log2FoldChange, x$padj), ]
    return(x)
}

res.lfcThreshold <- function(dds, fc, p.adj, annot_data, fpkmMatrix) {
    res.data <- results(dds, lfcThreshold = fc, alpha = p.adj, format = "DataFrame")
    temp = as.data.frame(res.data)
    temp <- subset(temp, (abs(log2FoldChange) >= fc & padj <= p.adj))
    x = rownames(temp)
    x = cbind(x, temp)
    colnames(x) = c("gene", colnames(temp))
    x = as.data.frame(x)
    x = merge(annot_data, x, by.x = "gene_id",by.y="gene", all.y = T)
    fpkmMatrix <- fpkmMatrix[, ! grepl("gene_name", colnames(fpkmMatrix)) ]
    x = merge(x, fpkmMatrix, by="gene_id", all.x=T)
    x = x [ order(x$log2FoldChange, x$padj), ]
    return(x)
}

medianNormCounts <- function(dds, colData){
    normCounts = counts(dds, normalize = T)
    conditions = unique(as.character(colData$condition))
    countMatrix = c()
    for(cond in conditions){
        reqCols = rownames(colData[colData$condition == cond,])
        temp = rowMedians(normCounts [, reqCols])
        countMatrix = cbind(countMatrix, temp)
    }
    rownames(countMatrix) = rownames(normCounts)
    colnames(countMatrix) = conditions
    return(countMatrix)
}

getDESeq2 <- function(counttable, sample.matrix, c1, c2){
    
    " Extract data corresponding to the conditions"
    colData = subset(sample.matrix, condition %in% c(c1,c2))
    colData$condition = droplevels(colData$condition, c1) # Drop unwanted levels
    colData$condition = relevel(colData$condition, ref = c1) # Keep first condition as the reference level
    
    " Filter count data for low count data, reads/gene > 0"
    countdata <- counttable
    message(rownames(colData))
    countdata = as.matrix(countdata [ , rownames(colData)])
   
    rmax = apply(countdata, 1, max)
    message("[ Info ] Total genes = ", nrow(countdata))
    countdata = countdata [ rmax >= 1, ]
    message("[ Info ] Genes with reads = ", nrow(countdata))
    
    
    " -- Check if the row order in contrast and columns in count data are in the same order"
    if( ! all(rownames(colData) == colnames(countdata)) ){
        message("[ Error ] Total rows in filtered contrast matrix != Total columns in count table")
        message("[ Error ] Aborting")
        save.image("deseq_error.RData")
        quite(save="no")
    } else {
        message("[ Info ] All seems correct. Proceeding forward ...")
    }
    
    " Create DESeq object"
    condition = colData$condition
    condition = relevel(condition, ref = c1)
    dds <- DESeqDataSetFromMatrix(countData = countdata, colData = colData, design = ~ condition)
    #dds$condition <- factor(dds$condition, levels=c(c1,c2))
    print(condition)
    
    " Perform DESeq "
    dds <- DESeq(dds)
    res <- results(dds, alpha=0.05, format = "DataFrame")
    head(res)
    medCounts = medianNormCounts(dds, colData)
    outpt = merge(res, medCounts, by = "row.names")
    #summary(res)
    " Extract results"
    genesInRes = as.character(rownames(res))
    n = paste0(c1,"__",c2)
    " -- Data from DESeq "
    return(outpt)
}

message("Here ....")

res.results = getDESeq2(counttable, sample.matrix, "Control", "BCOR_shRNA")




save.image("temp.RData")
q(status = 0)
library(WriteXLS)
WriteXLS(res.results, ExcelFileName=paste(title,".deseqAll.xlsx",sep=""), row.names=F)

# 
# " Write results "
# fold = paste(c1, "VS", c2, sep="_")
# system(paste("mkdir ", fold), intern = F)
# library(WriteXLS)
# WriteXLS(res.results, ExcelFileName=paste(fold,"/deseqAll.xlsx",sep=""), row.names=F)
# 
# " Plots "
# plot.volcano <- function(res, lfc, pval){
#     tab = data.frame(logFC = res$log2FoldChange, negLogPval = -log10(res$pvalue))
#     head(tab)
#     par(mar = c(5, 4, 4, 4))
#     plot(tab, pch = 16, cex = 0.6, xlab = expression(log[2]~fold~change), ylab = expression(-log[10]~pvalue))
#     signGenes = (abs(tab$logFC) > lfc & tab$negLogPval > -log10(pval))
#     points(tab[signGenes, ], pch = 16, cex = 0.8, col = "red") 
#     abline(h = -log10(pval), col = "green3", lty = 2) 
#     abline(v = c(-lfc, lfc), col = "blue", lty = 2) 
#     mtext(paste("pval =", pval), side = 4, at = -log10(pval), cex = 0.8, line = 0.5, las = 1) 
#     mtext(c(paste("-", lfc, "fold"), paste("+", lfc, "fold")), side = 3, at = c(-lfc, lfc), cex = 0.8, line = 0.5)
# }
# 
# pdf(paste(fold,"/deseq.pdf",sep=""),5,5)
# plotMA(res, frame=F, alpha=0.05)
# plot.volcano(res, 1, 0.05)
# plotDispEsts(dds)
# dev.off()
# 
# 
# " Save R image"
# save.image(paste(fold,"/deseq.RData",sep=""))
