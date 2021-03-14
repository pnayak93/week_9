R
library( "DESeq2" )
sampleInfo = read.table("shortRNAlist1.txt")
sampleInfo$FullSampleName = as.character(sampleInfo$FullSampleName)

countdata = read.table("fly_counts_2.txt", header=TRUE, row.names=1)
countdata = countdata[ ,6:ncol(countdata)]
temp=colnames(countdata)
temp = colnames(countdata)
temp = gsub("aligned_","",temp)
temp = gsub(".sort.bam","",temp)
colnames(countdata) = temp

cbind(temp,sampleInfo$FullSampleName,temp == sampleInfo$FullSampleName)

dds = DESeqDataSetFromMatrix(countData=countdata, colData=sampleInfo, design=~TissueCode)
dds <- DESeq(dds)
res <- results( dds )
res

#This step is to remove the "Error in plot.new() :  figure margins too large" Error
par(mar = rep(2, 4))

plotMA( res, ylim = c(-1, 1) )
plotDispEsts( dds )
hist( res$pvalue, breaks=20, col="grey" )
rld = rlog( dds )
head( assay(rld) )
sampleDists = dist( t( assay(rld) ) )
sampleDistMatrix = as.matrix( sampleDists )
rownames(sampleDistMatrix) = rld$TissueCode
colnames(sampleDistMatrix) = NULL
library( "gplots" )
library( "RColorBrewer" )

#generate distance matrix heatmap
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours)

#PCA plot by tissue type
print( plotPCA( rld, intgroup = "TissueCode"))
library( "genefilter" )


#Volcano Plot
alpha <- 0.01 # Threshold on the adjusted p-value
cols <- densCols(res$log2FoldChange, -log10(res$pvalue))
plot(res$log2FoldChange, -log10(res$padj), col=cols, panel.first=grid(),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.6)

abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")