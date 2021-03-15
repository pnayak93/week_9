# week_9

## RNA Seq

The subsample of RNA_seq files used in the sample data for analysis was the first 96 Paired end read files, corresponding to 48 total samples. We moved these files into the /subsamp directory within the RNA_seq directory. We can use the following code to give a list of all the subsampled fastq sample names to use in the DEseq2 pipeline.
While in the folder with the read fastq files we can do:

```
echo *_R1.fq.gz > RNAprefix.txt
for i in {1..48}; do cut -d " " -f $i RNAprefix.txt ; done > RNAprefix2.txt
cut -b 1-7 RNAprefix2.txt > RNAprefixlist.txt
```

This gives a list of all the subsampled fastq gzipped file sample names.

Then, we want to add in the correct Tissue code from the original RNAseq data text file (RNAseq384_SampleCoding.txt). I wrote a short python script to do this quickly. The code is included in this repository as well, named as recreate_tissuecode.py.

The code is as follows:

```
import pandas as pd
import os
import numpy as np

os.chdir("C:/Users/Pavan Nayak/Documents/UCI/Winter 2021 Q11/ECO_EVO_283")
shortRNA = pd.read_table('RNAprefixlist.txt',header=None,names=['FullSampleName'])
shortRNA["TissueCode"] = np.nan

for i in range(0,len(shortRNA)):
    pre = str(shortRNA.iloc[i,0])
    tcode=pre[5]
    shortRNA.iloc[i,1] = tcode

shortRNA.to_csv('shortRNAlist1.txt', sep='\t')

```
After this is finish, we transfer the counts table from the featurecounts.sh script created last week on the hpc3 to the local computer to begin the DESeq2 pipeline on Rstudio locally. To do this, I used the WinSCP program again and pulled the fly_counts_2.txt file to a local folder, started Rstudio, and setwd() to the folder containing the fly_counts_2.txt file and the shortRNAlist1.txt file.


from here, I did some data formatting steps in R to remove the "aligned_" and ".sort.bam" from the column names which were added when creating the bam files in the alignment steps. These are done locally on Rstudio.


```
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

```

Then we can check that the sample info names from the shortRNAlist1.txt file and the formatted sample names in the data frame are the same, and generate the DESeq Matrix using the following code block.

```
#check that the columns are the same
cbind(temp,sampleInfo$FullSampleName,temp == sampleInfo$FullSampleName)

#Generate the DESeq Matrix
dds = DESeqDataSetFromMatrix(countData=countdata, colData=sampleInfo, design=~TissueCode)
dds <- DESeq(dds)
res <- results( dds )
res
```

Then, we make the different plots:

MA plot:

```
#This step is to remove the "Error in plot.new() :  figure margins too large" Error
par(mar = rep(2, 4))

#MA plot
plotMA( res, ylim = c(-1, 1) )
```

![MA Plot](/Figures/MAPlot.png)

Gene Dispersion plot:

```
#Dispersion plot
plotDispEsts( dds )
```
![Gene Dispersion Plot](/Figures/Dispersion_Plot.png)

P_Value_Histogram:

```
#Histogram
hist( res$pvalue, breaks=20, col="blue" )
```
![P Value Histogram](/Figures/p_value_histogram.png)

Distance matrix heatmap:

```
#log transform
rld = rlog( dds )
head( assay(rld) )

#create distance matrix
sampleDists = dist( t( assay(rld) ) )
sampleDistMatrix = as.matrix( sampleDists )
rownames(sampleDistMatrix) = rld$TissueCode
colnames(sampleDistMatrix) = NULL



#generate distance matrix heatmap
library( "gplots" )
library( "RColorBrewer" )
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours)
```
![Distance Heatmap](/Figures/distance_heatmap.png)

PCA plot by tissue type:

```
#PCA plot by tissue type
print( plotPCA( rld, intgroup = "TissueCode"))
library( "genefilter" )
```
![PCA plot](/Figures/PCA_plot.png)

Differential Expression Gene Heatmap:

```
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 35 )
heatmap.2( assay(rld)[ topVarGenes, ], scale="row", trace="none", dendrogram="column", col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))
```

![Differential Expression Gene Heatmap](/Figures/DEG_heatmap.png)

Volcano Plot:

```
#Volcano Plot
 alpha <- 0.01 # Threshold on the adjusted p-value
 cols <- densCols(res$log2FoldChange, -log10(res$pvalue))
 plot(res$log2FoldChange, -log10(res$padj), col=cols, panel.first=grid(),
      main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
      pch=20, cex=0.6)
 
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")
```
![Volcano Plot](/Figures/Volcano_plot.png)


# ATAC Seq

To convert aligned bam files into bigwig files for use on ucsc genome browser tracks etc, we can use the following bash script (bigwigconvert.sh) and run as slurm job on hpc3
assuming the file is in the same folder as all the RG.bam files created in prior steps.

I first created an output directory using `mkdir bed` and then run the job using `sbatch bigwigconvert.sh`

It is not parallelized, but the job runs fast enough without parallelization.

```
#!/bin/bash
#SBATCH --job-name=bigwigconv      ## Name of the job.
#SBATCH -A ecoevo283         ## account to charge
#SBATCH -p standard          ## partition/queue name
#SBATCH --array=1         ## number of tasks to launch, given hint below wc -l $file is helpful
#SBATCH --cpus-per-task=12    ## number of cores the job needs, can the


module load bwa/0.7.17
module load samtools/1.10
module load bcftools/1.10.2
module load python/3.8.0
module load java/1.8.0
module load gatk/4.1.9.0
module load picard-tools/1.87
module load bamtools/2.5.1        # bamtools merge is useful
module load freebayes/1.3.2       # fasta_generate_regions.py is useful
module load vcftools/0.1.16
module load ucsc-tools/v393
module load bedtools2

ref=/data/class/ecoevo283/pnayak/seq_symlinks/ref/dmel-all-chromosome-r6.13.fasta
ATACfiles=/data/class/ecoevo283/pnayak/seq_symlinks/ATACseq/bam
ATACoutput=/data/class/ecoevo283/pnayak/seq_symlinks/ATACseq/bam/bed

cat $ref.fai | head -n 7 | awk '{printf("%s\t0\t%s\n",$1,$2)}' > $ATACoutput/major.bed

for i in *.RG.bam;
do
mybam=$i
name='echo $i | cut -b 11-17'
samtools view -b -L $ATACoutput/major.bed $mybam > ${ATACoutput}/$mybam.maj
Nreads=`samtools view -c -F 4 ${ATACoutput}/$mybam.maj`
Scale=`echo "1.0/($Nreads/1000000)" | bc -l`
samtools view -b ${ATACoutput}/$mybam.maj | genomeCoverageBed -bg -scale $Scale -ibam - > $ATACoutput/$mybam.coverage
bedSort $ATACoutput/$mybam.coverage $ATACoutput/$mybam.sort.coverage
bedGraphToBigWig $ATACoutput/$mybam.sort.coverage $ref.fai $ATACoutput/$mybam.bw
done

```

Then, we can use WinSCP to transfer the bigwig (.bw) files to the desktop so they can be uploaded to Cyverse Discovery Environment, as shown below:

![Cyverse Discovery Environment Bigwig File Upload](/ATACseq_genome_browser/Cyverse_upload.PNG)

After upload to Cyverse, we can generate links and use those links to upload the bigwig files as custom tracks to the UCSC genome browser, as shown below:

![UCSC genome browser bigwig upload](/ATACseq_genome_browser/genome_browser_upload.PNG)

Then, moving back to the genome browser, we can change the track representation from "dense" to "full" to visualize the read peaks on genomic regions:

![UCSC genome browser ATACseq peaks](/ATACseq_genome_browser/genome_browser_peaks.PNG)


