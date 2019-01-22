## EXERCISE 2: Differential gene expression analysis

Differential gene expression analysis allows us to test for statistically significant differences in gene epression levels among tissues, experimental treatments, or other designs. There are a number of popular methods, including `edgeR` and `DESeq2`. Here, we will use `DESeq2` with an example data set from [Gugger *et al.* (2017)](https://doi.org/10.1093/treephys/tpw122). A very [detailed tutorial](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) from the `DESeq2` developers is also available to help determine the best analysis for your own projects. 

The study for this workshop represents a small gene expression (RNA-Seq) experiment performed on oak (*Quercus*) seedlings in a greenhouse. The overall question is "what genes are involved in response to water stress?" An experiment was performed by exposing nine seedlings to water stress in the greenhouse. RNA was extracted from leaves before the treatment and again after 10 days of imposing water stress (no water at all). The idea will be to compare gene expression differences before versus after the water stress treatment. 

The data set consists of 12 FASTQ files containing raw RNA-Seq read data from an Illumina sequencer – one for each sample before treatment and one for each sample after. In addition, I have provided you with oak reference transcriptome and reference genome sequences in FASTA format (see "Genomes" folder). Much like a reference genome might have a representative sequence of each chromosome, a reference transcriptome contains representative sequences of each mRNA transcript.

### Exploring the data

Let's revisit how many reads are in each file. FASTQ files have four lines per read, so you could count the lines with `wc -l` and then divide by four. However, this task will be tedious with more than a few files, so it is time to introduce the concept of loops. A loop allows you to execute the same command(s) over and over on a series of files (or other inputs). 

Start by deciding what commands you want to run. Here, we might want to `zcat` the gzipped FASTQ file and then pipe to `wc -l` as we did yesterday, and then divide the result by four manually. But, let's be more efficient. For example, we could instead `zcat`, then `grep` only the header line (one of the four lines per read), and then count. `grep` is convenient here because all of the header lines in our example files start with "@SN603" and because `grep` has a count feature built in, `-c`. So, our command might look like this

	zcat MC511_after.fq.gz | grep -c '@SN603'

Try the above command on one file to make sure it works. Now, we can construct the loop. First, we need to identify the variable. We want to run the same command but each time do it on a different file, so the variable is the input file name. Here is the notation in single-line format

	for file in *.fq.gz; do zcat ${file} | grep -c '@SN603'; done

This command can be read as "for each file in the list of files ending in .fq.gz, `zcat` and then count the lines `-c` containing '@SN603' using `grep`." After it finishes with the first file in the list, it will go to the next file until it completes all files fitting that pattern. Note that the variable `file` is defined in the first part of the command. It is a word that has meaning to us, but you could put anything there, such as `for thing in *.fq.gz; do zcat ${thing} | grep -c '@SN603'; done`. Also notice how the variable `file` (or `thing`) is invoked in the command as `${file}` (can also be written as `$file`). This loop is simple so running it in a single line command is fine. 

For more complicated loops or ones you use often, you can modify the notation slightly and save it to file. For example, you could create a text file called count.sh and copy the following text into it.

	#!/bin/bash

	for file in *.fq.gz
	do 
		zcat ${file} | grep -c '@SN603'
	done

The first line is standard notation indicating that it is a Linux/Unix `bash` script. Then, the `for` loop is given on multiple lines instead of one, mostly for clarity. Once you save the file with this text, you have to make it executable

	chmod +x counts.sh
	
Then, to run it

	./counts.sh
	
After you run one version of the loop, you should see the read counts print to screen, in list order. I actually gave you only 10% of the reads in the published data set.

You will find `for` loops helpful for the rest of this workshop and in your own work.

### Aligning and counting reads

There are many short-read aligners available. When aligning to a transcriptome, standard short-read aligners, such as `bwa` and `bowtie2`, will work pretty well. When aligning to a genome, however, it is important to consider that RNA reads may represent sequence from more than one exon without the intervening intron sequence. Thus, there is another class of aligners that are splice aware, some of which utilize a more basic aligner under the hood. Popular splice-aware aligners include [STAR](https://github.com/alexdobin/STAR), [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) (successor to TopHat), and [Rsubread](https://bioconductor.org/packages/release/bioc/html/Rsubread.html), among others. Another consideration is the genomic "scale" of interest. For example, are you interested in differences in expression at the whole-transcript or whole-gene level or at the exon or "genomic feature" level? 

We will start with a simple case of generating read counts from alignment to a transcriptome and count expression level by transcript. Then, we will analyze those data with `DESeq2` ot test for differential expression. If time permits and the computational resources are available, then we will move to a more complicated case of aligning to a genome with a splice-aware aligner and generate counts by genomic feature. You can then use what you learned in DESeq2 to analyze those data as well.

### Alignment to transcriptome, counts by transcript

The advantage of this approach is that it is relatively simple, fast, and will not consume large amounts of computational resources that may be limited during the workshop. I should caution that most RNA-Seq analyses involve more nuance than premitted in this section, some of which we will get to later. 

Before we can align reads to the transcriptome, we need to "index" the reference sequence. 

	bowtie2-build ~/Workshop/Genomes/Qlobata_transcriptome/Qlobata_transcriptome.fasta ~/Workshop/Genomes/Qlobata_transcriptome/Qlobata_transcriptome

If it worked you will end up with a set of bt2 files in the folder with the reference FASTA. `bowtie2` uses these to speed up the alignments.

Generating alignments in the proper format has three main steps: (1) the basic alignment to produce a SAM-format alignment file (`bowtie2`), (2) conversion of the SAM file to a BAM file which is a more efficient binary version of the SAM file (`samtools view -b`), and (3) sorting the BAM file (`samtools sort`). I will lay out the three example commands below, then we will try to combine them into an efficient pipe (`|`) and wrap that in a `for` loop to work through all the samples efficiently. I suggest running through each command with one sample to see what each produces. Then, we can run everything with a loop. 

	#Align with bowtie2
	cd ~/Workshop/RNA-Seq_Data/
	bowtie2 -p 16 --no-unal -x ~/Workshop/Genomes/Qlobata_transcriptome/Qlobata_transcriptome -U MC511_after.fq.gz -S MC511_after.sam 
	
	#Convert SAM to BAM
	samtools view -b MC511_after.sam > MC511_after.bam 
	
	#Sort BAM
	samtools sort --threads 16 MC511_after.bam > MC511_after.sorted.bam
	
See the [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) and [samtools](http://www.htslib.org/doc/samtools.html) manuals for more details. Open the SAM file with `less -S` to look at its structure. The top contains a header, one line per contig in the reference transcriptome. Scroll down 100,000 lines by typing `100000` and then `ENTER` while in `less`. Each line is a read and information about where it aligns and how well. Refer to [this document](http://genome.sph.umich.edu/wiki/SAM) for details. Note that you could actually count how many reads align to a given transcript by `grep`ing the name of the transcript, *e.g.*, `grep -c 'm01oak00002CC' MC5-11_after.sam`. Of course, there is a better way to do so for all transcripts at once.

Notice in the above commands that we had to save to file with each command, but we don't really need any of these except the last sorted BAM. Instead, we could have piped each of the outputs to the next command and only save the final result:

	bowtie2 -p 16 --no-unal -x ~/Workshop/Genomes/Qlobata_transcriptome/Qlobata_transcriptome -U MC511_after.fq.gz | samtools view -b | samtools sort --threads 16 > MC511_after.bam

See how we do not direct the output of `bowtie2` to file using `-S MC511_after.sam`. Instead, it pipes to `samtools view -b`. Similarly, we do not direct those results to file using `> MC511_after.bam` and instead pipe them to `samtools sort`. The ouput is directed to file only once, avoiding the tedium of running each command separately and the mess of generating many intermediate files.

If we have more than a few files, it could still be tedious to run this command over and over, so let's use it to build another `for` loop. Try it on your own before looking at my version below. Don't run yet, because we will add a few more tricks.

Given what you know from above, you're loop might look something like this:

	#Do not run
	for file in *.fq.gz; do bowtie2 -p 16 --no-unal -x ~/Workshop/Genomes/Qlobata_transcriptome/Qlobata_transcriptome -U ${file} | samtools view -b | samtools sort --threads 16 > ${file}.bam; done

That will work fine but notice that the result will be a file with an extension of `.fq.gz.bam`. To avoid endlessly adding extensions, we can loop through the sample names instead of file names (do not run yet):

	#!/bin/bash

	samples="MC511_after
	MC511_before
	MC711_after
	MC711_before
	MK1015_after
	MK1015_before
	MK54_after
	MK54_before
	SVG10_after
	SVG10_before
	SVH10_after
	SVH10_before"

	for sample in $samples; do
		bowtie2 -p 16 --no-unal -x ~/Workshop/Genomes/Qlobata_transcriptome/Qlobata_transcriptome -U ${sample}.fq.gz | samtools view -b | samtools sort --threads 16 > ${sample}.bam
	done

Notice how this script would lead to the file extension being replaced, from `.fq.gz` to `.bam`. Also, note that `sample` and `samples` are different variables. `samples` is a list of the sample names. `sample` is an individual sample name from the list samples. So, the loop substitutes a different sample name each round in the order they appear in the list. 

The final steps to go from alignments to read counts per transcript are to index the BAM files (`samtools index`) and then use `samtools idxstats` to get counts. We can add those into our script too:

	#!/bin/bash

	samples="MC511_after
	MC511_before
	MC711_after
	MC711_before
	MK1015_after
	MK1015_before
	MK54_after
	MK54_before
	SVG10_after
	SVG10_before
	SVH10_after
	SVH10_before"

	for sample in $samples; do 
		bowtie2 -p 16 --no-unal -x ~/Workshop/Genomes/Qlobata_transcriptome/Qlobata_transcriptome -U ${sample}.fq.gz | samtools view -b | samtools sort --threads 16 > ${sample}.bam
		samtools index ${sample}.bam
		samtools idxstats ${sample}.bam > ${sample}.counts
	done

Save this script (*e.g.*, `align_count.sh`), make it executable (`chmod +x align_count.sh`), and run (`./align_count.sh`). It will take about 10-15 minutes to run and as it completes each sample it provides summary information about the fraction of reads that aligned, *etc*.

The results include a sorted BAM file, an associated BAI index file, and a text file with counts per transcript for each sample. View one of the count files with `less`. The first column is the transcript name, the second is the length of the transcript in the reference FASTA, and the third is the number of reads that map to each transcript, and the fourth is the number of unmapped reads (ignore).

The final step is to put them together into a single tab-separated table, with one row per transcript (as is) and columns as counts for each sample, e.g.

	Transcript	MC511_after MC511_before MC711_after	...
	m01oak00001CC           0            0           0
	m01oak00002CC         580         6809         694
	m01oak00003cm      153486        15223       41814
	...

Here is one way to do it while staying in the Linux terminal:
	
	#Cut the third column (counts) from each file and save temporarily
	for file in *.counts; do cut -f3 ${file} > ${file}.temp; done
	
	#Paste all the files containing column 3
	paste *.counts.temp > all.temp
	
	#Paste the names of the transcripts to the above
	paste <(cut -f1 MC511_after.counts) <(cat all.temp) > all2.temp
	
	#Add header row based on sample names in sample.info (assumes same order)
	paste <(echo -e "Transcript") <(cut -f1 ~/Workshop/sample.info | tail -n +2 | tr '\n' '\t') | cat - all2.temp > all.counts 
	
	#Confirm result
	less -S all.counts
	
	#Copy to DESeq2 folder
	cp all.counts ~/Workshop/DESeq2/
	
	#Delete temp files
	rm *.temp
	
There are many other ways to go about this, and you might find it easier to do it `R` or other software.

### Differential expression analysis

DESeq2 is one popular method of testing for differences in expression. It implements a negative binomial generalized linear model (appropriate for count data), which you can learn more about by consulting the original publication or manual. To start, enter into R, load DESeq2, and set your working directory:

	library("DESeq2")
	setwd("~/Workshop/DESeq2/")
	
Load the read count and sample data
	
	all.counts <- read.table("all.counts", header = T, row.names = 1)
	head(all.counts)

	sample.info <- read.table("~/Workshop/sample.info", header=T, row.names=1)                         
	sample.info

From `sample.info`, we need to extract only the variables that we want to test in our model. This depends on your question and design. Our main interest is in testing difference before versus after drought, so we definitely want column 4. We might want to control for sequencing lane (*i.e.* batch effects) and consider differences among populations, which means we would also want columns 2 and 3. Alternatively, we could control just at the individual level since we have repeated measures within individuals: columns 1 and 4. Choose whichever you prefer, but you will get errors if you try to choose all of them at once:

	exp.design <- sample.info[, 2:4]   #or try sample.info[,c(1,4)]
	exp.design

Now we can define the complete input data frame with a built-in DESeq2 function that requires the count data, experimental design, and model as inputs:

	all.input <- DESeqDataSetFromMatrix(countData = all.counts, colData = exp.design, design = ~Lane + Population + Drought)
	#Adjust the design part according to the variables you chose in exp.design

Running the DESeq2 model is very simple.

	output <- DESeq(all.input)

To extract the relevant results from the output, we specify the "contrast" of interest.

	results.drought <- results(output, contrast = c("Drought", "before", "after"))
	
	#You can also save to file
	write.table(as.data.frame(results.drought), file="deseq2.drought.results", sep="\t", quote=F)

Although less interesting, you can also extract other results, such as those comparing mean expression levels across populations (not related to interaction with drought).

	results.pop <- results(output, contrast = c("Population", "MC", "SV"))
	#...and so on for other pairs of populations

If you look at the first few lines of the main drought results (`head(results.drought)`), you will see that there are six columns of information for each transcript: mean baseline expression, log<sub>2</sub>-fold change in expression before versus after treatment (positive values are decreases from before to after, in this case), its standard error, a test statistic, the *p*-value, and finally the *p*-value adjusted for multiple testing. Focus on the log<sub>2</sub>-fold change and the adjusted *p*-value. 

You explore these results with simple queries. For example, how many transcripts significantly (α = 0.01) differ in expression due to the treatment? 

	length(which(results.drought$padj < 0.01))

What fraction of the total number of transcripts is that? Avoid counting the ones with NA for *p*-values because those didn’t have enough alignments to perform statistical tests. 

	length(which(results.drought$padj < 0.01)) / length(na.omit(results.drought$padj))

How many have log<sub>2</sub>-fold change > 3 (or < -3) and *p*<sub>adj</sub> < 0.01? Note that log<sub>2</sub>-fold change of 3 means a fold change of 8, *i.e.*, an change in expression by 8 times.

	#All genes with large, significant changes
	length(which(results.drought$padj < 0.01 & abs(results.drought$log2FoldChange) > 3))
	
	#Downregulated genes (before versus after drought) with large, significant changes
	length(which(results.drought$padj < 0.01 & results.drought$log2FoldChange > 3))
	
	#Upregulated genes (before versus after drought) with large, significant changes
	length(which(results.drought$padj < 0.01 & results.drought$log2FoldChange < -3))
	
Let's save a subset of the results table of the significant up-regulated genes, following the criteria as above. Note that `DESeq2` considers the after-drought samples as the reference in this case (probably because "after" is alphabetically before "before"), so a negative log<sub>2</sub>-fold change means a *decrease in expression from after to before treatment*, or more importantly, an *increase in expression from before to after treatment* (you can confirm this by referring back to your `all.counts` table. We will use these results for Exercise 4:

	drought.upreg <- results.drought[ which(results.drought$padj < 0.01 & results.drought$log2FoldChange < -3), ]
	
	write.table(as.data.frame(drought.upreg), file="drought.upreg", sep="\t", quote=F)
	
You can save other subsets following this format if you like.

Finally, we can make some plots to summarize the results. Here is one showing log<sub>2</sub>-fold change versus gene expression level with the significant ones colored red:

	pdf("DESeq2_BeforeAfter.pdf")
	plotMA(results.drought, ylim=c(-10,10), main="Before versus after drought", alpha=0.01)
	dev.off()
 
You can also make a principal components analysis plot of the overall expression patterns, but the data need to be transformed first to better meet the assumptions of the analysis. Regularized-log (rlog) is one popular transformation:

	transformed.counts <- rlog(output)

	pdf("DESeq2_PCA.pdf")
	plotPCA(transformed.counts, intgroup="Drought")  #Specify a different intgroup to color points by that instead
	dev.off()
	
Before we move on, let's save the transformed counts to file because we will use them tomorrow.

	write.table(assay(transformed.counts), "~/Workshop/WGCNA/transformed.counts", sep="\t", quote=F)

There are many ways you might want to further explore the data. Heatmaps are a popular choice, but it helps to think first about what you want to show. Here are a few possibilities:

	#Heatmap showing distances among samples
	library("gplots")
	library("RColorBrewer")
	dist.transformed.counts <- dist(t(assay(transformed.counts)))
	mat <- as.matrix(dist.transformed.counts)
	hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
	rownames(mat) <- colnames(mat) <- with(colData(all.input), paste(Drought, Population, sep=" : "))
	
	pdf("DESeq2_HeatmapSamples.pdf")
	heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
	dev.off()

	#Heatmap based on 30 most expressed genes
	top.exp <- order(rowMeans(counts(all.input, normalized=F)), decreasing=TRUE)[1:30]
	
	pdf("DESeq2_HeatmapMostExpressed.pdf")
	heatmap.2(assay(transformed.counts)[top.exp,], col = hmcol, Rowv = FALSE, Colv = FALSE, scale="none", dendrogram="none", trace="none", margin=c(10, 6))
	dev.off()

	#Heatmap based on 30 genes with most variable expression
	library("genefilter")
	topVarGenes <- head(order(rowVars(assay(transformed.counts)), decreasing=TRUE), 30)
	
	pdf("DESeq2_HeatmapMostVariable.pdf")
	heatmap.2(assay(transformed.counts)[topVarGenes, ], scale="row", trace="none", dendrogram="column", col = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255))
	dev.off()

#### Gene function

You might wonder what the functions of the differentially expressed genes are or whether certain gene functional classes are more highly represented than expected by change. The former might be as simple as looking up the relevant genes in the annotation files I provided. The latter requires more consideration and we will revisit these in Exercise 4.

#### Interpretation

What conclusion would you reach from the analysis so far? Are there any confounding factors that might also explain the patterns? What follow-up experiment might you want to do?

### Alignment to a genome and counting by feature

Now we will try to align to the genome and count reads by feature (*e.g.*, exon). You can do this whether using a reference genome or transcriptome, but the reference sequence must be annotated and have an annotation file in GTF format. If you have GFF-formatted annotations, you can [convert GFF to GTF](https://github.com/gpertea/gffread). We will use the `STAR` aligner, which is splice aware and finds novel splice junctions using two passes at alignment. Then, we will move to `R` and use the `Rsubread` package to count alignments per feature. `Rsubread` can actually also do alignments, but we will just use the `featureCounts` function.

First, we need to move all the BAM files that we created by mapping to the transcriptome above to make room for the new ones mapped to the genome. Navigate to `~/Workshop/RNA-Seq_Data/`, create a new directory (`mkdir`) to deposit the old files, and move them with `mv`.

	cd ~/Workshop/RNA-Seq_Data/
	mkdir transcriptome
	mv *.bam* transcriptome/
	mv *.counts transcriptome/
	
We will also want the reference genome annotation file (.gtf), so copy it to our working folder that has the FASTQ files.

	cp ~/Workshop/Genomes/Qlobata_genome_v05/Qlobata.reduced.subset.gtf ./
	
Take a look at the GTF file using `less`, so you can understand what it contains. See [here](http://mblab.wustl.edu/GTF22.html) for more details.

Now we are ready for the [STAR](https://github.com/alexdobin/STAR) aligner. First we need to index the genome for use with `STAR`.

	STAR --runThreadN 16 --runMode genomeGenerate --genomeDir /home/pgugger/Workshop/Genomes/Qlobata_genome_v05 --genomeFastaFiles /home/pgugger/Workshop/Genomes/Qlobata_genome_v05/Qlobata.reduced.fasta --sjdbOverhang 50 --sjdbGTFfile /home/pgugger/Workshop/Genomes/Qlobata_genome_v05/Qlobata.reduced.subset.gtf --sjdbGTFtagExonParentTranscript Parent --genomeSAindexNbases 12 --genomeChrBinNbits 11 

Then, we are ready to align reads. However, the STAR aligner requires decompressed FASTQ, so we need to `gunzip` them first. Here is how you can do all the files at once in parallel (rather than serially with `gunzip *.fq.gz`)

	for file in *.fq.gz; do gunzip "${file}" & done
	
Using what you learned about loops earlier today, create a `for` loop to align all the samples (FASTQ files) to the reference genome using the following command as a template
	
	STAR --runThreadN 16 --twopassMode Basic --twopass1readsN -1 --genomeDir ~/Workshop/Genomes/Qlobata_genome_v05 --sjdbOverhang 50 --readFilesIn MC5-11after.fq --outSAMattributes All --outSAMtype BAM SortedByCoordinate --outFileNamePrefix MC511_after

Each sample takes about 2 minutes to run so your loop with 12 samples will take over 20 minutes. In the meantime, read about the options I chose in my example command and think about how they might (or might not) differ for your study. When complete, you will see a number of outputs, which you can read about in the manual. The most important for us now are the BAM files.

When it finishes, enter `R` to generate counts per feature (and then test for differential expression among features) from the BAM files. 

	setwd("~/Workshop/RNA-Seq_Data/")
	library(Rsubread)
	
Load the list of BAM files that we just created with `STAR`.

	bam.files <- list.files(path = ".", pattern = ".bam$", full.names = TRUE)
	bam.files

Now we can use the `featureCounts` function in `Rsubread` to generate read counts per exon.	Notice in the command below that we need to indicate the annotation file (`annot.ext="Qlobata.reduced.subset.gtf"`), specify the type (`isGTFAnnotationFile=T`), specify which feature we are interested in (`GTF.featureType="exon"`), and within what attribute they are contained (`GTF.attrType="transcript_id"`). I also specified a minimum read mapping quality (`minMQS=20`).
	
	feature.counts <- featureCounts(bam.files, annot.ext="Qlobata.reduced.subset.gtf", isGTFAnnotationFile=T, minMQS=20, useMetaFeatures=F, GTF.featureType="exon", GTF.attrType="transcript_id")

Take a quick look at the structure of the output we saved in `feature.counts` when it finishes. 
	
	#See what components there are
	names(feature.counts)

	#View summary statistics
	feature.counts$stat
	
	#Display dimensions of the read count table: exons x samples
	dim(feature.counts$counts)
	
	#Display a piece of the count table
	head(feature.counts$counts)

The most relevant component is `feature.counts$counts`, which is analogous to the `all.counts` table that we generated earlier when working with the transcriptome. Let's extract those and save with the same name (you can choose a different name if you prefer not to overwrite what we did earlier).
	
	all.counts <- feature.counts$counts

Recall, however, that the sample names are the names of the BAM files `colnames(all.counts)` rather than the actual sample names we have been using. So, let's change the column names accordingly. We can easily do so by importing the `sample.info` and using the *row* names, assuming they are in the same order (here they are, which is very important).

	sample.info <- read.table("~/Workshop/sample.info", header=T, row.names=1)                         
	sample.info
	
	colnames(all.counts) <- rownames(sample.info)
	head(all.counts)

Now, I will leave it up to you to analyze the data in `DESeq2` following the example that I gave with the transcriptome-aligned data. The `all.counts` and `sample.info` files here are analagous to the ones with the same names earlier.
	
Any differences in your interpretation of this genomic anlaysis by exon compared to the transcriptome analysis by transcript? Any cases where only certain exons within transcripts/genes are differentially expressed?
