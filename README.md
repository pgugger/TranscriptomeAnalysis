# *Bioinformatics for Transcriptome Analysis*
Exercises for transcriptome analysis workshop taught at UNAM-Morelia, Mexico in January 2019. Topics include Linux basics, Illumina data quality assessment, quantifying gene expression levels with RNA-Seq, testing for differential expression among experimental treatments, and weighted gene coexpression network analysis. 

## Required software

### On local computer

Windows: [PuTTY](https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html) and FTP client (e.g., [FileZilla](https://filezilla-project.org/) or [WinSCP](https://winscp.net/eng/download.php))

MacOs: FTP client (e.g., [FileZilla](https://filezilla-project.org/))

### On server

[bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

[samtools](http://www.htslib.org/) 

[fastq_screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/)

[fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)  (see "Requirements" section for dependencies)

[R]( https://www.r-project.org/) (optionally with [RStudio](https://www.rstudio.com/))

Then in R, run the following commands to install necessary packages (and any dependencies):
	
	install.packages(c("gplots", "pheatmap", "RColorBrewer", "igraph"))
	if (!requireNamespace("BiocManager"))
	  install.packages("BiocManager")
	BiocManager::install()
	BiocManager::install(c("qvalue", "DESeq2", "WGCNA", "Rsubread", "genefilter"))

[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)

[STAR](https://github.com/alexdobin/STAR)
