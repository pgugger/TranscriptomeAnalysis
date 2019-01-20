## EXERCISE 4: Exploring gene annotations and functional enrichment analysis

After obtaining a list of interesting genes, such as differentially expressed genes or those in coexpression modules, you might be interested to know more about their functions. A common analysis is to test whether certain functions are more common in the list of interesting genes than you would expect from a random sample of all genes included in the analysis. In this exercise, we will explore gene annotation files and test for enrichment of gene ontology functions using the results of our differential expression analyses in Exercise 2. Specifically we will use the results that we aligned to the reference transcriptome using `bowtie2`, primarily because the reference transcriptome  annotation files are already in a useful format. 

### Understanding reference sequence annotation files

Navigate to `~/Workshop/Genomes/Qlobata_transcriptome/` and view each of the annotation files with `less -S`. The ones with "TAIR" in the name represent genes that are orthologous with genes in the model plant *Arabidopsis thaliana* and thus we can use the relevant gene names and functional annotations from that system. [TAIR](https://www.arabidopsis.org/) is a database of genetic molecular biology data for *Arabidopsis*, including gene annotations. `Annotations_TAIR.txt` is a table whose first column has the oak gene ID (*e.g.*, m01oak00002CC-t01.1), followed by the *Arabidopsis* gene ID (*e.g.*, AT2G38540), the length of the protein produced, a variety of details about the quality of the match between the oak and *Arabidopsis*, and finally at the end any details about the gene function. *Arabidopsis* genes also have [Gene Ontology](http://www.geneontology.org/) functional classifications associated with them, and there are in `Annotations_TAIR_GO.txt`. If you view that file you will see that most genes have more than one: the third column indicates the uppermost hierarchical GO classification (*e.g.*, biological_process) and the fifth and sixth columns give a finer annotation code and description. Similarly,  the [Pfam](https://pfam.xfam.org/) files have annotations of molecular functions at the level of protein domains and motifs and associated GO terms. The GTF (and related GFF) file that we used in a previous exercise is also an annotation file, but it provides structural annotation rather than funcitonal annotation. that we used More details about the transcriptome assembly and annotation can be found in the original [publication](https://doi.org/10.1186/s12864-015-1761-4) and [website](http://genomes.mcdb.ucla.edu/OakTSA/).

Copy (`cp`) the functional annotation files to the `~/Workshop/Enrichment/` to have all the files for this exercise in one place. Also, copy `drought.upreg` and `deseq2.drought.results` (based on alignemt to reference transcriptome) from the `DESeq2` folder to the `Enrichment` folder.

### Annotation of a single gene

From your list of significantly differentially expressed genes, which we saved as `drought.upreg` in Exercise 2, select one to explore in detail. Use `grep` to find it in each of the annotation files and learn more about it. For example,

	grep "^m01oak00032cC" Annotations_TAIR.txt

Note that not all genes have annotations, so some will be missing from these files. Also, notice that there are many lines of GO annotation for every one line of TAIR or Pfam annotation. 

After exploring the function of your gene in the annotation files, find it in the reference transcriptome FASTA file. You can use `grep` again but you will want to use the `-A` argument. Look it up with `man grep` to learn how.

Copy your sequence to the clipboard and search for matches with [blastn](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome), [blastx](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastx&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome), and [tblastx](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Translations&PROGRAM=tblastx&PAGE_TYPE=BlastSearch&BLAST_SPEC=). `blastn` searches for matches in nucleotide sequence databases, `blastx` is for identifying potential protein products encoded by a (translated) nucleotide query, and `tblastx` is for identifying nucleotide sequences similar to the query based on their coding potential. Scroll down, click on one of the top results, and see the alignment. From there you can click links to more information. How do the results compare between each method and with the actual annotation in the file?

### Enrichment analysis on the web

Gene set functional enrichment analysis (or singular enrichment analysis) is a way of testing whether your set of interesting genes is enriched for certain functions. In other words, we want to know if some functions are represented more than you would expect by chance. Therefore, we need to extract the functions from a list of interesting genes and compare it to the list of function in all analyzed genes.

If you work with a model organism or a well-studied organism that has abundant reference resources online and standardized gene expression assays, all you need to do is put the list of official gene IDs into the box in the upper left of [this page](http://www.geneontology.org/) (or any similar one) and select your organism. However, most of us work on non-model systems so we must take additional steps. First, we at least need some standardized annotation data, such as GO, which we might ultimately get via similarity with a model system, as we did with oak. Then, we need to be careful about how we define the background expectation. For example, we cannot take our list of *Arabidopsis* IDs and compare it to the full list of *Arabidopsis* IDs that exist for that species. We need to define the list of *Arabidopsis* IDs that were *present in oaks and analyzed in our specific study*. 

Let's try this exercise using a webpage that allows us to define the background list, assuming we have IDs from a model organism. We can do this on [agriGO](http://bioinfo.cau.edu.cn/agriGO/analysis.php) with the TAIR IDs, for example. First, we need to create the list of TAIR IDs for all the interesting genes and then for all the genes included in our `DESeq2` analysis. 

We can use a string of Linux commands to rapidly execute these tasks. Let's start with the interesting (upregulated) genes. Think about how to break down the problem. We want the TAIR IDs that correspond to oak IDs. Those are both in `Annotations_TAIR.txt` in columns one and two. The interesting oak IDs are in `drought.upreg` in the first column. So, we can `cut` the first column of `drought.upreg`, then remove the header by keeping all lines but the first (`tail -n+2`), then find (`grep`) the matches in `Annotations_TAIR.txt`, and finally keep just the second column with the TAIR IDs. We can pipe (`|`) from one step to the next so it can all happen in one line of code. Here is how it might look:
	
	cut -f1 drought.upreg | tail -n+2 | grep -w -f - Annotations_TAIR.txt | cut -f2 | less -S

`grep -w` looks for whole-word matches only, and `grep -f -` tells `grep` to search more than one string (`-f`) and that those strings are the list that we are piping to `grep` from `tail` (this is called "standard input") (`-`). After `grep` we keep only the second column with the TAIR IDs. Examine what each step of the code does as follows:

	cut -f1 drought.upreg | less -S
	
	cut -f1 drought.upreg | tail -n+2 | less -S
	
	cut -f1 drought.upreg | tail -n+2 | grep -w -f - Annotations_TAIR.txt | less -S
	
	cut -f1 drought.upreg | tail -n+2 | grep -w -f - Annotations_TAIR.txt | cut -f2 | less -S
	
Finally, save the result

	cut -f1 drought.upreg | tail -n+2 | grep -w -f - Annotations_TAIR.txt | cut -f2 > drought.upreg.TAIR

Now, **write your own command** to extract all the relevant TAIR IDs for all the genes analyzed in `DESeq2`. Recall that all the `DESeq2` results were saved as `deseq2.drought.results`. However, this time you have to add an initial step to your pipeline to ignore (remove) all the lines with "NA" in `deseq2.drought.results`. Those genes were not tested statistically and thus should not be included as part of our background expectation. Consider `grep -v` as one possible method. The result should be a substantially reduced list from the full set of over 83k genes (or gene fragments) in the reference transcriptome. Save the result to file.

When you are done, go to the [agriGO](http://bioinfo.cau.edu.cn/agriGO/analysis.php) singular enrichment analysis page and select `Arabidopsis thaliana TAIR10` from the species dropdown menu in part 2, and copy the list of upregulated TAIR IDs to the box. In part 3, select `Customized referece [ Example ]` and copy the list of all TAIR IDs in our `DESeq2` analysis (that were not "NA"). In part 4 (Advanced Options) select `Plant GO Slim`.

The results include a table at the bottom showing which GO functional classes are enriched in the differentially expressed genes. You can also generate various graphs showing the hierarchy and relationship of signifcantly enriched GO terms, as well as bar charts summarizing the differences between the GO terms in interesting genes versus all genes. Do these results make sense given the experiment?

### Enrichment analysis "by hand"

If you have non-standard annotations or just like doing things yourself, you can still test for enrichment by command line and `R`. First we need to choose an annotation category, then test whether it is more common than expected by chance in our upregulated genes than in the full list. Let's choose GO:0006950 which is defined as "response to stress". 

We need to calculate numbers to fill a 2x2 contingency table and its marginal totals, for example

			#Upregulated	#NotUpreg	#TOTAL
	#GO:0006950	x	m-x	m
	#NotGO:0006950	k-x	n-(k-x)	n
	#TOTAL		k	m+n-k	m+n

*k* is easy because in our case it will just be the number of genes in our list of TAIR IDs of upregulted genes, which we saved as `drought.upreg.TAIR` (you can count the lines with `wc -l`). *m*+*n* is also easy because that is the number of all genes with TAIR IDs in our study (which you should have created in the previous section). To get *x*, the number of upregulated genes with GO:0006950 annotation, we need to extract the subset of the `Annotation_TAIR_GO.txt` file that matches upregulated TAIR IDs (in `drought.upreg.TAIR`) and then count how many instances (lines) of "GO:0006950" occur. I suggest using `grep` for both steps (piping the results of the first `grep` to the second `grep`. Finally, you can get *m* by repeating this procedure considering all the TAIR IDs in our study. Generate commands to fill in the table and use those to calculate any missing values.

When finished you have all you need to perform the statistical test in `R`. [Fisher's exact test](https://en.wikipedia.org/wiki/Fisher%27s_exact_test) is the standard method. Underlying the test is a hypergeometric distribution. Here is how to generate a *p*-value for observing a value equal to or greater than *x*

	phyper(x-1, m, n, k, lower.tail=FALSE)

You may get a different result than you did on agriGO, but it appears to be due to updates in the GO annotations since we created ours for oak, or that we are using a slightly different version of GO. In any case, this example illustrates the underlying statistical approach that has been automated on various webpages and software packages. 

