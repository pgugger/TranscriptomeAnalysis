## EXERCISE 3: Weighted gene coexpression network analysis

[Weighted gene coexpression network analysis (WGCNA)](https://doi.org/10.1186/1471-2105-9-559) is a systems biology method for describing the patterns of gene expression correlation among genes. Clusters of highly correlated genes ("modules") are identified with hierarchical clustering combined with user-defined rules for defining discrete modules from the clustering dendrogram. These gene modules may reflect genes that interact in expression responses or genes that coincidentally respond or are expressed similarly. An advantage of the approach is that it does not require any prior annotation information to discover putative functional interactions among genes. The expression pattern captured by modules can then be summarized using principal components analysis to define the "eigengene" or hub gene. This summary allows the expression patterns of each module to be associated with traits or other variables of interest. Such an analyses might provide candidate genes for involvement in traits or responses to environmental variables. `WGCNA` provides an extensive package of R functions to perform nuanced analyses for diverse and large data sets. 

In this exercise, we will work through an example of WGCNA using the same RNA-Seq data from oak tree seedlings that we used yesterday for DESeq2 [(Gugger *et al.* 2017)](https://doi.org/10.1093/treephys/tpw122). This tutorial borrows heavily from the [excellent tutorials](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html) created by the software developers (Langfelder & Horvath). The `R` code can be complicated at times, but I hope that you will get a sense for the general workflow and important considerations when employing this method. For additional details or potential analyses, I highly recommend their tutorials. After identifying expression modules and summarizing them by eigengene, we will test for signifcant correlations with the drought treatment (encoded as a dummy variable) as well as latitude, longitude, elevation, and five climate variables measured at the seed source for the plants used in the experiment. Note that `WGCNA` does not demand a data set with experimental treatments as we have here.

### Preparing the data and environment

From the Linux terminal enter `R` to begin: `srun -p ALL --mem 20G --cpus-per-task 10 --pty R`. In `R`, set the working environment:

	setwd("~/Workshop/WGCNA/")
	library(WGCNA)
	enableWGCNAThreads(nThreads = 10)  #May not work with RStudio
	options(stringsAsFactors = FALSE)

We will use the rlog-transformed counts from `DESeq2` as input to `WGCNA`. However, it is best to restrict the data set to only the transcripts (of features) that have sufficient expression levels. We will define "sufficient" as >100 reads across all the samples, but this is not meant to be a universal suggestion and it may not even be "optimal" for this data set. The following code will subset the data accordingly and assumes that the tables are organized the same.

	transformed.counts <- read.table("transformed.counts", header=T, row.names=1)
	all.counts <- read.table("~/Workshop/DESeq2/all.counts", header=T, row.names=1)

	counts.high <- transformed.counts[ which( rowSums(all.counts) > 100) , ]
	dim(counts.high)
	head(counts.high)
	
Let's also reduce the amount of data to a random sample of 5000 genes to reduce computational time.

	counts.high.subset <- counts.high[sample(nrow(counts.high), 5000), ]

Finally, the count table needs to be transposed for use in `WGCNA`.

	counts.wgcna = as.data.frame(t(counts.high.subset))
	dim(counts.wgcna)

We also need the sample information, and this time, we will use the DroughtDummy variable, spatial, and environmental variables.

	sample.info <- read.table("../sample.info", header=T, row.names=1)
	sample.info

	sample.wgcna <- sample.info[ , 5:13]
	sample.wgcna

### Exploring the structure of samples 

Before we begin the core `WGCNA`, we may want to explore and graphically summarize our data at the level of samples. That is to say, we can cluster *samples* (not genes) by gene expression patterns and then examine what factors are related to the observed clustering. To perform the clustering, type

	sample.tree = hclust(dist(counts.wgcna), method = "average")

Next, we can create a color scale for each of the environmental variables in `sample.wgcna`, where white means low, red means high, and grey means missing data.

	env.colors = numbers2colors(sample.wgcna, signed = FALSE);

Finally, plot the sample dendrogram with the colors underneath.
	
	pdf("Dendrogram_bySample.pdf")
	plotDendroAndColors(sample.tree, env.colors, groupLabels = names(sample.wgcna), main = "Dendrogram of Samples")
	dev.off()

Notice that the samples collected before drought are all on one branch of the dendrogram and the after-drought samples are on the other main branch. It looks like the experimental treatment drives large differences in overall gene expression patterns - an observation that might be obvious given the results of `DESeq2` yesterday. Below we can also see how some of the environmental variables might be correlated. In relation to the dendrogram, the correlations do not show an obvious pattern but they may relate to some of the substructure.

### Defining gene expression modules

Now we will define gene expression modules based on hierarchical clustering of the *genes* (not samples) in a very particular way. Co-expression networks are "undirected, weighted gene networks. The nodes of such a network correspond to gene expression profiles, and edges between genes are determined by the pairwise correlations between gene expressions. By raising the absolute value of the correlation to a power [>1], the weighted gene co-expression network construction emphasizes high correlations at the expense of low correlations" [(Langfelder & Horvath 2008)](https://doi.org/10.1186/1471-2105-9-559). The [soft  thresholding  power](https://doi.org/10.2202/1544-6115.1128) is a part of what makes the co-expression network weighted, hence **W**GCNA. Thus, we first need to define this value. This weighted network is called the "adjacency matrix". The adjacency matrix is then converted to another similarity matrix called the "topological overlap matrix", which helps minimize noise and spurious correlations. This, in turn, is then converted to a dissimilarity measure for use in hierarchical clustering. The resulting clustering dendrogram provides a framework for defining modules with various criteria to define minimum module size (number of genes) and cutoff points on the dendrogram. More information on these steps below. Let's start with choosing the soft threshold power.

`WGCNA` includes a function to choose the "best" value for the soft threshold power based on the criterion of approximate "scale-free topology." Scale-free topology is an observed phenomenon in biological networks in which the probability that a node is connected with another node decays as a power law. The function evaluates a series of user-defined powers and returns a set of network indices to help choose the best value.

	powers = c(c(1:10), seq(from = 12, to=20, by=2))
	sft = pickSoftThreshold(counts.wgcna, powerVector = powers, verbose = 5)

	pdf("SoftThresholdingPower.pdf")
	par(mfrow = c(1,2))
	plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence"))
	text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 	labels=powers,cex=0.9,col="red")
	#abline(h=0.90,col="red") 
	plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
	text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.9,col="red")
	dev.off()

The resulting plot on the left shows how the fit to a "scale-free topology" model increases as the chosen power increases. We want to choose the value at the start of the plateau (I chose 9 in my test run).

#### One-step network construction and module detection

Once we choose the soft thresholding power, the entire network construction and module detection process can be performed in one step with the `blockwiseModules` function. Do not run this code because we will run everything step by step below.

	#Do not run
	network <- blockwiseModules(counts.wgcna, power = 9, TOMType = "unsigned", minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = TRUE, saveTOMFileBase = "AllTOM", verbose = 3, maxBlockSize = 5000, corType="bicor")

The function calculates the adjacency matrix using a kind of correlation called biweight midcorrelation (`bicor`, instead of Pearson) and the power we chose above, converts to topological overlap matrix (TOM), performs hierarchical clustering on dissimilarities based on the TOM, and finally identifies modules with a minimum size (number of genes) of 30 by cutting the dendrogram following the `mergeCutHeight` parameter. In the next section, we do this step by step to see the products and consider the various options that you can fine-tune in your analyses.

#### Step by step network construction and module detection

The first step is to calculate the adjacency matrix with the soft thresholding power we chose above and midweight bicorrelation (as suggested by the developers). You could use other correlation metrics, such as Pearson or Spearman, if you prefer.

	adjacency = adjacency(counts.wgcna, power = 9, corFnc = "bicor")

If you look at a piece (first five rows and five columns) of the adjacency matrix you can see it is essentially a symmetrical matrix of pairwise correlations.

	head(adjacency[1:5, 1:5])	

Then, we transform the adjacency matrix into a topological overlap matrix, which is more robust to noise. The idea is to consider the direct connection strengths along with the connection strengths via shared neighbors. In turn, we calculate dissimilarity from the TOM.

	TOM = TOMsimilarity(adjacency)
	head(TOM[1:5, 1:5])
	
	dissTOM = 1-TOM
	head(dissTOM[1:5, 1:5])

From the TOM dissimilarity, we can perform hierarchical clustering.

	gene.tree = hclust(as.dist(dissTOM), method = "average")

Plot the resulting clustering tree (dendrogram).
	
	pdf("Dendrogram_byGene.pdf")
	plot(gene.tree, xlab="", sub="", main = "Gene clustering on TOM-based 	dissimilarity", labels = FALSE, hang = 0.04);
	dev.off()

In the dendrogram, each of the tips is a gene. There are so many that they are hard to distinguish individual lines. To identify how many modules there are we will use the function `cutreeDynamic`. We have to define a minimum module size, which I arbitrarily set to 20.

	dynamicMods = cutreeDynamic(dendro = gene.tree, distM = dissTOM, deepSplit = 4, pamRespectsDendro = FALSE, minClusterSize = 30);
	table(dynamicMods)

To plot, we convert numeric labels into colors and then use a built-in plotting function

	dynamicColors = labels2colors(dynamicMods)
	table(dynamicColors)
	
	pdf("Dendrogram_byGene_withModules.pdf")
	plotDendroAndColors(gene.tree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
	dev.off()

Beneath the dendrogram we now see colors depicting which module each gene belongs to (gray represents unassigned genes). How many modules did you get? 

`cutreeDynamic` sometimes leads to "too many" modules, meaning that some modules may still be highly correlated with each other. We might want to combine those. To do so, we calculate the module "eigengene", which is a synthetic gene that captures the main axis of variation of all genes in the module. As the name implies, the eigengene is generated based on the first axis of principal components analysis. To calculate, we can use the built-in function.

	MEList = moduleEigengenes(counts.wgcna, colors = dynamicColors)
	MEs = MEList$eigengenes
	MEs

Then, we can build a dendrogram of the eigengenes and use it to combine modules based on a cutoff that we might choose. That means we need a dissimilarity matrix based on module eigengenes to then cluster.

	MEDiss = 1-cor(MEs)
	METree = hclust(as.dist(MEDiss), method = "average");
	
	#Plot the resulting eigengene dendrogram
	pdf("Dendrogram_byEigengene.pdf")
	plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
	abline(h=0.25, col = "red")
	dev.off()

In my dendrogram, several modules are defined by low heights less than 0.25, so we might consider merging these.
	
	#Choose a cutoff
	MEDissThres = 0.25

	#Call an automatic merging function
	merge = mergeCloseModules(counts.wgcna, dynamicColors, cutHeight = MEDissThres, verbose = 3)

	#The merged module colors
	mergedColors = merge$colors

	#Eigengenes of the new merged modules:
	mergedMEs = merge$newMEs

And, finally we can make a new dendrogram with colored membership bar based on the new merged modules to replace our old one.

	pdf("Dendrogram_byGene_withMergedModules.pdf")
	plotDendroAndColors(gene.tree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
	dev.off()

	#Then, reassign these colors and labels to the variables we were using before
	moduleColors = mergedColors
	colorOrder = c("grey", standardColors(50))
	moduleLabels = match(moduleColors, colorOrder)-1
	MEs = mergedMEs

We may also want to know what genes are in each module, and which gene within each module is most connected (*i.e.*, hub genes).

	#All genes in "blue" module. Modify color accordingly.
	gene.names = names(counts.wgcna)
	modGenes.blue = gene.names[moduleColors=="blue"]
	length(modGenes.blue)
	modGenes.blue
	
	#Hub genes
	chooseTopHubInEachModule(counts.wgcna, moduleColors)

The list of genes in each module will be useful for further exploration of functional annotation and enrichment tomorrow. The specific hub genes might also have interesting functions that interact with many other genes. Which genes are hub genes in your analysis? Look up their functions in the annotation files I included in the reference transcriptome folder. Anything interesting?

So far what we have learned is that there are several groups of genes that show similar expression patterns as represented by thier respective modules. Modules might make it easier to narrow down interesting genes for further analysis. Nonetheless, there is more we can do to see what might explain (or be correlated with) the co-expressed modules.

### Network visualization

Before we move on to testing what might explain module membership, let's briefly explore network visualization in the way most people have in mind, with nodes and edges connecting nodes. With the list of genes in a particular module, we can extract the relevant values from the topological overlap matrix. 
	
	inModule.blue = moduleColors=="blue"
	modTOM.blue = TOM[inModule.blue, inModule.blue]
	dim(modTOM.blue)
	head(modTOM.blue[1:5,1:5])
	dimnames(modTOM.blue) = list(modGenes.blue, modGenes.blue) #modGenes.blue defined above
	
To produce a clearer plot, we can subsample the top hub genes.

	nTop = 30
	genes = names(counts.wgcna)
	IMConn = softConnectivity(counts.wgcna[, genes[(moduleColors=="blue")]]);
	top = (rank(-IMConn) <= nTop)
	modTOM.blue[top, top]
	
To graph the network based on the TOM, we need to use other software. One option is to export for other software (more below). Another is to use an `R` package such as `igraph`. 
	
	library(igraph)
	graph <- graph.adjacency(modTOM.blue[top, top], weighted=TRUE, mode="undirected", diag=FALSE)
	
	pdf("Network.blue.pdf")
	plot(graph)
	dev.off()
	
The plot shows a big blob in this case but you might want to play around with setting threshold for showing connections, or weight the thickness of lines accorging to strength, or adjust the node size, or choose another module that might be clearer or more interesting, *etc*.

	#e.g., setting TO < 0.1 to 0
	modTOM.blue[ abs(modTOM.blue) < 0.1] <- 0
	
	#Then rerun graph.adjacency and plot
	
See `?graph.adjacency` and `plot` for more details on how to modify the visualization to suit your needs.

There is advanced external software for drawing and anlyzing networks, which you may also want to consider. `WGCNA` provides [functions and code](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-06-ExportNetwork.pdf) to output data for two such programs, [VisANT](http://visant.bu.edu/tutorials/index.html) and [Cytoscape](https://cytoscape.org/). We do not have an annotation file in the right format to create the output, but I encourage you to do so with your own data and explore other network visualizations.

### Module correlations with external variables

We might wonder whether certain measurements (*e.g.*, traits, environmental variables, *etc*.) taken on the samples correlate with each module. For example, in our case we might expect genes responding to the drought treatment to show similar expression patterns and thus form a module. We might also expect that the environment of the seed sources will relate to gene expression, for example, if they are locally adapted to their home environments. 

First, we define a few terms and reorganize a bit.

	nGenes = ncol(counts.wgcna)
	nSamples = nrow(counts.wgcna)
	
	#Recalculate module eigengenes with color labels
	MEs0 = moduleEigengenes(counts.wgcna, moduleColors)$eigengenes
	MEs = orderMEs(MEs0)

Then, we calculate the correlation between the environmental data (in `sample.wgcna`) and our module eigengenes and generate *p*-values.

	moduleTraitCor = cor(MEs, sample.wgcna, use = "p")
	moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

Finally, we can plot the results in heatmap form, showing the correlations and *p*-values.

	textMatrix = paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "")
	dim(textMatrix) = dim(moduleTraitCor)

	pdf("Heatmap_Module-EnvironmentCorrelations.pdf")
	par(mar = c(6, 8.5, 3, 3))
	labeledHeatmap(Matrix = moduleTraitCor,xLabels = names(sample.wgcna), yLabels = names(MEs), ySymbols = names(MEs), colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = textMatrix,setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1), main = paste("Module-trait relationships"))
	dev.off()

Which variables are most strongly correlated with each module? In mine, one or two modules are very highly correlated with the drought treatment, which is perhaps no surprise. Several modules have moderate to high correlations with variables related to the source environment. What might this mean?

We can look in a bit more detail at the intramodular correlations. I will choose the `DroughtDummy` since it is so highly correlated with one of the modules, but you should take a look at some of the others too. The details of the following code are not so important, but essentially we are looking at the correlation between module membership of a gene and the significance of that gene's association with drought, in this case. Gene signficance is the correlation between the gene and the environmental variable, while module membership is the correlation of the module eigengene with the gene expression profile. 

	drought = as.data.frame(sample.wgcna$DroughtDummy)
	names(drought) = "Drought"

	modNames = substring(names(MEs), 3)
	geneModuleMembership = as.data.frame(cor(counts.wgcna, MEs, use = "p"))
	MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
	
	names(geneModuleMembership) = paste("MM", modNames, sep="")
	names(MMPvalue) = paste("p.MM", modNames, sep="")
	
	geneTraitSignificance = as.data.frame(cor(counts.wgcna, drought, use = "p"))
	GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
	
	names(geneTraitSignificance) = paste("GS.", names(drought), sep="")
	names(GSPvalue) = paste("p.GS.", names(drought), sep="")

	module = "turquoise"
	column = match(module, modNames)
	moduleGenes = moduleColors==module
	
	pdf("Drought_MMvGS.pdf")
	par(mfrow = c(1,1));
	verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]), abs(geneTraitSignificance[moduleGenes, 1]), xlab = paste("Module Membership in", module, "module"), ylab = "Gene significance for drought", main = paste("Module membership vs. gene significance\n"), cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
	dev.off()           
                   
Clearly the drought experiment is a major predictor in this module. It is interesting how many genes are in that module and thus respond with the treatment. The effect seems to dominate any other patterns. If we wanted to dig deeper into background patterns of co-expression unrelated to the drought treatment, we might consider reanalyzing the data with only the pre-drought (before) samples.      

Tomorrow, we will explore how to test for functional enrichment, within modules from today or among differentially expressed genes from yesterday. 
