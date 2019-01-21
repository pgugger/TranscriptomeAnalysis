## EXERCISE 1: Preparing Illumina sequence data files

### Viewing sequence files

Log in to the cluster with `ssh` using the Terminal on a Mac or PuTTY on Windows. Type the following and enter your password

	ssh -p 60307 username@10.99.0.146
	
The workshop folder has been copied to your home folder already. View the contents of you home folder by typing `ls`, then navigate to the `Workshop` folder with `cd Workshop` and type `ls` again to view the contents.

	
Navigate to the folder called "RNA-Seq_Data", and list the files in the folder (`cd` and `ls`). When you get your Illumina data, they will likely come in FASTQ files, either in a series of files related to how they came off the sequencer, or as a set of files that have been demultiplexed, meaning each file represents one sample. If you used standard read indexing in the adapter, the sequencing center will usually demultiplex for you. Here, you will see that I have given you demultiplexed, compressed (.gz) FASTQ files, each named after a sample except one called `test.fq.gz`, which is what we will use for this exercise. FASTQ files have four lines per read: header starting with "@", sequence read, "+" (sometimes followed by additional info), and finally quality scores for each base.

We can look at a compressed file without extracting it:

	zcat test.fq.gz | less -S

`cat` is a command that will print text files to screen, while `zcat` does the same for compressed (.gz) files. The `|` then indicates to send, or "pipe", the results of `zcat` to `less` for more controlled viewing. The `-S` flag indicates not to wrap lines when longer than the screen is wide. If you just used `zcat` without `less`, millions of lines of Illumina data would start flashing across the screen uncontrollably. `less` allows you to control how you scroll and view the file. If you scroll down, you can see the file is very long. Press `q` to quit.

We can find out how many lines the file has using `wc`, which counts the number of lines, words, and characters by default. `wc -l` specifically counts the number of lines.

	zcat test.fq.gz | wc -l

How many reads are in this file?

### Assessing quality with FastQC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc) is a convenient program for summarizing base call quality scores, GC content, possible adapter contamination, and other information. It's easy to run but can be slow, so just try it on one file.

	fastqc test.fq.gz

When the program finishes, download the results files to your laptop using an FTP client (*e.g.*, FileZilla). Then open the HTML file by double-clicking. How do the Phred scores look at each base?

### Assessing possible contaminants with FastQ_Screen

Let's say you want to know how much contamination is in your data. You might try to search (`grep`) the FASTQ for part of the adapter sequence and count (`-c`) how many sequences have it. But, remember we need to start with `zcat` if the file is gz-compressed.

	zcat test.fq.gz | grep -c 'GATCGGAAGAGCGGTTCAGCAGG' 

However, what if there are errors, so you don't have an exact match to the adapter, and what if there are other contaminants like your own DNA or another organism? One way to get a rough estimate is to use [FastQ_Screen](http://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqscreen), which will align your reads to specified genomes or sequences).

First, we need to set up the configuration file, which also represents an opportunity to learn how to edit text files in the terminal. Open fastq_screen.conf with `vi` (or using your FTP client (*e.g.*, FileZilla) if you prefer not to learn the command line way). The top section indicates the path to the [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) aligner for aligning short reads to genomes. The "Threads" section indicates how many processors your computer has available to run the analysis. The last section, called "Database", give the path to different genomes or sequences that you want to align your data to. Note that the genomes/sequences must be indexed first (not FASTA) using `bowtie2-build`, which I have already done for you and placed in `~/Workshop/Genomes/`. Edit the path names to reflect your working directory (*e.g.*, in the "Database" section, change `/home/pgugger/*` to `/home/your_username/*`). To edit in `vi`, type `i` for "insert" then use the keyboard to edit; press `ESC` to stop editing; and then save and close the file by typing `:x` (or close without saving using `:q!`). You can now run the program on one example file.

	fastq_screen --aligner bowtie2 --conf ~/Workshop/RNA-Seq_Data/fastq_screen.conf test.fq.gz

When the program finishes, download the results files to your laptop, and open the HTML file. What fraction of reads align to the specified genomes? How many don't align anywhere? What could that mean!?

`fastq_screen` can also be used to filter reads that map to undesired organisms (see [documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/fastq_screen_documentation.html)), but we will ignore this for the workshop. Furthermore, it is unlikely that reads from non-target species will map to the reference sequences that we will use tomorrow. 
	
### Quality filtering and trimming reads

For many applications, you may want to remove low quality reads or adapter contamination, especially when they are found to be substantial in the analyses above. There are many software packages available to do this. Here, we will use [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) because it has many useful options for trimming and filtering. 

All of the files are already trimmed and filtered, except that I have provide you with one test file. Try running Trimmomatic on `test.fq.gz`.

	java -jar path/to/trimmomatic.jar SE -phred64 test.fq.gz test.trim.fq.gz TOPHRED33 ILLUMINACLIP:TruSeq3-SE:2:30:10 SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:40

While it is running read about the various options I chose. These settings follow the suggested starting point but can be adjusted as needed. This command clips adapter contamination (`ILLUMINACLIP`) following the specified rules, removes low quality bases using a sliding window approach (`SLIDINGWINDOW`), trims bases off the start and end if they are very poor (`LEADING` and `TRAILING`), and discards any read that winds up under 40 bp (`MINLEN`). Notice that I also had it convert from the old Phred64 quality encoding to the newer Phred33.

When you have finished this exercise, either delete (`rm`) `test.fq.gz` and `test.trim.fq.gz` or create a new directory (`mkdir test`) and move them there (`mv test.* test/`) so they are not in the way for tomorrow's exercise.

In Exercise 2, you will use the filtered real data to learn how to align sequences to a reference, generate read counts per transcript, and then test for differential expression among experimental treatments.
