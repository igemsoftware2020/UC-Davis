Tutorial:
===============

In this tutorial you will learn how to find the optimal parameters of the motif
 discovery tool MEME within the MEME Suite 5.1.1 software distribution. This
software was developed as a wrapper for the desktop download version of MEME in
order to automate the continuous running of the MEME tool with a variety of 
different parameters included within the meme software. 


This tutorial includes instructions for:
1. Setting up
2. Generating .fa files with mbed
	- Testing mbed
	- Testing markov analysis of mbed
3. Using memescape to find explore MEME capabilities
	- Testing memescape
	- Using markov analysis function of memescape
4. Understanding program outputs and finding optimal parameters
	- Notes on scoring method
	- Condensed data display
		- Notes on condensed data display
	- Uncondensed data display
		- Notes on uncondensed data display 
5. Finding binding sites using motifomatic

## Disclaimers
We make several assumptions about users of this software. If any of the 
following expectations are not met, users may discover challenges with using
 this software. 

+ It is assumed that the user is familiar with Unix/Linux, is comfortable 
navigating directories and running software programs from the command line.
+ It is assumed that the user knows how to create environment variables and 
edit their PATH and PYTHONPATH.
+ It is assumed that MEME suite version 5.1.1 is installed locally and can be 
run as an environment variable (see MEME suite software page for further 
information.
+ It is assumed that user possesses basic understanding of the MEME suite motif
 finder MEME and its search models (i.e. zoops, oops, anr)
+ It is assumed that the user understands the jaspar and fasta file formats. 
+ It is assumed the user has basic understanding of position weight matrices 
and Markov analysis. 
+ It is assumed that the user has installed Python>3.6
 
 
### 1. Setting up:

1. Before installing the UC Davis software package, it is important to first 
install the MEME suite and to make sure that it is working locally. 
2. Download the UC Davis repository and update `PATH` to include 
`UC-Davis/bin`.
	+ For example, if `UC-Davis` was downloaded to your Home directory the 
	PATH should now include `${HOME}/UC-Davis/bin`
		+ If UC-Davis is downloaded to a different directory the PATH
		MUST have the correct filepath
	+ Click [here](https://astrobiomike.github.io/unix/modifying_your_path)
	for more information on writing a PATH 
3. Update `PYTHONPATH` to include `UC-Davis/` 
	+ For example, if UC-Davis is downloaded to the home directory, 
	PYTHONPATH should now include `${HOME}/scratch/UC-Davis/`
	+ If not done so already, make sure both PATH and PYTHONPATH are 
	exported in the PATH file
4. After updating the `PATH` and `PYTHONPATH`, open a new terminal and check
 that the programs in bin are working as environment variables.
	+ Type a program name in the terminal. If UC-Davis/bin and
	 UC-Davis/sporecore are  set up correctly, the command line prompts 
	 for the program should appear.
5. Navigate to the UC-Davis directory and run the following command:
`python3 setup.py test`
	+ If all tests pass, you may proceed. 
	+ If errors, contact lilymaryam5@gmail.com
	
**Note that if a user cannot or wishes not to update their environment
variables as described above, programs in this software package can be run 
by navigating to the UC-Davis directory and typing `python3 bin/PROGRAM-NAME`.
If the user prefers this option, for all test commands written in this TUTORIAL
they will type the above command instead of only typing the program name.
	+However, if the user chooses to do this, they are still expected to 
	install MEME as directed (which requires updating the PATH with MEME 
	suite)**

 
### 2. Generating .fa files with mbed
mbed is a python program that will create a FASTA formatted file of synthetic 
promoter sequences. The program accepts several promoter formatting parameters
as well as a jaspar-formatted file containing the position weight matrix of a 
transcription factor binding site of interest. Promoter sequences are generated
and a motif binding sequence representing the the position weight matrix of the
binding motif is embedded in the promoter sequence at a known position. 
 
Promoter background sequences can be generated in one of two ways. Promoter 
sequence can be generated from the markov analysis of a background DNA 
sequence fasta file or with a user-specified nucleotide distribution frequency.
After promoter background sequence is created, a motif binding sequence is 
embedded and recorded at a known position.  


These promoters are inserted into a fasta file which is intended to be run 
through the motif finding software MEME to test its ability to find binding 
sites. Using a known motif and motif location allows for software performance
to be precisely analyzed and compared across different parameters.


To run mbed, there is one required command line argument: 
	
	`--jasparfile`


This argument specifies the file path to a jaspar file 
containing the position weight matrix of the binding motif to be embedded in 
the promoters. For this program, the user must have jaspar files available or 
format their intended position weight matrix in jaspar format. A testfile 
'testjaspar.1.jaspar' is provided for testing mbed, as well as other 
relevant programs mentioned in this tutorial.

To use the markov analysis function of mbed, the user must provide a fasta
formatted file containingDNA sequence. The purpose of the markov analysis
 option for mbed is to allow users to create promoters with higher complexity
 sequence. The optional command line arguments for markov analysis in 
mbed are as follows.

+ `--dnafile`: a string argument that specifies the file path to a fasta file 
that contains genomic or promoter sequence from actual species 
+ `--markov_order`: an integer argument (default 0) that indicates the markov 
order to be used for analyzing the DNA sequence
	+ note that if markov order specified exceeds complexity of given sequence, 
  an error will occur


  
To initiate markov analysis, `--dnafile` must be input with a valid file in
the command line. The default is a markov order of 0. User can specify higher 
markov orders which are greater than 0.

When `--dnafile` isn't input, the program will generate sequence using a 
nucleotide distribution frequency that can be modified in the command line.


mbed also has several optional arguments with adjustable default settings that 
can modify the final fasta output file. 

+ `--numseq`: an integer argument (default 10) that specifies how many promoter
sequences will be generated
+ `--seqlen`: an integer argument (default 100) that specifies how long each 
promoter in the fasta output file will be 
	+ note that --seqlen can't be shorter than the length of the binding motif
+ `--mps`: an integer argument (default 1) that determines how many motif 
sequences will be embedded in per promoter
+ `--freq`: a float argument (default 0.9) that specifies a probability at 
which motifs will be embedded in promoters (i.e. if --freq = 1.0, all promoters 
will have a motif, if --freq = 0.9, on average 9/10 promoters will have a 
motif)
	+ note that --freq must be a number between 0 and 1.0
+ `--PA`, `--PC`, `--PG`, `--PT`: four separate float arguments (each default 
0.25) that specify each nucleotide's frequency in the promoter background 
sequence
	+ note that the sum of these 4 values has to add up to 1.0
+ `--bothstrands`, `--negstrands`: two different switch options (default is 
positive strands only) for the strandedness of the generated promoters
	+ note that only one type of strandedness may be chosen at a time and 
   `--bothstrands` is the advised setting for more realistic results

  
  
To view these options and defaults again in command line, run the command:
 `mbed -h` 
 
#### Testing mbed

To test mbed it is recommended that the user try to create a fasta file with a
specific set of parameters. In this example it is assumed that the user is 
working in the UC-Davis directory (which was downloaded from the GitHub 
repository) and any paths are relative to that.

###### Example
If a user wanted to make 15 promoters that are each 500 base pairs long and 
embed a binding motif from testjaspar.1.jaspar (the test file mentioned above)
 they would need to first find the jaspar filepath 
('data/testjaspar.1.jaspar'). They could further specify that they wanted ~0.8
 of the total promoters to have a single embedded motif and that the background
A,C,G,T frequencies of the promoter files should be  0.3, 0.2, 0.2, and 0.3 
respectively. 

For this scenario, the command would be: 
 
`mbed --jasparfile data/testjaspar.1.jaspar --numseq 15 --seqlen 500 --mps 1
 --freq 0.8 --PA 0.3 --PC 0.2 --PG 0.2 --PT 0.3 `
 
This command will print the promoter sequences directly to the
Terminal in FASTA format. To create a file for future use, it is essential to 
select a directory and file name to save your file to. 

For example: '~/sporecorefiles/testpromoters123.fa'

To create this file, the command would be:

`mbed --jasparfile data/testjaspar.1.jaspar --numseq 15 --seqlen 500 --mps 1 
--freq 0.8 --PA 0.3 --PC 0.2 --PG 0.2 --PT 0.3 >  
~/sporecorefiles/testpromoters123.fa`

#### Testing markov analysis of mbed

To run the same scenario using markov analysis, the user would need a 
background DNA sequence. As a testfile, the complete sequence of 
 *Saccharomyces cerevisiae S288C* Chromosome I was downloaded from GenBank in
 FASTA format. This testfile is called “sequence.fasta” and can be found via 
 this filepath: 'data/sequence.fasta'.

Next the user must decide which markov order they want to use to train mbed. 
For the provided testfile sequence.fasta, the maximum order possible is 6. 
Selecting a higher number will result in an error message. To increase markov
 order one must also increase the amount (length) of DNA sequence in the file. 

It is again assumed that the user is working in the UC-Davis directory and any
paths are relative to that. To effectively run this program, a competent 
understanding of Markov models is recommended, however, a user can run this 
option if they simply know that higher complexity sequence will create more 
realistic promoters. 

###### Scenario:
If a user wanted to make 15 promoters 500 base pairs long and embed a binding
motif from testjaspar.1.jaspar (the test file mentioned above) they would need
the jaspar filepath: data/testjaspar.1.jaspar. They could further specify that 
they wanted ~80% of the total promoters to have a single embedded motif.

In this scenario they would like a more complex promoter sequence that more 
closely imitates the sequence of real promoters. For this example they must 
train the program with real sequence data. 

(Note that sequence data ideally should be extracted from promoter regions 
upstream of genes, however any sequence data is accepted.)

For this scenario provided, the command would be as follows: 
 
`mbed --jasparfile data/testjaspar.1.jaspar --dnafile data/seqeunce.fasta 
--markov_order 6 --numseq 15 --seqlen 500 --mps 1 --freq 0.8`
 
This command will print the promoter sequences directly to the terminal in 
fasta format. To use this file in future experiments, it is essential to select
a directory and file name to save your file to. 

For example: '~/sporecorefiles/testpromoters1234.fa'

To create a file, the command would be:

`mbed.py --jasparfile data/testjaspar.1.jaspar --dnafile data/seqeunce.fasta 
--markov_order 6 --numseq 15 --seqlen 500 --mps 1 --freq 0.8 > 
 ~/sporecorefiles/testpromoters1234.fa`

**A note about mbed:**
Though mbed has two slightly different promoter creation approaches, all 
outputs are identically formatted and thus files are indistinguishable unless
named accordingly.

(For file output format see fasta_output.fa)

At this point a fasta file containing synthetic promoters with clearly 
identified biding sites can be used for further experimentation. Motif
position and strand in a given promoter file is stored in the sequence name 
like so:

">seq-1 ['77 +']'" - indicating that the motif is at base pair position 77 on 
the positive strand


### 3. Using memescape to explore meme capabilities 

Although mbed is a standalone program, it was originally written with the 
intention of creating files to test the performance of the MEME suite 
software's motif-finder (called MEME) tool. By creating fasta 
files with exact motif location, strand, and position weight matrix known,
 MEME's ability to correctly locate motifs in a promoter sequence can be 
quantified and optimal running parameters can be designed. 

To execute this quantification and optimization, memescape was created. 
memescape is a wrapper program that automates the creation of synthetic 
promoter files (using mbed) and automates the analysis of those files by 
running MEME on them.

It is the intention of this program to run promoter files of various 
specifications through MEME with different motif-finding parameters in order to
determine when meme performs best. For each iteration in memescape, a fasta 
file with set parameters is created and run through MEME once with each model 
(oops,zoops, anr). For each run, memescape parses and records the fasta 
promoter file name, the file specifications for promoter number and size, the 
position and strand of the embedded motif in each sequences, and the position 
weight matrix of the embedded motif itself. 

After the MEME analyzes the file, memescape opens the MEME output and
extracts relevant information to compare to compare the expected motif 
information (from the initial FASTA file). After comparing and analysing the 
extracted information, memescape parses and prints all relevant performance 
data which can be saved to a .csv file for viewing. 

Using memescape for performance analysis is highly recommended over running 
meme manually multiple times for a single file because it allows for hundreds
of possible parameter combinations to be compared in a relatively short and 
unsupervised time frame.

memescape has the ability to find optimal parameters for motif finding as well
as the ability to definitively determine which scenarios will never 
successfully reveal binding motifs. In this way scientists can discover how 
the number and length of promoters that are needed and can attempt to adjust 
their data accordingly. All performance data generated creates guidelines that 
can increase confidence in results and reduce time spent running individual 
experimental datasets.  


memescape has 2 required arguments:
+ `--jasparfile`: this argument is exactly the same as the Mbed and markov_Mbed
argument(see above in #2) 
+ `--memepath`: a string argument that indicates the filepath of the MEME motif 
finder
	+ note that desktop version of MEME suite must be downloaded locally
	+ within the MEME suite directory path to meme finder tool is /bin/meme 

memescape has several optional arguments pertaining to the formatting of the 
fasta files. As with mbed, promoter sequence can either be generated with 
nucleotide distribution frequency or markov analysis of background DNA. There
are two arguments (the same as for mbed) associated with this decision.

+ `--dnafile`: a string argument that specifies the file path to a fasta file 
that contains genomic or promoter sequence from actual species 
+ `--markov_order`: an integer argument (default 0) that indicates the markov 
order to be used for analyzing the DNA sequence
	+ note that if markov order specified exceeds complexity of given sequence, 
	an error will occur

As in mbed, inputting a fasta-formatted DNA file will result in memescape using
markov analysis to create promoter sequences. When `--dnafile` is not input in 
the command line, a simple nucleotide distribution frequency is used for 
promoter generation. 

Because multiple promoter files are being tested, memescape expects the 
user to input a range of promoter lengths through three optional arguments. 
 
+ `--minpromoterlength`: an integer argument (default 100) that specifies the
shortest promoter length to be tested
	+ must be greater than or equal to the length of the binding motif
+ `--maxpromoterlength`: an integer argument (default 400) that specifies the 
longest promoter length to be tested
	+ must be greater than or equal to the shortest promoter length
+ `--promoterlengthstep`: an integer argument (default 100) indicating the step
value for the range of promoters 
	+ note that smaller step will increase program run time
	+ note that if previous promoter length in the range + promoterlengthstep > 
	maxpromoterlength, maxpromoterlength will not be tested

If only a single promoter length is desired for testing, maximum and minimum
promoter length can be set equal. This will shorten program runtime.

A similar setup is employed for testing a range of promoter numbers. 

+ `--minnumseq`: an integer argument (default 5) indicating the minimum number 
of promoters in a file 
	+ note that minimum must be greater than 1
+ `--maxnumseq`: an integer argument (default 20) indicating the maximum number
of promoters in a file
	+ note that maximum must be greater than or equal to minimum numseq
+ `--numseqstep`: an integer argument (default 5) indicating the step value for
the range of promoter numbers in a file
	+ note that smaller step will increase program run time
	+ note that if “previous promoter number in the range” + “numseqstep” > 
	maxnumseq -- maxnumseq will not be tested
  
These arguments are intended to give user control over how many and which 
promoter lengths they would like to test. It is recommended to create a 
realistic range for each variable (i.e. promoters are usually 400-1000bps 
long), and to keep step values high to reduce parameter redundancy and runtime. 

memescape also has five optional arguments that specify promoter file 
conditions: 

+ `--motiffrequency`: corresponds to `--freq` argument of mbed (see #2 above)
+ `--PA`, `--PC`, `--PG`, `--PT`: corresponds to arguments in mbed (see #2)
	+ when memescape isn't using markov-analysis for promoter generation,
	these arguments determine nucleotide frequency in the promoter background
	sequence

The defaults for the above arguments are set equal to the defaults in mbed.
When these arguments are not called, user is electing to use mbed defaults. 



memescape has 2 optional arguments that control the output display:

+ `--condenseddata`: a switch argument that when True creates an output that 
scores the performance of binding motifs found by meme. It finds a success 
rate by comparing the number of accurate sites found to the total number of
sites found. Each motif performance is based on the combined performance of
each promoter sequence.
	+ Note that the default (when condenseddata is not True) analyzes and 
	displays each sequence in the fasta file separately and how meme results 
	compared to expected results
	+ note that `--condenseddata` is recommended for analyzing optimal 
	parameters
+ `--numiterations`: an integer argument (default 1) that specifies how many 
iterations of a given set of parameters are desired. When numiterations is
greater than 1, multiple promoter files for each combination of paramenters
will be generated. 
	+ when `--condenseddata` is True and numiterations > 1, motif performance 
	data for that given set of parameters is averaged to provide more accurate
    insights for motif-finder performance at those parameters
	+ when `--condenseddata` is False and numiterations > 1, each file 
	generated with a given set of parameters will be assessed as usual and all 
	sequence performance data will be displayed 
       + note that numiterations must be greater than 0
 
There are three optional arguments that may be used to modify meme search 
conditions:

+ `--nummotifs`: an integer argument (default 1) that specifies how many motifs
meme should find
	+ note that the first motif identified by MEME is the best scoring one
	which means that when `--condenseddata` is True and `--numiterations`>1 
	similarly scoring motifs will have their performances averaged with other 
	motifs with the most similar performance (this is best scenario for
	averaging)
+ `--min_w`: an integer argument (default 6) that specifies the shortest motif
length meme should find 
+ `--max_w`:an integer argument (default 15) that specifies the longest motif
length meme should find 

(Note that there is no option for choosing the strandedness of motifs on the 
promoters. This is because by default, all files have motifs on bothstrands and
meme looks for motifs on both strands.)

#### Testing memescape:

With all these options in mind, if a user wanted to test the effectiveness of 
meme on promoter sizes from 200bp, 300bp, 400bp, and 500bp, and on fasta files
containing 10, 15, and 20 promoters they would use the following arguments:

`--minpromoterlength 200 --maxpromoterlength 500 --promoterlengthstep 100
 --maxnumseq 20 --minnumseq 10 --numseqstep 5`

If they wanted to change the motiffrequency in the fasta file to 0.8 of all 
promoter sequences and change the background nucleotide frequencies of A, C, G,
and T to 0.3,0.2,0.2, and 0.3 they would add arguments:
 
`--motiffrequency 0.8 --PA 0.3 --PC 0.2 --PG 0.2 --PT 0.3` 

If they wanted to understand the performance of motifs found by meme rather
than performance of individual promoter sequences, they could add the argument
`--condenseddata`. If they wanted an average performance of the motifs found by
the motif-finder for a given set of memeOP parameters they could increase the 
number of iterations to 3 (for example) with the argument `--numiterations 3`.

If all of these options are chosen for the test jasparfile, the command line 
will look like so:

`memescape --jasparfile data/testjaspar.1.jaspar --memepath meme 
--minpromoterlength 200 --maxpromoterlength 500 --promoterlengthstep 100 
--maxnumseq 20 --minnumseq 10 --numseqstep 5 --motiffrequency 0.8 --PA 0.3 
--PC 0.2 --PG 0.2 --PT 0.3 --numiterations 3 --condenseddata`

As is, this command will print all memescape output to the terminal. To save 
output data into files for further analysis, user needs to choose a 
directorypath and filename ending in .csv:

`memescape --jasparfile data/testjaspar.1.jaspar --memepath meme 
--minpromoterlength 200 --maxpromoterlength 600 --promoterlengthstep 100 
--maxnumseq 20 --minnumseq 10 --numseqstep 5 --motiffrequency 0.8 --PA 0.3
 --PC 0.2 --PG 0.2 --PT 0.3 --condenseddata > {directorypath}/output.csv`
 


#### Using markov analysis function of memescape:

To adapt the above scenario for markov generated promoters, first a DNA 
background file must be chosen. Using the testfile: “sequence.fasta”, memescape 
can run the same parameters with higher complexity promoters.  

The purpose of the markov analysis is to create more realistic promoters that 
challenge MEME to discover binding sites in higher complexity sequence. 
memescape is intended to give a more realistic idea of the boundaries of meme's
ability to find binding sites. Our goal is to find how many promoters are 
needed, how long promoters must be, and what meme parameters work best to find
motifs. 

To adapt the situation described in above, we can use the test dnafile
sequence.fasta provided in the data directory. The path to this file is 
'data/sequence.fasta'. We can use the markov order 6 to obtain high 
complexity sequence. This time we assume that the output file is saved
as .csv to a desired directory like so: 
`> {directorypath}/markovmemeoutput.csv`

The resulting command would be as follows:

`memescape --jasparfile data/testjaspar.1.jaspar --memepath meme --dnafile
 data/sequence.fasta --minpromoterlength 200 --maxpromoterlength 600 
 --promoterlengthstep 100 --maxnumseq 20 --minnumseq 10 --numseqstep 5 
 --motiffrequnecy 0.8 --condenseddata > {directorypath}/markovmemeoutput.csv`
 
The output of memescape is a large, comma-separated table 
with many values to observe and understand. 


### 4. Understanding program outputs and finding optimal parameters

memescape outputs are intended to be saved in comma separated value (.csv) file
format. As such, at the end of the written command, it is encouraged that the 
user write outputs to a file like so: 
`> {intended_directory}/{unique_filename}.csv`

Following program execution, the saved file should be opened in a spreadsheet 
display such as excel or google sheets.

It is important to understand that all fasta files generated in memescape are
saved into the computer in the /tmp directory. Each file was named and saved 
in the following format:

'/tmp/testmotif{processID}\_{promoterlength}\_{promoter-number}_{iteration}.fa'

where processID is a single value associated with the command run, the promoter
length indicates how long each promoter is, the promoter number tells how many 
sequences are in the file, and the iteration indicates which file containing 
the same parameters it is.

The purpose of saving these files was for manual inspection of promoters and
ability to confirm accuracy of the programs. These files are intended to be 
deleted or saved somewhere permanently shortly after use. Periodically deleting
files from /tmp is recommended for data management and storage. Selecting,
moving, and renaming promoter files of interest is also recommended for data
storage. 

Many of the data sections listed below contain information on how each value 
was calculated and what it indicates. Before delving further into all of the 
columns available to analyze, the main performance indicator-score-should be 
explained.

#### Notes on scoring method:
The score for this output describes the similarity between the original 
position weight matrix from the jaspar file and the position weight matrix
extracted from the meme output file. There were several possibilities
for scoring accuracy between these two data structures. The method chosen was 
to aligned the two motifs for any possible overlap. 

Each position in a given alignment was scored for similarity by adapting the 
concept of Manhattan distance. For each position of the aligned position weight 
matrices, the difference between each nucleotide frequency was added, the 
absolute value was taken, and the sum was subtracted from 2 (the maximum 
possible similarity score for a single position). The similarity score of each
position was summed to find the total score. The highest score of all the 
alignments was selected and returned


The maximum possible score for a given alignment is was 2* length of the 
shortest pwm.  
The score divided by the maximum score is returned as the value of the
percentage score.

#### Condensed data display
For both programs, if the `--condenseddata` output is selected, each row will
represent a single motif from a given file. The data itself will be stored in
the following columns:
 
1. promoter\_file: the name of the fasta promoter file which can be found in 
   the /tmp directory 
2. motif: the name of the motif discovered by meme
	 - note that there will be a single row for every motif discovered
3. number\_of\_sites: the number of sites where the specified motif was
   discovered in the fasta file 
4. jaspar\_file\_info: describes the information content (in bits) of the given
   motif from the jaspar file
	- will help determine what the information content of discoverable binding 
   sites are 
5. jaspar\_info\_percentage: describes information content of jaspar motif as a
   percentage of the maximum possible information content
	- this is useful for quickly determining how informative a motif is
6. e\_val: the evalue of the motif found by meme 
	- lower e-value indicates a more significant motif discovery
7. motif\_score: a similarity score between the given jaspar file motif and the 
   discovered meme motif
	- see **Notes on scoring method**
	- higher score indicates better performance
8. motif\_score\_percentage: the motif-scores's percentage of the maximum 
   possible similarity score
	- see **Notes on scoring method**
	- higher percentage indicates better performance 
9. failure\_rate: the number of incorrect sites divided by the number of sites
   found for a given motif
	- lower failure rate indicates better performance
10. success\_rate: the number of sites found that overlap the expected sites 
    divided by the total number of sites
	- high success rate indicates better performance
11. false\_positive\_rate: the number of sites where meme found a false 
    positive (no jaspar site was embedded) divided by the total number of sites
	- low false positive rate indicates better performance
12. false\_negative: the number of embedded sites in the fasta file that were
    not discovered by meme divided by the total number of embedded sites in the
    fasta file
	- low false negative rate indicates better performance
13. promoter\_size: the length of the promoters in the fasta file
14. number\_of\_sequences: the number of promoters in the fasta file
15. model: the model used by meme to find motifs ('oops', 'zoops', or 'anr')
16. markov\_order: The markov order meme uses to analyze the sequences
	- will always be 0 in memeOP, will be dependent on user in markov\_memeOP
17. iterations: The number of iterations used to calculate the average 
    performance of each motif in the condenseddata format 
	- if iterations = 1, the performance is not an average
18. stdev(avg\_eval): if the number of iterations is greater than 1, will 
    calculate the standard deviation of the e value
19. stdev(avg\_score): if the number of iterations is greater than 1, will
    calculate the standard deviation of the similarity score
20. stdev(avg\_pscore): if the number of iterations is greater than 1, will
    calculate the standard deviation of the percentage similarity score
21. stdev(avg\_frate): if the number of iterations is greater than 1, will
    calculate the standard deviation of the failure rate 
22. stdev(avg\_srate): if the number of iterations is greater than 1, will
    calculate the standard deviation of the success rate
23. stdev(avg\_fprate): if the number of iterations is greater than 1, will
    calculate the standard deviation of the false positive rate 
24. stdev(avg_fn): if the number of iterations is greater than 1, will
    calculate the standard deviation of the false negative score.


#### Some notes on the above output

The `--condenseddata` argument is the recommended display output for both 
memescape. This version focuses more on the overall accuracy of 
the motif found and quickly quantifies the success of the experimental results.
When reading the results in this output, the best indicators of performance 
are:

	- score
	- percentage score
	- success and/or failure rate
	- number of sites. 


#still relevant?
*By scanning the column of percentage scores, it is simple to quickly identify 
promising motifs by high percentage scores. Scanning these rows with a high 
percentage score, it can further be determined if a motif is high performing 
if the failure rate is low, and the success rate is high. In a given row, if 
standard deviation results are available, large values indicate less 
reliability in performance and may be discounted or subject to more iterations.
Users may determine their own threshold values for score and success, however
generally higher than 75% score percentage is recommended* 

If a motif is deemed 'high-performing' following a run, the initial 
run parameters must be observed and recorded. Motif 
information content, number of promoters, length of promoters, markov model, 
and MEME model all provide information to the user on how to ensure successful 
motif finding and high confidence data. 

The purpose of determining MEME performance confidence is so MEME results can 
be trusted in experimental situations with no verifiable results. 

#### Uncondensed data display
When the `--condenseddata` argument is not called in command line, the data 
is displayed in the following columns:

1. file: name of the promoter file (stored in the /tmp directory)
2. sequence: name of the promoter sequence (from fasta file) where the site is
   found
3. motif: name of the motif found by meme
4. number\_of\_sites: number of sites found for the given motif
5. meme\_position: postion (from the start of the sequence) where motif is
   found by meme
6. meme\_strand: DNA strand (+ or -) that meme finds motif on
7. meme\_width: width of motif found by meme
8. j\_position: actual position of motif in promoter file
9. j\_strand: actual strand of motif in promoter file
10. j\_width: width of actual motif from jaspar file 
11. p\_value: p\_value of the site found by meme
12. score: similarity score calculated between meme and jaspar motifs
    - see **Notes on scoring method**
13. score\_percentage: percentage score is of maximum possible similarity score
	- see **Notes on scoring method**
14. e\_value: e value of meme motif
15. false\_pos: indicator of false postive by meme 
    - 0 means not false positive
    - 1 means it is a false positive)
16. false\_neg: indicator of false negative by meme (0 means not false positive,
    1 means it is a false positive)
17. positional\_distance: distance in base pairs between meme position and 
    jaspar position 
    - if 98: means that no overlap exists between motifs
    - if 99: means that meme had a false negative
    - if 100: means that meme found a false positive
18. fail: an indicator if the site found failed to meet criteria of successful 
    discovery 
    - 0 = met criteria
    - 1 = fail
19. overlap: a value indicating the number of overlapping positions between the
    expected and observed motif
20. overlap percentage: overlap divided by the lenght of expected motif 
21. promoter length: sequnece length of each promoter in the dasta file
22. number\_of\_sequences: number of promoters in the fasta file 
23. model: MEME model used by the motif finder to dsicover motifs
24. markov_order: markov order used by meme to analyze promoters
25. iterations: number of iterations for each set of file and meme parameters
    - only 1 necessary for uncondensed data output
    
#### Notes on uncondensed display
This method is not recommended for optimal parameter discovery. The feature 
allows for direct comparison of meme output and expected output and as such,
data analysis is quite self explanatory. Using this display helps confirm that 
correct information is recorded and displayed as expected. This method can be 
used for any purpose in which sequence information is more important that 
motif information or for scenarios in which the number of parameter 
combinations is more narrow. Please note that each row in this display 
represents a sequence analysis (not motif like above) and there will be
significantly more rows in this set up. 

### 5. Finding binding sites using motifomatic

The program motifomatic was written to automate the running of meme on real
promoter files. 

The program takes one required argument: `--fastafile` 
`--fastafile` accepts a string that encodes the file pathway for a fasta file of
real promoter sequences. 

It accepts multiple optional arguments:

+ `--memepath`: a string argument that specifies the path to the meme motif 
finder (default is 'meme')
+  `--min_w`: an integer argument that specifies the minimum width for motifs
found by meme (default 6)
+ `--max_w`: an integer argument that specifies the maximum width for motifs
found by meme (default 15)  
+ `--trim`: an integer amount specifying the length of promoter sequence 
trimmed from upstream of the gene start site and used to search for motifs 
(default 200)                      
+ `--mememodel`: string argument specifying model used by motif-finder (default
 zoops)
+ `--markov_order`: integer argument specifying markov model used by 
motif-finder (default 0)
+ `--nummotifs`: integer argument specifying the number of motifs that should
be found by the meme motif-finder (default 1)
+ `--jasparfile`: a string specifying the file path to the expected PWM for
verification experiments (no default)
+ `--show_sites`: switch option that displays each site, position, strand and 
p-value found by meme

This program is intended to be a simple way to compile and display meme results
in the terminal. It displays the position weight matrix, its e-value, the 
number of sites found and other values relevant to appraising a proposed 
binding site. 

In scenarios where the binding site is known, but the user wants to 
confirm its presence in the promoters bioinformatically, motifomatic helps to
quantify the search success. By providing a jaspar file containing the position
 weight matrix of the expected binding site, motifomatic can compare and score
  the similarity of the found site and the expected site. This quantification 
can serve as proof of meme's ability to find binding sites for a given scenario
 or as proof that the promoter regions offered contain the expected DNA motifs.
 
When a jasparfile is supplied, it is important that the user knows that their 
promoter sequences are indeed regulated by the expected transcription factor 
and its corresponding binding site. 

The parameter options for motifomatic allow for the user to select preferred
 parameters that have been proven successful in memescape. motifomatic has 
 specifically been used for finding proof of concept by compiling promoters 
 with an expected binding site and comparing the found PWM to the known position
 weight matrix with a high level of accuracy. With the right parameters, 
 motifomatic can be used to search for unknown binding sites in co-regulated 
 promoters.
 
 
 To test-run motifomatic, a sample file of promoter sequences 3694.fa has been
 provided in the data directory. To search for a motif of width 10-15 base pairs
 In a 300bp region of the promoters provided, the user would select the following 
 parameters:
 
 `--fastafile 3694.fa --min_w 10 --trim 300`
 
 If they knew the expected binding motif for the promoters was provided by the
 testjaspar.1.jaspar file they could add the argument `--jasparfile 
 data/testjaspar.1.jaspar` 

As a whole, the command would be:
`motifomatic --fastafile 3694.fa --min_w 10 --trim 300 --jasparfile 
 data/testjaspar.1.jaspar`
 
 This would print all of the relevant information from meme into the terminal. 
 Saving this information to a textfile would be the following command:
 `motifomatic --fastafile 3694.fa --min_w 10 --trim 300 --jasparfile 
 data/testjaspar.1.jaspar > motifomatic_out.txt`
 
 Though the program itself is simple, it allows for the compilation of numerous
 quantitative measurements about the performance of meme on real promoter 
 sequence. It creates the potential for a future wrapper program to automate 
 the running of meme and the collection of meme outputs into a more permanent 
 location.

