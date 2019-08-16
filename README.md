# **BLENDER**
## About

BLENDER is a companion program to the DISCOVER-Seq assay to identify off-target editing sites from MRE11 ChIP-Seq experiments.  It takes aligned bamfiles from the IP experiment and (optionally) a control bamfile and identifies locations with stacks of reads at putative cutsites. PAM sequences can be provided by the user as well as a guide sequence. BLENDER makes use of the ENCODE ChIP-Seq Blacklists for human and mouse which are lists of regions in the genome that have artifactual large numbers of reads. These lists and the control bam plus PAM sequences and the guide are used to filter out false positives.  BLENDER runs on mouse mm10 and human hg38 genomes (blacklists coordinates are for these genomes). 
BLENDER has a driver script run_blender.sh that takes 5 (optionally 6) argmuments and runs identification of putative hits, filtering of those hits, and then creates an SVG of the aligned hits. Files are stored in the output directory given as a parameter to the bash script.
Alternatively, each step can be run separately. 

![schematic](https://github.com/staciawyman/blender/blob/master/BLENDER_schematic.png)

## Running the bash script:

        sh run_blender.sh <path to reference genome> \
            <path to IP bamfile> \
            <path to control bamfile> \
            <guide sequence> <output directory> ["options"]

This will run blender with option c set to 3 (details below), this means the program will run quickly, but may miss some very sparsely covered off target sites. I recommend running it initiall with c set to 3 (the default) and then running again with c set to 2 (it may take several days, depending on the number of reads in your bamfile). 

## To run the blender script directly:

        perl blender.pl [options] <reference genome> <guide sequence> <IP bamfile> <control bamfile>  > unfiltered_output.txt

        perl blender.pl [options] <reference genome> <guide sequence> <IP bamfile> <control bamfile>  | perl filter.pl > output.txt

        perl blender.pl [options] <reference genome> <guide sequence> <IP bamfile> <control bamfile>  | perl filter_pool.pl > pooled_output.txt

BLENDER can be run with or without being piped through the filtering script. There are two filtering scripts provided; the standard filter.pl script that implements the standard scoring scheme, and the filter_pool.pl script that implements the more stringent scoring scheme for pooled samples.
<CENTER>

![scoring scheme](https://github.com/staciawyman/blender/blob/master/scoring_scheme.png)

</CENTER>

## Input:

`reference genome`	Path to reference genome. If reference has "mm10" in it, then the mouse blacklist coordinates will be used. Otherwise, human is assumed and the hg38 blacklist coordinates will be used.

`guide sequence`	Guide sequence should be provided 5'-> 3' without the PAM sequence.

`IP bamfile`	This is the aligned bamfile for the MRE11 pulldown of ChIP-Seq of a Cas9 edited sample. BLENDER will extract the reference sequence fromthis file for use in the analysis. I typically use BWA for alignment, but bowtie2 can be used as well. BLENDER has not been tested with bamfiles from other aligners.

`control bamfile`	This is a ChIP-Seq for MRE11 pulldown from either unedited cells or cells that have been edited with a non-targeting gRNA. If there are greater than 10 reads in the control sample, the hit in the edited sample is filtered out. This option can be set by the user.

`output directory`


## Options:

**-p**	List of 2 nucleotide PAM sequences, separated by commas, in quotes. The default is "GG,AG".

**-c**	Cutoff threshold for number of read ends at a putative cut site. Default is 3. For maximum sensitivity, this can be set to 2 and the filtering scheme applied. BEWARE that this dramatically slows down running time. It can increase runtime from ~30min to 24hrs, depending on the guide.

**--verbose** This flag will turn on output of filtered out candidates while running if filtered out for more than maximum mismatches (8) in the guide sequence, or the hit occurs in a blacklist region or it is in a very deep region and thus likely an artifact.


## Output:

run_blender.sh outputs unfiltered_blender_hits.txt and filtered_blender_hits.txt and blender_hits.svg to the output directory provided by the user. This raw unfiltered output can be used for exploring bamfiles to assess whether adjustments might be needed for the scoring scheme. The output text files have the following columns: 

`Chr:Start-End`  Genomic coordinates of the putative guide

`Cutsite`  Where the cutsite is within the guide

`Disco score` Score given to the hit. Essentially summing a window of read ends around the cut site

`Cutsite Ends` This is the number of read ends that pile up at the cutsite. When you set the 'c' parameter, it is this value that is set. 

`Strand/PAM`

`Guide sequence`

`Mismatches`

BLENDER also creates an svg image of the hits something like this.

![blender sample output](https://github.com/staciawyman/blender/blob/master/blender_hits.png)


## Demo

Demo data is included in the DEMO file. The DEMO directory contains data for running BLENDER on a single chromosome (chr19)
for DISCOVER-Seq from MRE11 ChIP-Seq of the VEGFA site 2 guide in K562 cells.

Bam files for the IP file and the control file (BFP ChIP) are provided in the bwa
directory and expected output can be found in the expected_output directory. Chromosome 19
has three off-target hits.

To run BLENDER on the demo data, use the following command in the main blender directory:

	sh run_blender.sh <path/to/reference/genome> \
    		DEMO/bwa/BW43_VEFGA.chr19.bam \
    		DEMO/bwa/BW44_BFP_control.chr19.bam \
    		GACCCCCTCCACCCCGCCTC DEMO/blender

This will take approximately a minute or less to run, and the three output files can be found
in the DEMO/blender directory.

![blender sample output](https://github.com/staciawyman/blender/blob/master/blender_hits_demo.png)

The fastq files to run this data on all chromosomes can be found in the NCBI Short Read Archive
with BioProject Accession PRJNA509652. The links to the fastq files are (click on "Data Access"
tab to get fastq download):

VEGFA IP Fastqs (BW43): https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR8550675

BFP control fastqs (BW44): https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR8550676

Bam files for the data will be available in the near future on the BioProject.

## Citing: 
*Wienert, B., *Wyman, S. K., Richardson, C. D., Yeh, C. D., Akcakaya, P., Porritt, M. J., Morlock, M., Vu, J. T., Kazane, K. R., Watry, H. L., Judge, L. M., Conklin, B. R., Maresca, M. and Corn, J. E. (2019). Unbiased detection of CRISPR off-targets in vivo using DISCOVER-Seq. Science. *contributed equally

