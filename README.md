# **BLENDER**
## About

BLENDER is a companion program to the DISCOVER-Seq assay to identify off-target editing sites from MRE11 ChIP-Seq experiments.  It takes aligned bamfiles from the IP experiment and (optionally) a control bamfile and identifies locations with stacks of reads at putative cutsites. PAM sequences can be provided by the user as well as a guide sequence. BLENDER makes use of the ENCODE ChIP-Seq Blacklists for human and mouse which are lists of regions in the genome that have artifactual large numbers of reads. These lists and the control bam plus PAM sequences and the guide are used to filter out false positives.  BLENDER runs on mouse mm10 and human hg38 genomes (blacklists coordinates are for these genomes). 

## Usage:
`perl blender.pl [options] guide_sequence IP_bamfile control_bamfile  > unfiltered_output.txt`
`perl blender.pl [options] <guide sequence> <IP bamfile> <control bamfile>  | perl filter.pl > output.txt`

## Input:
	guide sequence	Guide sequence should be provided 5'-> 3' without the PAM sequence.
	IP bamfile	This is the aligned bamfile for the MRE11 pulldown ChIP-Seq sample. BLENDER will extract the reference sequence fromthis file for use in the analysis. I typically use BWA for alignment, but bowtie2 can be used as well. BLENDER has not been tested with bamfiles from other aligners.

	control bamfile	This is a bamfile from ChIP-Seq with BFP or GFP or similar pulldown and is used as to filter false positives. If there are greater than 10 reads in the control sample, the hit is filtered out. This option can be set by the user.


## Options:
	-p	List of 2 nucleotide PAM sequences, separated by commas, in quotes. The default is "GG,AG".
	-t	Threshold for number of read ends at a putative cut site. Default is 3. For maximum sensitivity, this can be set to 2 and the filtering scheme applied. BEWARE that this dramatically slows down running time. It can increase runtime from ~30min to 24hrs, depending on the guide.


## Output:
	The output of blender.pl is unfiltered. This output can be used for exploring bamfiles to assess whether adjustments might be needed for the scoring scheme. Alternatively, the output of blender.pl can be directly piped into filter.pl to apply scoring scheme and get a list of filtered results. The output has the following columns: Chr:Start-End	Cutsite	Discoscore	Cutsite Ends	Strand/PAM	Guide sequence	Mismatches

## Citing: TBD
