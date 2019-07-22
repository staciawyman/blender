# **BLENDER**
## About

BLENDER is a companion program to the DISCOVER-Seq assay to identify off-target editing sites from MRE11 ChIP-Seq experiments.  It takes aligned bamfiles from the IP experiment and (optionally) a control bamfile and identifies locations with stacks of reads at putative cutsites. PAM sequences can be provided by the user as well as a guide sequence. BLENDER makes use of the ENCODE ChIP-Seq Blacklists for human and mouse which are lists of regions in the genome that have artifactual large numbers of reads. These lists and the control bam plus PAM sequences and the guide are used to filter out false positives.  BLENDER runs on mouse mm10 and human hg38 genomes (blacklists coordinates are for these genomes). 
BLENDER has a driver script run_blender.sh that takes 5 (optionally 6) argmuments and runs identification of putative hits, filtering of those hits, and then creates an SVG of the aligned hits. Files are stored in the output directory given as a parameter to the bash script.
Alternatively, each step can be run separately. 

![schematic](https://github.com/staciawyman/blender/blob/master/BLENDER_schematic.png)

## Running the bash script:

>sh run_blender.sh <path to reference genome> <path to IP bamfile> \
>	<path to control bamfile> \
>       <guide sequence> <output directory> ["options"]

## Usage:

`perl blender.pl [options] <reference genome> <guide sequence> <IP bamfile> <control bamfile>  > unfiltered_output.txt`

`perl blender.pl [options] <reference genome> <guide sequence> <IP bamfile> <control bamfile>  | perl filter.pl > output.txt`

`perl blender.pl [options] <reference genome> <guide sequence> <IP bamfile> <control bamfile>  | perl filter_pool.pl > pooled_output.txt`

BLENDER can be run with or without being piped through the filtering script. There are two filtering scripts provided; the standard filter.pl script that implements the standard scoring scheme, and the filter_pool.pl script that implements the more stringent scoring scheme for pooled samples.
<CENTER>

![scoring scheme](https://github.com/staciawyman/blender/blob/master/scoring_scheme.png)

</CENTER>

## Input:

`genome`	Path to reference genome. If reference has "mm10" in it, then the mouse blacklist coordinates will be used. Otherwise, human is assumed and the hg38 blacklist coordinates will be used.

`guide sequence`	Guide sequence should be provided 5'-> 3' without the PAM sequence.

`IP bamfile`	This is the aligned bamfile for the MRE11 pulldown of ChIP-Seq of a Cas9 edited sample. BLENDER will extract the reference sequence fromthis file for use in the analysis. I typically use BWA for alignment, but bowtie2 can be used as well. BLENDER has not been tested with bamfiles from other aligners.

`control bamfile`	This is a ChIP-Seq for MRE11 pulldown from either unedited cells or cells that have been edited with a non-targeting gRNA. If there are greater than 10 reads in the control sample, the hit in the edited sample is filtered out. This option can be set by the user.


## Options:

**-p**	List of 2 nucleotide PAM sequences, separated by commas, in quotes. The default is "GG,AG".

**-c**	Cutoff threshold for number of read ends at a putative cut site. Default is 3. For maximum sensitivity, this can be set to 2 and the filtering scheme applied. BEWARE that this dramatically slows down running time. It can increase runtime from ~30min to 24hrs, depending on the guide.

**--verbose** This flag will turn on output of filtered out candidates while running if filtered out for more than maximum mismatches (8) in the guide sequence, or the hit occurs in a blacklist region or it is in a very deep region and thus likely an artifact.


## Output:

blender.pl outputs to stdout and the output is unfiltered. This raw output can be used for exploring bamfiles to assess whether adjustments might be needed for the scoring scheme. Alternatively, the output of blender.pl can be directly piped into filter.pl to apply scoring scheme and get a list of filtered results. The output has the following columns: 

`Chr:Start-End`

`Cutsite`

`Discoscore`

`Cutsite Ends`

`Strand/PAM`

`Guide sequence`

`Mismatches`

## Citing: 
*Wienert, B., *Wyman, S. K., Richardson, C. D., Yeh, C. D., Akcakaya, P., Porritt, M. J., Morlock, M., Vu, J. T., Kazane, K. R., Watry, H. L., Judge, L. M., Conklin, B. R., Maresca, M. and Corn, J. E. (2019). Unbiased detection of CRISPR off-targets in vivo using DISCOVER-Seq. Science. *contributed equally

