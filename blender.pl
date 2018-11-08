#!/usr/bin/perl

my %both_starts;

$debug = 0;
$verbose = 0;
$control_bam = "";
$gM = "GGCTGATGAGGCCGCACATG";
$gP = "AGCAGCAGCGGCGGCAACAG";
$vegfa = "GACCCCCTCCACCCCGCCTC";
if ($#ARGV < 2) { print "Missing argumnet $#ARGV \n"; exit; }
if (length($ARGV[0]) != 20) { print "Please provide a 20bp guide sequence.$ARGV[0]\n"; exit; }
$input_guide = $ARGV[0];
$edited_bam = $ARGV[1];
$control_bam = $ARGV[2];
@pamlist = ( "GG","AG" );
$cas12a = 0;
# Test params
$check_guide = 1;
$max_mismatches = 8;
$min_discoscore = 3;
$threshold = 3;

# Get organism/location of reference genome from @PG line in bamfile header
$PG = `samtools view -H $edited_bam | grep PG`;
if ($PG =~ /sampe\s+(.+)\s+/) {
    #print "Genome is $1\n";
    $genome = $1;
}

if ($genome =~ /mm10/) {
#    $genome = "/data/genomes/mouse/mm10/mm10.fa";
    $blacklist_file = "mm10.blacklist.bed";
    @chroms = (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,"X","Y");
} else { # assume human
#    $genome = "/data/genomes/GRCh38noM.mfa"; # get genome from PG line in bam file
    $blacklist_file = "hg38.blacklist.bed";
    @chroms = (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y");
}

print "Chr:Start-End	Cutsite	Discoscore	Cutsite Ends	Strand/PAM	Guide sequence	Mismatches\n";
# Get species/num chrom from bam file
# for (my $i = 16; $i <= 16; $i++) { # narrow to one chrom for quick debugging
foreach $i (@chroms) {
    undef %rev_starts;
    undef %for_starts;
    undef %both_starts;
    $c = "chr".$i;
    open(F,"samtools view $edited_bam $c |");
    #print "Processing chrom $c\n";
    while (<F>) {
	chomp;
        my ($read,$flag,$chr,$start_loc,$MQ,$cigar,$a,$b,$tlen,$seq) = split(/\t/);
	my $seqlen = length($seq) ;
	if ($tlen == 0 && $flag & 0x08 && $MQ > 25) {
	    if (!($flag & 0x10)) { 
                $for_starts{$start_loc}++;
	        $both_starts{$start_loc}++;
	    } else {
                $rev_starts{$start_loc+$seqlen-1}++;
	        $both_starts{$start_loc+$seqlen-1}++;
	    }
	}
	if ($MQ > 25 && $tlen > 0) {
            $for_starts{$start_loc}++;
	    $both_starts{$start_loc}++;
        }
	if ($MQ > 25 && $tlen < 0) {
            $rev_starts{$start_loc+$seqlen-1}++;
	    $both_starts{$start_loc+$seqlen-1}++;
        }
    }
    close(F);
    # Combine forw and rev blunt ends at cut site
    foreach my $start (sort { $a <=> $b } keys(%for_starts) ) {
	if (($rev_starts{$start-1} + $for_starts{$start}) >= $threshold) {
	    $both_starts{$start-1} += $both_starts{$start};
	    #$both_starts{$start} = 0;
	}
    }

    foreach my $start (sort { $a <=> $b } keys(%both_starts) ) {
	$l = "$c:$start-$start";
        if ($both_starts{$start} >= $threshold) {
        #if ($both_starts{$start} >= $threshold || 
	#($both_starts{$start} == 2 && ($both_starts{$start-1}+$both_starts{$start}+$both_starts{$start+2}) >= $threshold)) {
	    if ($debug) {print "Check both_starts $c:$start, $both_starts{$start}\n";}
	    if (blacklist($c,$start)) { 
    	        if ($verbose) { print "$c:$start\t$sum\t$both_starts{$start}\tWARNING:blacklisted\n"; }
		next; 
	    }
	    if ($control_bam) {
                $d_control = `samtools depth -r $l $control_bam`;
	        chomp($d_control);
	        ($x,$y,$depth_con) = split(/\t/,$d_control);
	        if ($depth_con > 10) { 
		    if ($verbose) {print "CONTROL Skipping $l, $depth_con, $both_starts{$start} \n";}
		    next;
	        }
	    } 

	    $over_max = 0;

	    $pamleft = check_pam_left($c,$start);
	    $pamright = check_pam_right($c,$start);

    	    if ($pamleft) { 
                $d = `samtools depth -r $l $edited_bam`;
	        chomp($d);
	        ($x,$y,$depth) = split(/\t/,$d);
	 	if ($depth > 0 ) {
	            if ($both_starts{$start}/$depth < .25) { 
    	                if ($verbose){ print "$c:$s-$e\t$start\t$sum\t$both_starts{$start}\tPAM-left($pamleft)\t$guide\tWARNING:deep area\n";}
			$over_max = 1;
		    }
		}
		if (!$over_max) {
    	        $s = $start - 2; $e = $start + 17;
    		$guide = get_guide($c,$s,$e);
    		$sum = add_window($start,5);
		if ($sum < $min_discoscore) { next;} 
		my $mm = guide_mm($input_guide,revcomp($guide));
		if ($check_guide)  {
			if ($mm <= $max_mismatches) {
    	                    print "$c:$s-$e\t$start\t$sum\t$both_starts{$start}\tPAM-left($pamleft)\t$guide\t$mm\n";
			} else {
    	                    if ($verbose) {print "$c:$s-$e\t$start\t$sum\t$both_starts{$start}\tPAM-left($pamleft)\t$guide FILTERED: $mm mismatches\n";}
			}
			
    		} else {
    	                print "$c:$s-$e\t$start\t$sum\t$both_starts{$start}\tPAM-left($pamleft)\t$guide\t$mm\n";
		}
		}
    	    }
    	    if ($pamright) { 
                $d = `samtools depth -r $l $edited_bam`;
	        chomp($d);
	        ($x,$y,$depth) = split(/\t/,$d);
	 	if ($depth > 0 ) {
	            if ($both_starts{$start}/$depth < .25) { 
    	                if ($verbose) {print "$c:$s-$e\t$start\t$sum\t$both_starts{$start}\tPAM-right($pamright)\t$guide\tWARNING:deep area\n";}
			$over_max = 1;
			next;
		    }
		}
		if (!$over_max) {
    	            $e = $start + 3; $s = $start - 16;
    		    $guide = get_guide($c,$s,$e);
    		    $sum = add_window($start,5);
		    if ($sum < $min_discoscore) { next;} 
		    my $mm = guide_mm($input_guide,$guide);
		    if ($check_guide)  { # guide was given as input parameter
			if ($mm <= $max_mismatches) {
    	                    print "$c:$s-$e\t$start\t$sum\t$both_starts{$start}\tPAM-right($pamright)\t$guide\t$mm\n";
			} else {
    	                    if ($verbose) {print "$c:$s-$e\t$start\t$sum\t$both_starts{$start}\tPAM-right($pamright)\t$guide FILTERED: $mm mismatches\n";}
		 	    next;
			}
		    } else {
    	                print "$c:$s-$e\t$start\t$sum\t$both_starts{$start}\tPAM-right($pamright)\t$guide\t$mm\n";
		    }
		}
    	    }
	}
    }
}

sub check_pam_left {
    my ($chr,$x) = @_;
    my $ref_pam = "";
    $s = $x-5; $e = $x-4;
    my $coords = $chr.":".$s."-".$e;
    
    open(F,"samtools faidx $genome $coords  | ");
    while (<F>) {
	if (/^>/) { next; }
	chomp;
	$ref_pam = uc($_);
    }
    close(F);
   
    foreach my $pam (@pamlist) {
      if ($ref_pam eq revcomp($pam)) { return $pam; }
    }
    return "";
}
	

# alt pam is AG
sub check_pam_right {
    my ($chr,$x) = @_;
    
    my $ref_pam = "";
    $s = $x+5; $e = $x+6;
    $coords = $chr.":".$s."-".$e;
    open(F,"samtools faidx $genome $coords  | ");
    while (<F>) {
	chomp;
	if (/^>/) { next; }
	$ref_pam = uc($_);
    }
    close(F);
    foreach my $pam (@pamlist) {
      if ($ref_pam eq $pam) { return $pam; }
    }
    return "";
}

sub get_guide {
    my ($chr,$start,$end) = @_;

    $coords = $chr.":".$start."-".$end;
    open(F,"samtools faidx $genome $coords | ");
    while (<F>) {
	chomp;
	if (/^>/) { next; } 
	return uc($_);
    }
}

sub get_guide_from_bam {
    my ($chr,$start,$end) = @_;

    $coords = $chr.":".$start."-".$end;
    open(F,"samtools mpileup -r $coords -f $genome $edited_bam 2>/dev/null | ");
    $str = "";
    while (<F>) {
	chomp;
	($chr,$loc,$nt,$depth) = split(/\t/);
	$str .= $nt;
    }
    close(F);
    return (uc($str));
}


# add window plus 1bp over cutsite for forward and for reverse
sub add_window {
    my ($start,$window_size) = @_;
    $sum = 0;
    for (my $i = $start-($window_size-1); $i <= $start+1; $i++) {
	$sum += $rev_starts{$i};
    } 
    for (my $i = $start; $i <= $start+$window_size; $i++) {
	$sum += $for_starts{$i};
    } 
    return $sum;
}

sub blacklist {
    my ($chr,$start) = @_;
    open(F,$blacklist_file) || die "couldn't open $blacklist_file";
    while (<F>) {
	chomp;
	($b_chr,$b_start,$b_end) = split(/\t/);
	if (($chr eq $b_chr) && ($start >= $b_start) && ($start <= $b_end)) {
	    #print "FILTERED BLACKLIST $chr;$start\n";
	    return 1;
	}
    }
    return 0;
}

# Check how many mismatches between input guide and guide sequence at target site
sub guide_mm {
    my ($input_guide,$guide) = @_;
    my $mm = 0;
    (@input) = split(//,uc($input_guide));
    (@check) = split(//,uc($guide));


    for (my $i = 0; $i <= $#input; $i++) {
	if ($input[$i] ne $check[$i]) { $mm++; }
    }
    return $mm;
}

sub revcomp {
    my ($dna) = @_;
    my $rc = reverse($dna);
    $rc =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
    return $rc;
}
