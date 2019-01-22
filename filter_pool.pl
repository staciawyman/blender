#!/usr/bin/perl
# 11/11/18 Stacia Wyman staciawyman@berkeley.edu
# Use this filter when samples have been pooled together
# It is more strict than the single sample filter

while(<>) {
    if ($_ !~ /^chr/) { print; } # print header line 
    ($coord,$cut,$disco,$ends,$pam,$guide,$mm) = split(/\t/);
    if ($mm <= 7 && $disco >= 5) { print; next;}
    if ($mm <= 5 && $disco >= 4) { print; next;}
    if ($mm <= 3 && $disco >= 3) { print; next;}
}
