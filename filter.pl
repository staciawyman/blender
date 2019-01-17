#!/usr/bin/perl
# 11/11/18 Stacia Wyman staciawyman@berkeley.edu

while(<>) {
    if ($_ !~ /^chr/) { print; } # print header line 
    ($coord,$cut,$disco,$ends,$strand,$pam,$guide,$mm) = split(/\t/);
    if ($mm <= 7 && $disco >= 4) { print; next;}
    if ($mm <= 5 && $disco >= 3) { print; next;}
    if ($mm <= 3 && $disco >= 2) { print; next;}
}
