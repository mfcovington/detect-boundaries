#!/usr/bin/env perl
# FILE_NAME.pl
# Mike Covington
# created: 2013-08-09
#
# Description:
#
use strict;
use warnings;
use autodie;
use feature 'say';

open my $bins_fh, "<", "sample-file/bins-snp.mid";
<$bins_fh>;
my %mids;
while (<$bins_fh>) {
    chomp;
    my ( $chr, $pos ) = (split)[0, 3];
    next if $pos eq 'NA';
    $mids{$chr}{$pos} = 1;
}

open my $snps_fh, "<", "sample-file/RIL_93.A01.filtered.snps";
my %snps;
while (<$snps_fh>) {
    chomp;
    my ( $chr, $pos, $genotype ) = (split)[0, 1, 3];
    $snps{$chr}{$pos} = $genotype;
}

# my $chr = 'A01';
for my $chr (sort keys %mids ) {
    for my $pos (sort { $a <=> $b } keys $mids{$chr} ) {
        $mids{$chr}{$pos} = $snps{$chr}{$pos} // 'NA';
    }
}



use Data::Printer;
p %mids;




