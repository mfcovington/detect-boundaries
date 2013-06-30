#!/usr/bin/env perl
# FILE_NAME.pl
# Mike Covington
# created: 2013-06-30
#
# Description:
#
use strict;
use warnings;
use autodie;
use feature 'say';

my $id = shift;
my @genotyped_files = @ARGV;

my %snps;
for my $file ( @genotyped_files ) {
    open my $snps_fh, "<", $file;
    while (<$snps_fh>) {
        my ( $chr, $pos, $par1, $par2, $tot ) = split /\t/;
        $snps{$chr}{$pos}{par1} += $par1;
        $snps{$chr}{$pos}{par2} += $par2;
        $snps{$chr}{$pos}{tot} += $tot;
    }
    close $snps_fh;
}

for my $chr (sort keys %snps) {
    open my $merged_fh, ">", "$id.$chr.genotyped.merged";
    for my $pos (sort {$a <=> $b} keys $snps{$chr} ) {
        my $par1 = $snps{$chr}{$pos}{par1};
        my $par2 = $snps{$chr}{$pos}{par2};
        my $tot = $snps{$chr}{$pos}{tot};
        say $merged_fh join "\t", $chr, $pos, $par1, $par2, $tot;
    }
    close $merged_fh;
}
