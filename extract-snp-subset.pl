#!/usr/bin/env perl
# extract-snp-subset.pl
# Mike Covington
# created: 2013-08-08
#
# Description:
#
use strict;
use warnings;
use autodie;
use feature 'say';
use Data::Printer;

my %ratios;

my $flank       = 0.2;
my $flank_extra = $flank + 0.1;

my $obs_rat_min       = 0.9;
my $obs_rat_min_extra = 0.8;

my $vcf_summary_file = "/Users/mfc/git.repos/snps-from-rils/merged.EC50.minflter.vcf/summaries.het_ratio_0_1.alt_ratmin_0_05.filters/merged.A01.EC50.minflter.vcf.summary";
open my $vcf_summary_fh, "<", $vcf_summary_file;
<$vcf_summary_fh>
while (<$vcf_summary_fh>) {
    next unless m/^\./;
    my ( $chr, $pos, $obs_rat ) = (split)[ 1, 2, 10 ];
    $ratios{$chr}{$pos} = $obs_rat;
}
close $vcf_summary_fh;

open my $bins_fh, "<", "sample-file/bins.tsv";
<$bins_fh>;
while (<$bins_fh>) {
    chomp;
    my ( $chr, $start, $end ) = split;

    my %obs_ratios;
    for my $pos (sort { $a <=> $b } keys $ratios{$chr} ) {
        next if $pos < $start;
        next if $pos > $end;
        $obs_ratios{$pos} = $ratios{$chr}{$pos};
    }

    my @snps   = sort { $a <=> $b } keys %obs_ratios;
    my $left   = $snps[0];
    my $right  = $snps[-1];
    my $length = $right - $left + 1;

    my $flank_adj       = $length * $flank;
    my $flank_extra_adj = $length * $flank_extra;

    my $left_flank        = sprintf '%0.0f', $left  + $flank_adj;
    my $right_flank       = sprintf '%0.0f', $right - $flank_adj;
    my $left_flank_extra  = sprintf '%0.0f', $left  + $flank_extra_adj;
    my $right_flank_extra = sprintf '%0.0f', $right - $flank_extra_adj;

    my $left_best;
    my $left_best_extra;
    my @left_max = ( 0, 0 );
    for my $pos (@snps) {
        my $obs_rat = $obs_ratios{$pos};

        next if $pos > $left_flank_extra;

        if ($left_max[1] < $obs_rat) {
            @left_max = ($pos, $obs_rat);
        }

        if ( $obs_rat > $obs_rat_min ) {
            $left_best = $pos;
            last;
        }

        if ( $obs_rat > $obs_rat_min_extra ) {
            $left_best_extra //= $pos;
        }

        if ( $pos > $left_flank ) {
            if ( $obs_rat > $obs_rat_min ) {
                $left_best = $pos;
                last;
            }

            if ( $obs_rat > $obs_rat_min_extra ) {
                $left_best_extra //= $pos;
            }
        }
    }    # for my $pos (@snps)

    # say "$left_best:$left_best_extra";
    say "@snps";
    say "$left to $right";
    say "@left_max";
    say $length;
    say join "\t", $left_flank, $left_flank_extra, $right_flank, $right_flank_extra;
    exit;
}
close $bins_fh;


