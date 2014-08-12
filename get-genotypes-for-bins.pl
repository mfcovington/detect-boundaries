#!/usr/bin/env perl
# Mike Covington
# created: 2014-07-25
#
# Description:
#
use strict;
use warnings;
use Log::Reproducible;
use autodie;
use feature 'say';
use File::Basename;
use List::Util 'first';

use Data::Printer;

my $usage = "usage: $0 [OUTPUT_FILE] [BINS_FILE] [BOUNDARIES_FILE_1 BOUNDARIES_FILE_2 ...]";
die "$usage\n" unless scalar @ARGV >= 3;

my ( $out_file, $bins_file, @boundaries_files_list ) = @ARGV;
# my $bins_file
#     = "/Users/mfc/Dropbox/best.files/boundaries.best.2014-06-27/boundaries.chr-converted.2014-07-22/bins.tsv";
# my @boundaries_files_list = (
#     "/Users/mfc/Dropbox/best.files/boundaries.best.2014-06-27/boundaries.chr-converted.2014-07-22/BIL_003.boundaries",
#     "/Users/mfc/Dropbox/best.files/boundaries.best.2014-06-27/boundaries.chr-converted.2014-07-22/BIL_005.boundaries",
#     "/Users/mfc/Dropbox/best.files/boundaries.best.2014-06-27/boundaries.chr-converted.2014-07-22/BIL_007.boundaries",
#     "/Users/mfc/Dropbox/best.files/boundaries.best.2014-06-27/boundaries.chr-converted.2014-07-22/BIL_008.boundaries",
#     "/Users/mfc/Dropbox/best.files/boundaries.best.2014-06-27/boundaries.chr-converted.2014-07-22/BIL_009.boundaries"
# );

my $bins = get_bins($bins_file);
# p $bins;
my $boundaries = {};
get_boundaries( $_, $boundaries ) for @boundaries_files_list;

# my $genotypes = {};
my %samples;
for my $chr ( sort keys %$bins ) {
    # say $chr;
    for my $bin_mid ( sort { $a <=> $b } keys %{ $$bins{$chr} } ) {
        # say $bin_mid;
        for my $sample_id ( sort keys %$boundaries ) {
            $samples{$sample_id} = 1;
            # say $sample_id;
            my $geno
                = get_genotype( $sample_id, $chr, $bin_mid, $boundaries );
            # $$genotypes{$sample_id}{$chr}{$bin_mid} = $geno;
            $$bins{$chr}{$bin_mid}{'genotypes'}{$sample_id} = $geno;
        }
    }
}

open my $out_fh, ">", $out_file;
say $out_fh join "\t", 'chr', 'bin-mid', 'bin-start', 'bin-end', sort keys %samples;
for my $chr ( sort keys %$bins ) {
    for my $bin_mid ( sort { $a <=> $b } keys %{ $$bins{$chr} } ) {
        my $start = $$bins{$chr}{$bin_mid}{'start'};
        my $end = $$bins{$chr}{$bin_mid}{'end'};
        my @genotypes = map { $$bins{$chr}{$bin_mid}{'genotypes'}{$_} }
            sort keys %{ $$bins{$chr}{$bin_mid}{'genotypes'} };
        say $out_fh join "\t", $chr, $bin_mid, $start, $end, @genotypes;
    }
}
close $out_fh;

# p $genotypes;
# p $bins;
sub get_genotype {
    my ( $sample_id, $chr, $bin_mid, $boundaries ) = @_;

    my $start = first { $_ <= $bin_mid }
    sort { $b <=> $a } keys %{ $$boundaries{$sample_id}{$chr} };
    my $end = $$boundaries{$sample_id}{$chr}{$start}{'end'};

    my $genotype
        = $bin_mid <= $end
        ? $$boundaries{$sample_id}{$chr}{$start}{'geno'}
        : 'NA';

# say join "\t", $sample_id, $chr, $start, $bin_mid, $end, $genotype;
die unless defined $start;
# exit;
    # my $genotype;
    return $genotype;
}

sub get_bins {
    my $bins_file = shift;

    my %bins;
    open my $bins_fh, "<", $bins_file;
    <$bins_fh>;
    while (<$bins_fh>) {
        chomp;
        my ( $chr, $start, $end ) = split;
        my $mid = int( ( $start + $end ) / 2 );
        $bins{$chr}{$mid}{'start'} = $start;
        $bins{$chr}{$mid}{'end'} = $end;
    }
    close $bins_fh;

    return \%bins;
}

sub get_boundaries {
    my ( $bounds_file, $boundaries ) = @_;

    my $sample_id = get_sample_id($bounds_file);

    open my $bounds_fh, "<", $bounds_file;
    while (<$bounds_fh>) {
        chomp;
        my ( $chr, $start, $end, $genotype ) = split;
        $$boundaries{$sample_id}{$chr}{$start}{'end'} = $end;
        $$boundaries{$sample_id}{$chr}{$start}{'geno'} = $genotype;
    }
    close $bounds_fh;
}

sub get_sample_id {
    my $bounds_file = shift;
    my ($sample_id) = fileparse( $bounds_file, ".boundaries" );
    return $sample_id;
}

# sub find_bin {
#     my ( $chr, $start, $bins ) = @_;

#     my $bin_start = first { $_ <= $start } sort {$a <=> $b } keys %{ $$bins{$chr} };

#     my $bin_mid;# = $$bins{$chr}{$bin_start}{'mid'};

#     return $bin_mid;
# }
