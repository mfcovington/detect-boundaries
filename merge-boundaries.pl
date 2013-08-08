#!/usr/bin/env perl
# merge-boundaries.pl
# Mike Covington
# created: 2013-08-07
#
# Description:
#
use strict;
use warnings;
use autodie;
use feature 'say';
use Number::RangeTracker;
use Data::Printer;

my @boundary_files = @ARGV; # || qw(sample-file/RIL_300.boundaries sample-file/RIL_300.boundaries);

# TODO: Extract chromosome lengths from bam file
# TODO: GetOptions: Which chromosomes we want to use (if a subset of .bam header)
my %chr_lengths = (
    A01 => 28608137,
    A02 => 27848129,
    A03 => 31716688,
    A04 => 18967243,
    A05 => 23941934,
    A06 => 26273242,
    A07 => 22587824,
    A08 => 21596550,
    A09 => 37123981,
    A10 => 17595035
);

my @chromosomes = sort keys %chr_lengths;

my %range;
for my $chr (@chromosomes) {
    $range{$chr} = RangeTracker->new();
}
p @boundary_files;
for my $file (@boundary_files) {
    open my $boundary_fh, "<", $file;
    my %boundaries;
    while (<$boundary_fh>) {
        my ($chr, $start, $end ) = split;
        $boundaries{$chr}{$start} = $end;
    }

    for my $chr (@chromosomes) {
        my $end = 0;
        for my $start ( sort { $a <=> $b } keys $boundaries{$chr} ) {
            $range{$chr}->add_range($end + 1, $start - 1);
            $end = $boundaries{$chr}{$start};
        }
        $range{$chr}->add_range($end + 1, $chr_lengths{$chr});
    }

}

my %borders;
my %bins;
for my $chr (@chromosomes) {
    %{$borders{$chr}} = $range{$chr}->output_ranges;
    %{$bins{$chr}} = $range{$chr}->inverse;
}
p %borders;
p %bins;

# TODO: Output min/max/mean border/bin size
# TODO: Write bin boundaries to file
say "BORDER SIZES:";
for my $chr (@chromosomes) {
    say "# of borders in $chr: " . scalar keys $borders{$chr};
    for my $start ( sort { $a <=> $b } keys %{$borders{$chr}} ) {
            say "\t"  . $borders{$chr}->{$start} - $start + 1;
    }
}

say "BIN SIZES:";
for my $chr (@chromosomes) {
    say "# of bins in $chr: " . scalar keys $bins{$chr};
    for my $start ( sort { $a <=> $b } keys %{$bins{$chr}} ) {
            say "\t"  . $bins{$chr}->{$start} - $start + 1;
    }
}









exit;
