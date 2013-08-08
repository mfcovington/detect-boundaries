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

use lib "/Users/mfc/git.repos/ranges/lib";

use Number::RangeTracker;
use Data::Printer;

# my $range = RangeTracker->new();

my @boundary_files = @ARGV; # || qw(sample-file/RIL_300.boundaries sample-file/RIL_300.boundaries);

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
    # my $range = RangeTracker->new();
    my %boundaries;
    while (<$boundary_fh>) {
        my ($chr, $start, $end ) = split;
        $boundaries{$chr}{$start} = $end;
    }
# p %boundaries;


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

say "BORDER SIZES:";
for my $chr (@chromosomes) {
    say "# of borders in $chr: " . scalar keys $borders{$chr};
    for my $start ( sort { $a <=> $b } keys %{$borders{$chr}} ) {
        # p $range_ref;
        # for my $start (keys sort { $a <=> $b } %$range_ref ) {
            say "\t"  . $borders{$chr}->{$start} - $start + 1;
        # }
    }
}

say "BIN SIZES:";
for my $chr (@chromosomes) {
    say "# of bins in $chr: " . scalar keys $bins{$chr};
    for my $start ( sort { $a <=> $b } keys %{$bins{$chr}} ) {
        # p $range_ref;
        # for my $start (keys sort { $a <=> $b } %$range_ref ) {
            say "\t"  . $bins{$chr}->{$start} - $start + 1;
        # }
    }
}



# for my $chr (@chromosomes) {

# }

__END__

# $range->add_range( 1, 10 );
# $range->add_range( 5, 15 );
# $range->add_range( 50, 150 );

my %borders = $range->output_ranges;

use Data::Printer;
p %borders;



