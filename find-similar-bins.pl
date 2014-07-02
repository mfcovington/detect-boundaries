#!/usr/bin/env perl
# Mike Covington
# created: 2014-07-02
#
# Description:
#
use strict;
use warnings;
use Log::Reproducible;
use autodie;
use feature 'say';

my ( $chr, $start, $end, $buffer, @boundaries_files ) = @ARGV;

for my $file (@boundaries_files) {
    open my $boundaries_fh, "<", $file;
    while (<$boundaries_fh>) {
        chomp;
        my ( $current_chr, $current_start, $current_end, $genotype ) = split;
        next unless $chr eq $current_chr;
        my $start_diff = abs( $start - $current_start );
        my $end_diff   = abs( $end - $current_end );
        if ( $start_diff < $buffer && $end_diff < $buffer ) {
            say join "\t", $file, $current_chr, $current_start, $current_end,
                $genotype;
        }
    }
    close $boundaries_fh;
}

