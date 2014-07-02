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

my ( $query_genotype, $target_genotype, $buffer, @boundaries_files ) = @ARGV;

for my $query_file (@boundaries_files) {
    open my $query_fh, "<", $query_file;
    while (<$query_fh>) {
        chomp;
        my ( $chr, $start, $end, $genotype ) = split;
        next unless $genotype eq $query_genotype;
        find_similar_bins( $chr, $start, $end, $target_genotype, $buffer,
            $query_file, \@boundaries_files );
    }
}

sub find_similar_bins {
    my ( $query_chr, $query_start, $query_end, $target_genotype, $buffer,
        $query_file, $boundaries_files )
        = @_;

    for my $file (@$boundaries_files) {
        open my $boundaries_fh, "<", $file;
        while (<$boundaries_fh>) {
            chomp;
            my ( $chr, $start, $end, $genotype ) = split;
            next unless $chr      eq $query_chr;
            next unless $genotype eq $target_genotype;
            my $start_diff = abs( $start - $query_start );
            my $end_diff   = abs( $end - $query_end );
            if ( $start_diff < $buffer && $end_diff < $buffer ) {
                say join "\t", $query_file, $query_chr, $query_start,
                    $query_end, $query_genotype, $file, $chr, $start, $end,
                    $genotype;
            }
        }
        close $boundaries_fh;
    }
}

