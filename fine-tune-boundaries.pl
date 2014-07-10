#!/usr/bin/env perl
# Mike Covington
# created: 2014-07-09
#
# Description:
#
use strict;
use warnings;
# use Log::Reproducible;
use autodie;
use feature 'say';
use Getopt::Long;
use File::Basename;
use Data::Printer;

my $genotyped_dir;
my $context = 10;

my $options = GetOptions(
    "genotyped_dir=s" => \$genotyped_dir,
    "context=i"       => \$context,
);

my @boundaries_files = @ARGV;

for my $bounds_file (@boundaries_files) {
    my $sample_id  = get_sample_id($bounds_file);
    my $boundaries = get_boundaries($bounds_file);

    # say $sample_id;
    # p $boundaries;
}

sub get_sample_id {
    my $bounds_file = shift;
    my ($sample_id) = fileparse( $bounds_file, ".boundaries" );
    return $sample_id;
}

sub get_boundaries {
    my $bounds_file = shift;
    my %boundaries;
    open my $bounds_fh, "<", $bounds_file;
    while (<$bounds_fh>) {
        chomp;
        my ( $chr, $start, $end, $geno ) = split;
        push @{ $boundaries{$chr} },
            { start => $start, end => $end, geno => $geno };
    }
    close $bounds_fh;
    return \%boundaries;
}
