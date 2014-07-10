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
    my $genotypes  = get_genotypes( $genotyped_dir, $sample_id );

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
            { 'start' => $start, 'end' => $end, 'geno' => $geno };
    }
    close $bounds_fh;
    return \%boundaries;
}

sub get_genotypes {
    my ( $genotyped_dir, $sample_id ) = @_;
    my %genotypes;
    my @genotyped_files = glob "$genotyped_dir/$sample_id.*.genotyped*";
    for my $geno_file (@genotyped_files) {
        open my $geno_fh, "<", $geno_file;
        while ( my $line = <$geno_fh> ) {
            chomp $line;
            my ( $chr, $pos ) = split /\t/, $line;
            $genotypes{$chr}{$pos} = $line;
        }
        close $geno_fh;
    }
    return \%genotypes;
}
