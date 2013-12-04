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
use List::Util qw(min max sum);
use Getopt::Long;

my @boundary_files = @ARGV; # || qw(sample-file/RIL_300.boundaries sample-file/RIL_300.boundaries);

# Temporary defaults values:
my $chr_list = "A01,A02,A03,A04,A05,A06,A07,A08,A09,A10";
my $bam_file = "~/git.repos/sample-files/bam/IMB211.good.bam";

my $options = GetOptions(
    "chr_list=s" => \$chr_list,
    "bam_file=s" => \$bam_file,
);

my @chromosomes = split /,/, $chr_list;
my $chr_lengths = get_chr_lengths( $bam_file, \@chromosomes );

my %range;
for my $chr (@chromosomes) {
    $range{$chr} = RangeTracker->new();
}

for my $file (@boundary_files) {
    open my $boundary_fh, "<", $file;
    my %boundaries;
    while (<$boundary_fh>) {
        my ( $chr, $start, $end ) = split;
        $boundaries{$chr}{$start} = $end;
    }

    for my $chr (@chromosomes) {
        my $end = 0;
        for my $start ( sort { $a <=> $b } keys $boundaries{$chr} ) {
            $range{$chr}->add_range( $end + 1, $start - 1 );
            $end = $boundaries{$chr}{$start};
        }
        $range{$chr}->add_range( $end + 1, $chr_lengths->{$chr} );
    }

}

my %borders;
my %bins;
for my $chr (@chromosomes) {
    %{ $borders{$chr} } = $range{$chr}->output_ranges;
    %{ $bins{$chr} }    = $range{$chr}->inverse;
}

# TODO: Write bin boundaries to file
say "BORDER SIZES:";
for my $chr (@chromosomes) {
    my @lengths;
    for my $start ( sort { $a <=> $b } keys %{ $borders{$chr} } ) {
        push @lengths, $borders{$chr}->{$start} - $start + 1;
    }
    summarize_ranges( $chr, \@lengths );
}

say "BIN SIZES:";
open my $bin_fh, ">", "bins.tsv";
say $bin_fh join "\t", "chr", "start", "end";
for my $chr (@chromosomes) {
    my @lengths;
    for my $start ( sort { $a <=> $b } keys %{ $bins{$chr} } ) {
        my $end = $bins{$chr}->{$start};
        push @lengths, $end - $start + 1;
        say $bin_fh join "\t", $chr, $start, $end;
    }
    summarize_ranges( $chr, \@lengths );
}
close $bin_fh;

exit;

sub get_chr_lengths {
    my ( $bam_file, $chromosomes ) = @_;

    my %chr_lengths = map { $_ => 0 } @{$chromosomes};
    my $chr_count = scalar @{$chromosomes};

    open my $bam_head_fh, "-|", "samtools view -H $bam_file";
    while (<$bam_head_fh>) {
        next unless /^\@SQ/;

        my ( $seq_id, $seq_len ) = $_ =~ m/\tSN:([^\t]+)\tLN:(\d+)/;

        next    # only skip 'nonexistent' chromosomes if custom list is supplied
          if $chr_count > 0
          && !exists $chr_lengths{$seq_id};

        $chr_lengths{$seq_id} = $seq_len;
    }
    close $bam_head_fh;

    return \%chr_lengths;
}

sub summarize_ranges {
    my ( $chr, $lengths ) = @_;

    my $count = @$lengths;
    my $min   = min @$lengths;
    my $mean  = sprintf( '%.0f', sum(@$lengths) / $count );
    my $max   = max @$lengths;

    my $summary = <<END_SUMMARY;
$chr
  bins: $count
  min:  $min
  mean: $mean
  max:  $max
END_SUMMARY

    say $summary;
}
