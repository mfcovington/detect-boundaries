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
use Data::Printer;

my @boundary_files = @ARGV; # || qw(sample-file/RIL_300.boundaries sample-file/RIL_300.boundaries);

my @chr_list = qw( A01 A02 A03 A04 A05 A06 A07 A08 A09 A10 );


my $bam_file = "~/git.repos/sample-files/bam/IMB211.good.bam";
my $chr_lengths = get_chr_lengths( $bam_file, \@chr_list );

p $chr_lengths;



sub get_chr_lengths {
    my ( $bam_file, $chr_list_ref ) = @_;

    my %chr_lengths = map { $_ => 0 } @{$chr_list_ref};

    open my $bam_head_fh, "-|", "samtools view -H $bam_file";
    while (<$bam_head_fh>) {
        next unless /^\@SQ/;
        chomp;
        my ( $seq_id, $seq_len ) = $_ =~ m/\tSN:([^\t]+)\tLN:([^\t]+)/;
        next unless exists $chr_lengths{$seq_id};
        $chr_lengths{$seq_id} = $seq_len;
    }
    close $bam_head_fh;

    return \%chr_lengths;
}


exit;
__END__

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

# TODO: Write bin boundaries to file
say "BORDER SIZES:";
for my $chr (@chromosomes) {
    my @lengths;
    for my $start ( sort { $a <=> $b } keys %{$borders{$chr}} ) {
        push @lengths, $borders{$chr}->{$start} - $start + 1;
    }
    summarize_ranges( $chr, \@lengths );
}

say "BIN SIZES:";
open my $bin_fh, ">", "bins.tsv";
say $bin_fh join "\t", "chr", "start", "end";
for my $chr (@chromosomes) {
    my @lengths;
    for my $start ( sort { $a <=> $b } keys %{$bins{$chr}} ) {
        my $end = $bins{$chr}->{$start};
        push @lengths, $end - $start + 1;
        say $bin_fh join "\t", $chr, $start, $end;
    }
    summarize_ranges( $chr, \@lengths );
}
close $bin_fh;

sub summarize_ranges {
    my ( $chr, $lengths_ref ) = @_;

    my $count = @$lengths_ref;
    my $min   = min @$lengths_ref;
    my $mean  = sprintf( '%.0f', sum(@$lengths_ref) / $count );
    my $max   = max @$lengths_ref;

    say "$chr";
    say "  bins: $count";
    say "  min:  $min";
    say "  mean: $mean";
    say "  max:  $max";
}

exit;
