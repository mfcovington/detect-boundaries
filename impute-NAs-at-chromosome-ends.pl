#!/usr/bin/env perl
# Mike Covington
# created: 2014-08-12
#
# Description:
#
use strict;
use warnings;
use Log::Reproducible;
use autodie;
use feature 'say';
use List::MoreUtils qw(firstidx lastidx);

my $bin_geno_file = $ARGV[0];  # e.g., "test-NAs/bin-genotype.BILs.2014-08-11"

die "Usage: $0 bin-genotype-file\n" unless defined $bin_geno_file;

open my $bin_geno_fh, "<", $bin_geno_file;
my $header_line   = <$bin_geno_fh>;
my $sample_ids    = get_sample_ids($header_line);
my $bin_genotypes = get_bin_genotypes( $bin_geno_fh, $sample_ids );
close $bin_geno_fh;

impute_genotypes( $bin_genotypes, $sample_ids );
output_bin_genotypes( "$bin_geno_file.imputed-NAs", $bin_genotypes,
    $header_line );

merge_adjacent_like_bins($bin_genotypes);
output_bin_genotypes( "$bin_geno_file.imputed-NAs.merged-like",
    $bin_genotypes, $header_line );

exit;


sub get_sample_ids {
    my $line = shift;
    my @elements = split /\t/, $line;
    chomp @elements;
    return [ @elements[ 4 .. $#elements ] ];
}

sub get_bin_genotypes {
    my ( $bin_geno_fh, $sample_ids ) = @_;

    my %bin_genotypes;
    while (<$bin_geno_fh>) {
        chomp;
        my ( $chr, $mid, $start, $end, @genotypes ) = split;
        $bin_genotypes{$chr}{$mid}
            = { 'start' => $start, 'end' => $end, 'merge' => 0 };
        for my $i ( 0 .. $#$sample_ids ) {
            $bin_genotypes{$chr}{$mid}{'genotypes'}{ $$sample_ids[$i] }
                = $genotypes[$i];
        }
    }

    return \%bin_genotypes;
}

sub impute_genotypes {
    my ( $bin_genotypes, $sample_ids ) = @_;

    for my $sample (@$sample_ids) {
        for my $chr ( sort keys %$bin_genotypes ) {

            my @genotypes;
            for my $bin_mid ( sort { $a <=> $b } keys %{ $$bin_genotypes{$chr} } ) {
                push @genotypes,
                    $$bin_genotypes{$chr}{$bin_mid}{'genotypes'}{$sample};
            }

            if ( $genotypes[0] eq 'NA' ) {
                my $first_geno_idx = firstidx { $_ ne 'NA' } @genotypes;
                $genotypes[$_] = $genotypes[$first_geno_idx]
                    for 0 .. $first_geno_idx - 1;
            }

            if ( $genotypes[-1] eq 'NA' ) {
                my $last_geno_idx  = lastidx  { $_ ne 'NA' } @genotypes;
                $genotypes[$_] = $genotypes[$last_geno_idx]
                    for $last_geno_idx + 1 .. $#genotypes;
            }

            for my $bin_mid ( sort { $a <=> $b } keys %{ $$bin_genotypes{$chr} } ) {
                $$bin_genotypes{$chr}{$bin_mid}{'genotypes'}{$sample}
                    = shift @genotypes;
            }
        }
    }
}

sub merge_adjacent_like_bins {
    my $bin_genotypes = shift;

    for my $chr ( sort keys %$bin_genotypes ) {
        my $previous_geno = '';
        my $previous_pos  = '';
        for my $bin_mid ( sort { $a <=> $b } keys %{ $$bin_genotypes{$chr} } ) {
            my $genotypes = join ",",
                extract_genotypes_for_position( $chr, $bin_mid,
                $bin_genotypes );

            if ( $genotypes eq $previous_geno ) {
                $$bin_genotypes{$chr}{$previous_pos}{'merge'} = 1;
            }

            $previous_geno = $genotypes;
            $previous_pos  = $bin_mid;
        }
    }
}

sub output_bin_genotypes {
    my ( $out_file, $bin_genotypes, $header_line ) = @_;

    open my $out_fh, ">", $out_file;
    print $out_fh $header_line;
    for my $chr ( sort keys %$bin_genotypes ) {
        my @merged_start;
        for my $bin_mid ( sort { $a <=> $b } keys %{ $$bin_genotypes{$chr} } ) {
            my $start = $$bin_genotypes{$chr}{$bin_mid}{'start'};

            if ( $$bin_genotypes{$chr}{$bin_mid}{'merge'} ) {
                push @merged_start, $start;
                next;
            }

            my $end = $$bin_genotypes{$chr}{$bin_mid}{'end'};
            my @genotypes = extract_genotypes_for_position( $chr, $bin_mid,
                $bin_genotypes );

            if ( scalar @merged_start > 0 ) {
                $start = $merged_start[0];
                $bin_mid = int( ( $start + $end ) / 2 );
            }

            say $out_fh join "\t", $chr, $bin_mid, $start, $end, @genotypes;

            @merged_start = ();
        }
    }
    close $out_fh;
}

sub extract_genotypes_for_position {
    my ( $chr, $bin_mid, $bin_genotypes ) = @_;

    return map { $$bin_genotypes{$chr}{$bin_mid}{'genotypes'}{$_} }
        sort keys %{ $$bin_genotypes{$chr}{$bin_mid}{'genotypes'} };
}
