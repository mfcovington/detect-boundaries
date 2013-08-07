#!/usr/bin/env perl
# FILE_NAME.pl
# Mike Covington
# created: 2013-06-28
#
# Description:
#
use strict;
use warnings;
use autodie;
use feature 'say';
use Data::Printer;
use List::Util qw(min max sum);

# get rid of lonely snps
my $id              = shift;
my @genotyped_files = @ARGV;

my $par1_id = "R500";
my $par2_id = "IMB211";

my $min_cov      = 10;
my $min_momentum = 10;
my $min_ratio    = 0.9;
my $het_offset   = 0.2;

my $het_max = min( $min_ratio, 0.5 + $het_offset );

my %snps;
for my $file (@genotyped_files) {
    my @recent;
    my $momentum    = $min_momentum;
    my $last_parent = '';
    my $monitor     = 1;
    my $chr;

    open my $snp_fh, "<", $file;
    while (<$snp_fh>) {
        ( $chr, my ( $pos, $par1, $par2, $tot ) ) = split /\t/;
        next if $par1 + $par2 < $min_cov;

        my $ratio = max( $par1, $par2 ) / ( $par1 + $par2 );
        my $cur_parent;
        if ( $ratio >= $min_ratio ) {
            $cur_parent = $par1 > $par2 ? $par1_id : $par2_id;
            $snps{$chr}{$pos} = {
                score  => max( $par1, $par2 ),
                par_id => $cur_parent
            };
        }
        elsif ( $ratio < $het_max ) {
            $cur_parent = 'HET';
            $snps{$chr}{$pos} = {
                score  => $par1 + $par2,
                par_id => $cur_parent,
                hetrat => $par1 / ($par1 + $par2)
            };
        }
        else { next }

        if ( $cur_parent ne $last_parent ) {
            if ( $momentum < $min_momentum ) {
                delete $snps{$chr}{$_} for @recent;
                @recent = ();
            }
            $momentum = 0;
        }
        else {
            $momentum++;
        }

        push @recent, $pos;
        @recent = () if $momentum >= $min_momentum;
        $last_parent = $cur_parent;
    }
    close $snp_fh;
    delete $snps{$chr}{$_} for @recent;
}

open my $boundaries_out_fh, ">", "$id.boundaries";
for my $chr ( sort keys %snps ) {
    my $last_parent = '';
    my $last_pos    = 0;

    open my $snp_out_fh, ">", "$id.$chr.filtered.snps";
    for my $pos ( sort { $a <=> $b } keys $snps{$chr} ) {
        my $score  = $snps{$chr}{$pos}{score};
        my $par_id = $snps{$chr}{$pos}{par_id};

        # output snps
        say $snp_out_fh join "\t", $chr, $pos, $score, $par_id;

        # output boundaries
        if ( $last_parent ne $par_id ) {
            say $boundaries_out_fh join "\t", $last_pos, $last_parent
              unless $last_pos == 0;
            print $boundaries_out_fh join "\t", $chr, $pos, "";
            $last_parent = $par_id;
        }
        $last_pos = $pos;
    }
    say $boundaries_out_fh join "\t", $last_pos, $last_parent;

    close $snp_out_fh;
}
close $boundaries_out_fh;
