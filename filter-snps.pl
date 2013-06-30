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
use List::Util qw(min max);

# get rid of lonely snps
my $id              = shift;
my @genotyped_files = @ARGV;

my $par1_id = "R500";
my $par2_id = "IMB211";

my $min_cov      = 3;
my $min_momentum = 5;
my $min_ratio    = 0.9;

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
        next if max( $par1, $par2 ) / ( $par1 + $par2 ) < $min_ratio;

        my $cur_parent = $par1 > $par2 ? $par1_id : $par2_id;

        $snps{$chr}{$pos} = {
            score  => max( $par1, $par2 ),
            par_id => $cur_parent
        };

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

for my $chr ( sort keys %snps ) {
    open my $out_fh, ">", "$id.$chr.filtered.snps";
    for my $pos ( sort { $a <=> $b } keys $snps{$chr} ) {
        my $score  = $snps{$chr}{$pos}{score};
        my $par_id = $snps{$chr}{$pos}{par_id};
        say $out_fh join "\t", $chr, $pos, $score, $par_id;
    }
    close $out_fh;
}
