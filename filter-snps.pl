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
# my $file = $ARGV[0];
my $file = 'original/RIL_1.12.A02.genotyped';

my $par1_id = "R500";
my $par2_id = "IMB211";

my $min_cov      = 3;
my $min_momentum = 5;

open my $snp_fh, "<", $file;

my %snps;
my @recent;
my $momentum = $min_momentum;
my $count;
my $last_parent = '';
my $monitor     = 1;
while (<$snp_fh>) {
    my ( $chr, $pos, $par1, $par2, $tot ) = split /\t/;
    next if $par1 && $par2;
    next if $par1 + $par2 < $min_cov;

    my $cur_parent = $par1 > $par2 ? $par1_id : $par2_id;

    $snps{$pos} = {
        chr    => $chr,
        score  => max( $par1, $par2),
        par_id => $cur_parent
    };
    # $snps{$pos} = $cur_parent;

    if ( $cur_parent ne $last_parent ) {
        if ( $momentum < $min_momentum ) {
            delete $snps{$_} for @recent;
            @recent = ();
        }
        $momentum = 0;
    }
    else {
        $momentum++; # = min( $min_momentum, ++$momentum );
        # $momentum = min( $min_momentum, ++$momentum );
    }
    push @recent, $pos;
    @recent = () if $momentum >= $min_momentum;
    # @recent = () if $momentum == $min_momentum;
    say "$pos:$cur_parent:$momentum";
    $last_parent = $cur_parent;
    $count++;
}
close $snp_fh;
delete $snps{$_} for @recent;

# p %snps;
say scalar keys %snps;
say $count;

open my $out_fh, ">", "snps.out";
for ( sort { $a <=> $b } keys %snps ) {
    say $out_fh join "\t", $snps{$_}{chr}, $_, $snps{$_}{score}, $snps{$_}{par_id};
}
close $out_fh;
