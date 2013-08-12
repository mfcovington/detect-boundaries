#!/usr/bin/env perl
# get-genotypes-for-subset.pl
# Mike Covington
# created: 2013-08-09
#
# Description:
#
use strict;
use warnings;
use autodie;
use feature 'say';

open my $bins_fh, "<", "sample-file/bins-snp.mid";
<$bins_fh>;
my %mids;
while (<$bins_fh>) {
    chomp;
    my ( $chr, $pos ) = (split)[0, 3];
    next if $pos eq 'NA';
    $mids{$chr}{$pos} = ();
}

# ls sample-file | grep snps | grep A01 | cut -d. -f1
my @ids = qw(RIL_1 RIL_103 RIL_104 RIL_113 RIL_115 RIL_12 RIL_123 RIL_124 RIL_131 RIL_136 RIL_143 RIL_146 RIL_147 RIL_15 RIL_150 RIL_154 RIL_155 RIL_16 RIL_164 RIL_166 RIL_171 RIL_174 RIL_175 RIL_176 RIL_182 RIL_183 RIL_184 RIL_187 RIL_190 RIL_193 RIL_198 RIL_199 RIL_2 RIL_201 RIL_204 RIL_205 RIL_206 RIL_207 RIL_208 RIL_21 RIL_211 RIL_212 RIL_213 RIL_215 RIL_222 RIL_223 RIL_225 RIL_228 RIL_229 RIL_23 RIL_232 RIL_234 RIL_235 RIL_240 RIL_242 RIL_243 RIL_248 RIL_25 RIL_250 RIL_251 RIL_253 RIL_255 RIL_256 RIL_259 RIL_264 RIL_265 RIL_267 RIL_268 RIL_270 RIL_277 RIL_281 RIL_282 RIL_284 RIL_285 RIL_288 RIL_289 RIL_290 RIL_294 RIL_30 RIL_300 RIL_301 RIL_303 RIL_308 RIL_31 RIL_311 RIL_318 RIL_325 RIL_329 RIL_332 RIL_337 RIL_339 RIL_340 RIL_341 RIL_344 RIL_346 RIL_347 RIL_353 RIL_354 RIL_355 RIL_357 RIL_359 RIL_36 RIL_360 RIL_363 RIL_373 RIL_376 RIL_380 RIL_39 RIL_42 RIL_46 RIL_53 RIL_60 RIL_63 RIL_65 RIL_66 RIL_69 RIL_7 RIL_70 RIL_73 RIL_76 RIL_80 RIL_89 RIL_9 RIL_93);
@ids = sort @ids;
my @chromosomes = qw(A01 A02 A03 A04 A05 A06 A07 A08 A09 A10);
for my $id (@ids) {
    for my $chr (@chromosomes) {
        open my $snps_fh, "<", "sample-file/$id.$chr.filtered.snps";
        my %snps;
        while (<$snps_fh>) {
            chomp;
            my ( $pos, $genotype ) = (split)[1, 3];
            $snps{$chr}{$pos} = $genotype;
        }
        close $snps_fh;

        for my $pos (sort { $a <=> $b } keys $mids{$chr} ) {
            $mids{$chr}{$pos}{$id} = $snps{$chr}{$pos} // 'NA';
        }
    }
}

my %polydb;
for my $chr (@chromosomes) {
    open my $polydb_fh, "<", "sample-file/polyDB.$chr";
    <$polydb_fh>;
    $polydb{$chr} = { map { (split)[1] => 1 } <$polydb_fh> };
    close $polydb_fh;
}

open my $genotype_fh, ">", "bin-genotype";
say $genotype_fh join "\t", "chr", "pos", @ids;
for my $chr ( sort keys %mids ) {
    for my $pos ( sort {$a <=> $b} keys $mids{$chr} ) {
        next unless exists $polydb{$chr}{$pos};
        my $line = "$chr\t$pos";
        for my $id ( sort keys $mids{$chr}{$pos} ) {
            $line .= "\t$mids{$chr}{$pos}{$id}";
        }
        say $genotype_fh $line;
    }
}
close $genotype_fh;

