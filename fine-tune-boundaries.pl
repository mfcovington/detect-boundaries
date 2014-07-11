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
use File::Path 'make_path';
use List::Util qw(min max);
use List::MoreUtils 'first_index';
use Scalar::Util 'looks_like_number';
use Term::ANSIColor;

use Data::Printer;

# TODO: Change temporary fix for chromosome name mismatch

my $genotyped_dir;
my $context = 10;
my $out_dir = ".";

my $options = GetOptions(
    "genotyped_dir=s" => \$genotyped_dir,
    "context=i"       => \$context,
    "out_dir=s"       => \$out_dir,
);

my @boundaries_files = @ARGV;

make_path $out_dir;

for my $bounds_file (@boundaries_files) {
    my $sample_id = get_sample_id($bounds_file);

    my $bounds_out_file = "$out_dir/$sample_id.boundaries";
    die
        "ERROR: Output file ($bounds_out_file) already exists, choose a clean output directory.\n"
        if -e $bounds_out_file;

    my $boundaries = get_boundaries($bounds_file);
    my $genotypes = get_genotypes( $genotyped_dir, $sample_id );
    compare_chromosome_counts( $boundaries, $genotypes, $sample_id );

    my %corrected_boundaries;
    for my $chr ( sort keys %$boundaries ) {
        $corrected_boundaries{$chr} = analyze_bins_for_chr( $sample_id, $chr, $$boundaries{$chr},
            $$genotypes{"SL2.40$chr"} ); #temporary fix for chromosome name mismatch
    }

    my $warnings = validate_boundaries( \%corrected_boundaries, $genotypes );

    if ( scalar @$warnings > 0 ) {
        issue_warnings( $warnings, $sample_id );
    }
    else {
        write_boundaries( \%corrected_boundaries, $bounds_out_file );
    }
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

sub compare_chromosome_counts {
    my ( $boundaries, $genotypes, $sample_id ) = @_;
    my @b_chromosomes = sort keys %$boundaries;
    my @g_chromosomes = sort keys %$genotypes;
    die
        "$sample_id: Inconsistent chromosome count for boundaries vs genotypes!\n"
        if scalar @b_chromosomes != scalar @g_chromosomes;
}

sub analyze_bins_for_chr {
    my ( $sample_id, $chr, $bins, $geno_scores ) = @_;
    my $bin_count = scalar @$bins;

    my @geno_positions = sort { $a <=> $b } keys %$geno_scores;

    my $old_pos  = $geno_positions[0];
    my $old_geno = 'START OF CHROMOSOME';

    my $previous_start;
    my %corrected_bins;
    my $corrected_breakpoints;
    for my $current_bin (@$bins) {
        my $new_pos  = $$current_bin{'start'};
        my $new_geno = $$current_bin{'geno'};

        display_breakpoint( $old_geno, $old_pos, $new_geno, $new_pos,
            \@geno_positions, $geno_scores, $sample_id, $chr );

        if ( !defined $previous_start ) {
            $corrected_breakpoints = is_breakpoint_good($new_geno);
        }
        else {
            $corrected_breakpoints = is_breakpoint_good( $old_geno, $new_geno );

            if ( exists $$corrected_breakpoints{$old_geno} ) {
                $corrected_bins{$previous_start}{'end'} = $$corrected_breakpoints{$old_geno};
            }
            else {
                $corrected_bins{$previous_start}{'end'} = $old_pos;
            }
        }

        $new_pos = $$corrected_breakpoints{$new_geno}
            if exists $$corrected_breakpoints{$new_geno};
        $corrected_bins{$new_pos}{'geno'} = $new_geno;

        $previous_start = $new_pos;

        $old_pos  = $$current_bin{'end'};
        $old_geno = $new_geno;
    }

    display_breakpoint( $old_geno, $old_pos, 'END OF CHROMOSOME',
        $geno_positions[-1], \@geno_positions, $geno_scores, $sample_id,
        $chr );

    $corrected_breakpoints = is_breakpoint_good( $old_geno );

    if ( exists $$corrected_breakpoints{$old_geno} ) {
        $corrected_bins{$previous_start}{'end'} = $$corrected_breakpoints{$old_geno};
    }
    else {
        $corrected_bins{$previous_start}{'end'} = $old_pos;
    }

    return \%corrected_bins;
}

sub display_breakpoint {
    my ( $old_geno, $old_pos, $new_geno, $new_pos, $geno_positions,
        $geno_scores, $sample_id, $chr )
        = @_;

    my $error_msg
        = "ERROR: Overlapping bins for $sample_id on $chr ($old_pos should be less than $new_pos)\n";
    my $terminal_merge = 0;
    if ( $old_pos == $new_pos ) {

        if (   $old_geno eq 'START OF CHROMOSOME'
            || $new_geno eq 'END OF CHROMOSOME' )
        {
            $terminal_merge++;
        }
        else {
            die $error_msg;
        }
    }
    elsif ( $old_pos > $new_pos ) {
        die $error_msg;
    }

    my $old_idx = first_index { $_ == $old_pos } @$geno_positions;
    my $new_idx = first_index { $_ == $new_pos } @$geno_positions;

    my $pre = max( $old_idx - $context, 0 );
    my $post = min( $new_idx + $context, $#$geno_positions );

    if ( $old_idx > $pre ) {
        say $$geno_scores{$_} for @$geno_positions[ $pre .. $old_idx - 1 ];
    }

    if ($terminal_merge) {
        say colored ['bright_white on_magenta'],
            "$$geno_scores{$$geno_positions[$old_idx]}\t$old_geno --- $new_geno ";
    }
    else {
        say colored ['bright_white on_red'],
            "$$geno_scores{$$geno_positions[$old_idx]}\t$old_geno ";
        say $$geno_scores{$_}
            for @$geno_positions[ $old_idx + 1 .. $new_idx - 1 ];
        say colored ['bright_white on_blue'],
            "$$geno_scores{$$geno_positions[$new_idx]}\t$new_geno ";
    }

    if ( $new_idx < $post ) {
        say $$geno_scores{$_} for @$geno_positions[ $new_idx + 1 .. $post ];
    }
}

sub is_breakpoint_good {
    my @genotypes = @_;

    my $yes_no;
    my $input_valid = 0;
    while ( !$input_valid ) {
        print colored ['bold bright_cyan on_black'],
            "Does this breakpoint look good? (y/n) ";
        $yes_no = <STDIN>;
        $input_valid++ if $yes_no =~ /^[yn]$/i;
    }
    return if $yes_no =~ /^y$/i;

    my %corrected_breakpoints;
    enter_new_breakpoint( $_, \%corrected_breakpoints ) for @genotypes;
    return \%corrected_breakpoints;
}

sub enter_new_breakpoint {
    my ( $genotype, $corrected_breakpoints ) = @_;

    my $position;
    my $input_valid = 0;
    while ( !$input_valid ) {
        print colored ['bold bright_cyan on_black'],
            "Enter new position for $genotype (leave blank if OK): ";
        chomp( $position = <STDIN> );
        $input_valid++ if looks_like_number $position || $position eq '';
    }
    $$corrected_breakpoints{$genotype} = $position unless $position eq '';
}

sub validate_boundaries {
    my ( $corrected_boundaries, $genotypes ) = @_;

    my @warnings;
    for my $chr ( sort keys %$corrected_boundaries ) {

        my $previous_geno = "";
        my $previous_end  = 0;

        my @geno_positions
            = sort { $a <=> $b } keys $$genotypes{"SL2.40$chr"}; #temporary fix for chromosome name mismatch
        my $first_pos = $geno_positions[0];
        my $last_pos  = $geno_positions[-1];

        for my $start ( sort { $a <=> $b } keys $$corrected_boundaries{$chr} )
        {
            my $end = $$corrected_boundaries{$chr}{$start}{'end'};
            my $current_genotype
                = $$corrected_boundaries{$chr}{$start}{'geno'};

            push @warnings,
                "Current bin start ($chr:$start) <= previous bin end ($chr:$previous_end)"
                if $start <= $previous_end;
            push @warnings, "Bin start ($chr:$start) >= bin end ($chr:$end)"
                if $start >= $end;
            push @warnings,
                "Bin start ($chr:$start) < first genotyped position ($first_pos)"
                if $start < $first_pos;
            push @warnings,
                "Bin end ($chr:$end) > last genotyped position ($last_pos)"
                if $end > $last_pos;
            push @warnings, "Genotypes don't alternate on $chr"
                if $current_genotype eq $previous_geno;

            $previous_end = $end;
        }
    }

    return \@warnings;
}

sub issue_warnings {
    my ( $warnings, $sample_id ) = @_;

    say colored [ 'bright_red on_black' ], "ERROR: Boundaries file not written for $sample_id";
    say colored [ 'bright_red on_black' ], "  $_" for @$warnings;
}

sub write_boundaries {
    my ( $boundaries, $bounds_out_file ) = @_;

    open my $bounds_out_fh, ">", $bounds_out_file;
    for my $chr ( sort keys %$boundaries ) {
        for my $start ( sort { $a <=> $b } keys $$boundaries{$chr} ) {
            say $bounds_out_fh join "\t", $chr, $start,
                $$boundaries{$chr}{$start}{'end'},
                $$boundaries{$chr}{$start}{'geno'};
        }
    }
    close $bounds_out_fh;
}
