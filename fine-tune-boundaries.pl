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
use Term::ANSIScreen 'cls';
use Term::ReadKey;
use v5.10.1;

use Data::Printer;

# TODO: Change temporary fix for chromosome name mismatch

my $sample_id;

sub safe_exit {
    ReadMode 0;
    say "\n", colored ['bright_blue on_bright_yellow'],
        "  * You stopped while working on $sample_id *  ";
    exit;
}
$SIG{INT} = \&safe_exit;

my ( $genotyped_dir, $seq_list, $par1, $par2 );
my $context = 10;
my $out_dir = ".";

my $options = GetOptions(
    "genotyped_dir=s" => \$genotyped_dir,
    "context=i"       => \$context,
    "out_dir=s"       => \$out_dir,
    "seq_list=s"      => \$seq_list,
    "par1=s"          => \$par1,
    "par2=s"          => \$par2,
);

my @boundaries_files = @ARGV;

my $sequences;
$sequences = { map { $_ => 1 } split /,/, $seq_list }
    if defined $seq_list;

my $parents;
$parents = { '1' => $par1, '2' => $par2, } if defined $par1 && defined $par2;

make_path $out_dir;

my $counter = 1;
my $total = scalar @boundaries_files;
SAMPLE: for my $bounds_file (@boundaries_files) {
    $sample_id = get_sample_id($bounds_file);
    my $bounds_out_file = "$out_dir/$sample_id.boundaries";
    die
        "ERROR: Output file ($bounds_out_file) already exists, choose a clean output directory.\n"
        if -e $bounds_out_file;

    my $boundaries = get_boundaries( $bounds_file, $sequences );
    my $genotypes = get_genotypes( $genotyped_dir, $sample_id, $sequences );
    compare_chromosome_counts( $boundaries, $genotypes, $sample_id );

    my %corrected_boundaries;
    for my $chr ( sort keys %$boundaries ) {
        my $redo_sample = 0;
        $corrected_boundaries{$chr} = analyze_bins_for_chr( $sample_id, $chr, $$boundaries{$chr},
            $$genotypes{$chr}, \$redo_sample, $parents ); #temporary fix for chromosome name mismatch
        redo SAMPLE if $redo_sample;
    }

    my $warnings = validate_boundaries( \%corrected_boundaries, $genotypes );

    if ( scalar @$warnings > 0 ) {
        issue_warnings( $warnings, $sample_id );
        redo if redo_or_continue();
    }
    else {
        write_boundaries( \%corrected_boundaries, $bounds_out_file );
        motivate();
    }
    $counter++;
}

sub get_sample_id {
    my $bounds_file = shift;
    my ($sample_id) = fileparse( $bounds_file, ".boundaries" );
    return $sample_id;
}

sub get_boundaries {
    my ( $bounds_file, $sequences ) = @_;
    my %boundaries;
    open my $bounds_fh, "<", $bounds_file;
    while (<$bounds_fh>) {
        chomp;
        my ( $chr, $start, $end, $geno ) = split;
        if ( defined $sequences ) {
            next unless exists $$sequences{$chr};
        }
        push @{ $boundaries{$chr} },
            { 'start' => $start, 'end' => $end, 'geno' => $geno };
    }
    close $bounds_fh;

    if ( defined $sequences ) {
        for my $chr ( sort keys %$sequences ) {
            push @{ $boundaries{$chr} },
                { 'start' => 0, 'end' => 0, 'geno' => 'NA' }
                unless defined $boundaries{$chr};
        }
    }

    return \%boundaries;
}

sub get_genotypes {
    my ( $genotyped_dir, $sample_id, $sequences ) = @_;
    my %genotypes;
    my @genotyped_files = glob "$genotyped_dir/$sample_id.*.genotyped*";
    for my $geno_file (@genotyped_files) {
        open my $geno_fh, "<", $geno_file;
        while ( my $line = <$geno_fh> ) {
            chomp $line;
            my ( $chr, $pos ) = split /\t/, $line;
            if ( defined $sequences ) {
                next unless exists $$sequences{$chr};
            }
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
    my ( $sample_id, $chr, $bins, $geno_scores, $redo_sample, $parents ) = @_;
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
            $corrected_breakpoints
                = is_breakpoint_good( $redo_sample, $parents, $old_pos, $new_geno );
            ($new_geno) = keys %$corrected_breakpoints if $new_geno eq 'NA';
        }
        else {
            $corrected_breakpoints
                = is_breakpoint_good( $redo_sample, $parents, undef, $old_geno, $new_geno );
            ($new_geno) = keys %$corrected_breakpoints if $new_geno eq 'NA';

            if ( exists $$corrected_breakpoints{$old_geno} ) {
                $corrected_bins{$previous_start}{'end'} = $$corrected_breakpoints{$old_geno};
            }
            else {
                $corrected_bins{$previous_start}{'end'} = $old_pos;
            }
        }
        return if $$redo_sample;

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

    $corrected_breakpoints
        = is_breakpoint_good( $redo_sample, $parents, $geno_positions[-1],
        $old_geno // 'NA' );

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
    $old_geno //= 'NA';

    print cls();
    print "\n";
    say colored ['bright_blue on_bright_yellow'],
        "  * Processing $sample_id ($counter/$total) *  ";
    print "\n";

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
            die $error_msg unless $new_geno eq 'NA';
        }
    }
    elsif ( $old_pos > $new_pos ) {
        die $error_msg unless $new_geno eq 'NA';
    }

    my $old_idx = first_index { $_ == $old_pos } @$geno_positions;
    my $new_idx = first_index { $_ == $new_pos } @$geno_positions;

    my $pre = max( $old_idx - $context, 0 );
    my $post = min( $new_idx + $context, $#$geno_positions );

    if ( $old_idx > $pre ) {
        highlight_zeroes( $$geno_scores{$_} )
            for @$geno_positions[ $pre .. $old_idx - 1 ];
    }

    if ($terminal_merge) {
        say colored ['bright_white on_magenta'],
            "$$geno_scores{$$geno_positions[$old_idx]}\t$old_geno --- $new_geno ";
    }
    elsif ( $new_idx == -1 ) {
        say colored ['bright_white on_red'],
            "$$geno_scores{$$geno_positions[$old_idx]}\t$old_geno ";
        say colored ['bright_white on_blue'],
            "$new_geno\t$new_geno\t$new_geno\t$new_geno\t$new_geno ";
        highlight_zeroes( $$geno_scores{$_} )
            for @$geno_positions[ $old_idx + 1 .. $post ];
    }
    elsif ( $old_idx == -1 ) {
        highlight_zeroes( $$geno_scores{$_} )
            for @$geno_positions[ $new_idx - 9 .. $new_idx - 1 ];
        say colored ['bright_white on_red'],
            "$old_geno\t$old_geno\t$old_geno\t$old_geno\t$old_geno ";
        say colored ['bright_white on_blue'],
            "$$geno_scores{$$geno_positions[$new_idx]}\t$new_geno ";
    }
    else {
        say colored ['bright_white on_red'],
            "$$geno_scores{$$geno_positions[$old_idx]}\t$old_geno ";
        highlight_zeroes( $$geno_scores{$_} )
            for @$geno_positions[ $old_idx + 1 .. $new_idx - 1 ];
        say colored ['bright_white on_blue'],
            "$$geno_scores{$$geno_positions[$new_idx]}\t$new_geno ";
    }

    if ( $new_idx >= 0 && $new_idx < $post ) {
        highlight_zeroes( $$geno_scores{$_} )
            for @$geno_positions[ $new_idx + 1 .. $post ];
    }
}

sub highlight_zeroes {
    my $geno_line = shift;

    say join "\t",
        map { $_ = $_ eq 0 ? colored( 0, 'black on_bright_yellow' ) : $_ }
        split /\t/,
        $geno_line;
}

sub is_breakpoint_good {
    my ( $redo_sample, $parents, $default_pos, @genotypes ) = @_;

    my $yes_no;
    my $input_valid = 0;
    while ( !$input_valid ) {
        print colored ['bold bright_cyan on_black'],
            "\nDoes this breakpoint look good? (y/n/r/p/x/?) ";
        ReadMode 3;
        while ( not defined( $yes_no = ReadKey(-1) ) ) { }
        ReadMode 0;

        if    ( $yes_no =~ /^[yn]$/i ) { $input_valid++ }
        elsif ( $yes_no =~ /^r$/i )    { $$redo_sample++; return }
        elsif ( $yes_no =~ /^p$/i )    { take_a_break() }
        elsif ( $yes_no =~ /^x$/i )    { safe_exit() }
        elsif ( $yes_no =~ /^\?$/i )   { help() }
    }
    return if $yes_no =~ /^y$/i;

    print "\n";
    my %corrected_breakpoints;
    enter_new_breakpoint( $_, $parents, $default_pos, \%corrected_breakpoints )
        for @genotypes;
    return \%corrected_breakpoints;
}

sub take_a_break {
    my $pause_msg = <<EOF;

Analysis paused. This eliminates unnecessary CPU usage. Enjoy your break!
Press enter to continue.
EOF
    chomp $pause_msg;
    print colored ['bold bright_cyan on_black'], $pause_msg;
    <STDIN>;
}

sub help {
    my $help_message = <<EOF;

Y: Yes, breakpoint looks good.
N: No, breakpoint needs to be adjusted.
R: Make a mistake? Redo the current sample.
P: Pause analysis. This eliminates unnecessary CPU usage.
X: Exit.
?: You are here.
EOF
    chomp $help_message;
    say colored ['bright_blue'], $help_message;
}

sub enter_new_breakpoint {
    my ( $genotype, $parents, $default_pos, $corrected_breakpoints ) = @_;

    my $genotype_response;
    if ( $genotype eq 'NA' ) {
        print colored ['bold bright_cyan on_black'],
            "Enter correct genotype (leave blank if OK): ";
        chomp( $genotype_response = <STDIN> );
        $genotype = $genotype_response eq '' ? $genotype : $genotype_response;

        $genotype = $$parents{$genotype}
            if defined $parents && exists $$parents{$genotype};
    }

    my $position;
    my $input_valid = 0;
    while ( !$input_valid ) {
        print colored ['bold bright_cyan on_black'],
            "Enter new position for $genotype (leave blank if OK): ";
        chomp( $position = <STDIN> );
        $position = $default_pos
            if defined $default_pos && $position =~ /^d$/i;
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
            = sort { $a <=> $b } keys %{ $$genotypes{$chr} }; #temporary fix for chromosome name mismatch
        my $first_pos = $geno_positions[0];
        my $last_pos  = $geno_positions[-1];

        for my $start ( sort { $a <=> $b }
            keys %{ $$corrected_boundaries{$chr} } )
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

sub redo_or_continue {
    my $response;
    my $input_valid = 0;
    while ( !$input_valid ) {
        print colored ['bold bright_cyan on_black'],
            "\nDo you want to (r)edo this sample or (c)ontinue to the next? (r/c) ";
        ReadMode 3;
        while ( not defined( $response = ReadKey(-1) ) ) { }
        ReadMode 0;
        $input_valid++ if $response =~ /^[rc]$/i;
    }
    return 1 if $response =~ /^r$/i;
}

sub write_boundaries {
    my ( $boundaries, $bounds_out_file ) = @_;

    open my $bounds_out_fh, ">", $bounds_out_file;
    for my $chr ( sort keys %$boundaries ) {
        for my $start ( sort { $a <=> $b } keys %{ $$boundaries{$chr} } ) {
            say $bounds_out_fh join "\t", $chr, $start,
                $$boundaries{$chr}{$start}{'end'},
                $$boundaries{$chr}{$start}{'geno'};
        }
    }
    close $bounds_out_fh;
}

sub motivate {
    my @colors = qw(clear bright_red bright_green bright_yellow bright_blue
        bright_magenta bright_cyan);
    my $motivation = <<'EOF';

   ______                __       __      __    __ __ __
  / ____/___  ____  ____/ /      / /___  / /_  / // // /
 / / __/ __ \/ __ \/ __  /  __  / / __ \/ __ \/ // // /
/ /_/ / /_/ / /_/ / /_/ /  / /_/ / /_/ / /_/ /_//_//_/
\____/\____/\____/\__,_/   \____/\____/_.___(_)(_)(_)

EOF

    print cls();
    print colored [ $colors[ rand @colors ] ], $motivation;
}
