#!/usr/bin/perl

################################
#Compute all promoter positions#
#and promoter sequences        #
#of a contig.                  #
################################

use List::Util ("max");

#assumed start and end of promoters related to transcription-start
my $from = 1000;
my $to   = 50;

#contig of interest
my $chromosome = 7;
my $contig     = "supercont2.7";

#parameters of annotation file
my $tabsBeforeLocus  = 0;
my $tabsBeforeStart  = 4;
my $tabsBeforeStop   = 5;
my $tabsBeforeStrand = 6;

#information to extract from annotation file
my @loci;
my @starts;
my @stops;
my @strands;

#information to extract from genome file
my $sequence;
my $nucleotides;

#open genome file
open( my $genome, "<", "aspergillus_flavus_2_supercontigs.fasta" ) || die "genome input error: " . $!;
my $found = 0;

#we are not at the end of the file and
#not at a line after the selected contig
while ( !eof($genome) && !$found ) {
	my $line = <$genome>;

	#found the header of the selected contig
	if ( index( $line, ">" ) != -1 && index( $line, $contig ) != -1 ) {
		$found = 1;
		$line  = <$genome>;

		#read the genome line by line
		#until we get to the next contig > next header
		while ( !eof && index( $line, ">" ) == -1 ) {
			$sequence .= $line;
			chomp($sequence);
			$line = <$genome>;
		}
	}
}
close $genome;
$sequence    = lc($sequence);
$nucleotides = length($sequence);
print "Contig " . $contig . " has " . $nucleotides . " nucleotides";

#open annotation file
open( my $annotation, "<", "aspergillus_flavus_2_genome_summary_per_gene.txt" )
  || die "annotation input error: " . $!;

#analyse each line of the annotation file
<$annotation>;    #skip first line (header)
while (<$annotation>) {
	my $line = $_;

	#check if the current line belongs to the selected contig
	if ( index( $line, "\t" . $chromosome . "\t" ) != -1 ) {
		my $locus;
		my $startPos;
		my $stopPos;
		my $strandPos;

		#extract locus, start, stop and strand
		my $noOfTab = 1;
		my $tabPos  = 0;
		do {
			$tabPos = index( $line, "\t", $tabPos + 1 );
			if ( $noOfTab == $tabsBeforeLocus + 1 ) {
				$locus = substr( $line, 0, $tabPos );
			}
			elsif ( $noOfTab == $tabsBeforeStart ) {
				$startPos = $tabPos + 1;
			}
			elsif ( $noOfTab == $tabsBeforeStop ) {
				$stopPos = $tabPos + 1;
			}
			elsif ( $noOfTab == $tabsBeforeStrand ) {
				$strandPos = $tabPos + 1;
			}
			++$noOfTab;
		} while ( $tabPos != -1 );
		my $start  = substr( $line, $startPos,  $stopPos - $startPos - 1 );
		my $stop   = substr( $line, $stopPos,   $strandPos - $stopPos - 1 );
		my $strand = substr( $line, $strandPos, 1 );
		push( @loci,    $locus );
		push( @starts,  $start );
		push( @stops,   $stop );
		push( @strands, $strand );
	}
}
close $annotation;
print " with " . @loci . " loci (genes).\n";

#compute borders of promoters
my @positions;
my $skipNextLocus = 0;
for ( 0 .. $#loci ) {

	#two genes share the same promoter > special case 9
	#did computation with first gene > skip second gene
	if ($skipNextLocus) {
		$skipNextLocus = 0;
	}

	#first gene of the contig
	#AND NOT special case 9
	#[ $_ -1 ] not present
	elsif (
			$_ == 0
			&& !(    ( ( $strands[$_] eq '-' ) && ( $strands[ $_ + 1 ] eq '+' ) )
				  && ( $stops[$_] + $from >= $starts[ $_ + 1 ] - $from ) )
	  )
	{
		if ( $strands[$_] eq '+' ) {
			if ( ( ( $starts[$_] - $from ) >= 0 ) && ( $stops[$_] > $starts[$_] + $to ) ) {    #1
				push( @positions, $starts[$_] - $from );
				push( @positions, $starts[$_] + $to );
				push( @names,     "around_ts_of_" . $loci[$_] );
			}
			elsif ( ( ( $starts[$_] - $from ) < 0 ) && ( $stops[$_] > $starts[$_] + $to ) ) {    #2
				push( @positions, 1 );
				push( @positions, $starts[$_] + $to );
				push( @names,     "around_ts_of_" . $loci[$_] );
			}
			elsif ( ( ( $starts[$_] - $from ) >= 0 ) && ( $starts[$_] + $to >= $stops[$_] ) ) {    #3
				push( @positions, $starts[$_] - $from );
				push( @positions, $stops[$_] );
				push( @names,     "around_ts_of_" . $loci[$_] );
			}
			elsif ( ( ( $starts[$_] - $from ) < 0 ) && ( $starts[$_] + $to >= $stops[$_] ) ) {     #7
				push( @positions, 1 );
				push( @positions, $stops[$_] );
				push( @names,     "around_ts_of_" . $loci[$_] );
			}
			else {
				print "problem with promoter at gene " . $loci[$_] . "!\n";
			}
		}
		elsif ( $strands[$_] eq '-' ) {
			if ( ( $starts[$_] < $stops[$_] - $to ) && ( $starts[ $_ + 1 ] > $stops[$_] + $from ) ) {    #4
				push( @positions, $stops[$_] - $to );
				push( @positions, $stops[$_] + $from );
				push( @names,     "around_ts_of_" . $loci[$_] );
			}
			elsif ( ( $starts[$_] < $stops[$_] - $to ) && ( $starts[ $_ + 1 ] <= $stops[$_] + $from ) ) {   #5
				push( @positions, $stops[$_] - $to );
				push( @positions, $starts[ $_ + 1 ] - 1 );
				push( @names,     "around_ts_of_" . $loci[$_] );
			}
			elsif ( ( $starts[$_] >= $stops[$_] - $to ) && ( $starts[ $_ + 1 ] > $stops[$_] + $from ) ) {   #6
				push( @positions, $starts[$_] );
				push( @positions, $stops[$_] + $from );
				push( @names,     "around_ts_of_" . $loci[$_] );
			}
			elsif ( ( $starts[ $_ + 1 ] <= $stops[$_] + $from ) && ( $starts[$_] >= $stops[$_] - $to ) ) {  #8
				push( @positions, $starts[$_] );
				push( @positions, $starts[ $_ + 1 ] - 1 );
				push( @names,     "around_ts_of_" . $loci[$_] );
			}
			else {
				print "problem with promoter at gene " . $loci[$_] . "!\n";
			}
		}
		else {
			print "strand-problem at gene " . $loci[$_] . "!\n";
		}
	}

	#last gene of contig
	#[ $_ +1 ] not present
	elsif ( $_ == $#loci && !$skipNextLocus ) {
		if ( $strands[$_] eq '+' ) {
			if ( ( $stops[ $_ - 1 ] < $starts[$_] - $from ) && ( $stops[$_] > $starts[$_] + $to ) ) {    #1
				push( @positions, $starts[$_] - $from );
				push( @positions, $starts[$_] + $to );
				push( @names,     "around_ts_of_" . $loci[$_] );
			}
			elsif ( ( $stops[ $_ - 1 ] >= $starts[$_] - $from ) && ( $stops[$_] > $starts[$_] + $to ) ) {   #2
				push( @positions, $stops[ $_ - 1 ] + 1 );
				push( @positions, $starts[$_] + $to );
				push( @names,     "around_ts_of_" . $loci[$_] );
			}
			elsif ( ( $stops[ $_ - 1 ] < $starts[$_] - $from ) && ( $starts[$_] + $to >= $stops[$_] ) ) {   #3
				push( @positions, $starts[$_] - $from );
				push( @positions, $stops[$_] );
				push( @names,     "around_ts_of_" . $loci[$_] );
			}
			elsif ( ( $stops[ $_ - 1 ] >= $starts[$_] - $from ) && ( $starts[$_] + $to >= $stops[$_] ) ) {  #7
				push( @positions, $stops[ $_ - 1 ] + 1 );
				push( @positions, $stops[$_] );
				push( @names,     "around_ts_of_" . $loci[$_] );
			}
			else {
				print "problem with promoter at gene " . $loci[$_] . "!\n";
			}
		}
		elsif ( $strands[$_] eq '-' ) {
			if ( ( $starts[$_] < $stops[$_] - $to ) && ( ( $stops[$_] + $from ) <= $nucleotides ) ) {       #4
				push( @positions, $stops[$_] - $to );
				push( @positions, $stops[$_] + $from );
				push( @names,     "around_ts_of_" . $loci[$_] );
			}
			elsif ( ( $starts[$_] < $stops[$_] - $to ) && ( ( $stops[$_] + $from ) > $nucleotides ) ) {     #5
				push( @positions, $stops[$_] - $to );
				push( @positions, $nucleotides );
				push( @names,     "around_ts_of_" . $loci[$_] );
			}
			elsif ( ( $starts[$_] >= $stops[$_] - $to ) && ( ( $stops[$_] + $from ) <= $nucleotides ) ) {   #6
				push( @positions, $starts[$_] );
				push( @positions, $stops[$_] + $from );
				push( @names,     "around_ts_of_" . $loci[$_] );
			}
			elsif ( ( $starts[ $_ + 1 ] <= $stops[$_] + $from ) && ( ( $stops[$_] + $from ) > $nucleotides ) )
			{                                                                                               #8
				push( @positions, $starts[$_] );
				push( @positions, $nucleotides );
				push( @names,     "around_ts_of_" . $loci[$_] );
			}
			else {
				print "problem with promoter at gene " . $loci[$_] . "!\n";
			}
		}
		else {
			print "strand-problem at gene " . $loci[$_] . "!\n";
		}
	}

	#special-case 9
	elsif (    ( ( $strands[$_] eq '-' ) && ( $strands[ $_ + 1 ] eq '+' ) )
			&& ( $stops[$_] + $from >= $starts[ $_ + 1 ] - $from ) )
	{
		if ( ( $stops[$_] > $starts[$_] + $to ) && ( $starts[$_] < $stops[$_] - $to ) ) {    #9.1+4
			push( @positions, $stops[$_] - $to );
			push( @positions, $starts[ $_ + 1 ] + $to );
			push( @names,     "around_ts_of_" . $loci[$_] . "_and_" . $loci[ $_ + 1 ] );
		}
		elsif ( ( $starts[$_] < $stops[$_] - $to ) && ( $starts[ $_ + 1 ] + $to >= $stops[ $_ + 1 ] ) )
		{                                                                                    #9.3+4
			push( @positions, $stops[$_] - $to );
			push( @positions, $stops[ $_ + 1 ] );
			push( @names,     "around_ts_of_" . $loci[$_] . "_and_" . $loci[ $_ + 1 ] );
		}
		elsif ( ( $starts[$_] >= $stops[$_] - $to ) && ( $stops[ $_ + 1 ] > $starts[ $_ + 1 ] + $to ) )
		{                                                                                    #9.1+6
			push( @positions, $starts[$_] );
			push( @positions, $starts[ $_ + 1 ] + $to );
			push( @names,     "around_ts_of_" . $loci[$_] . "_and_" . $loci[ $_ + 1 ] );
		}
		elsif ( ( $starts[$_] >= $stops[$_] - $to ) && ( $starts[ $_ + 1 ] + $to >= $stops[ $_ + 1 ] ) )
		{                                                                                    #9.3+6
			push( @positions, $starts[$_] );
			push( @positions, $stops[ $_ + 1 ] );
			push( @names,     "around_ts_of_" . $loci[$_] . "_and_" . $loci[ $_ + 1 ] );
		}
		else {
			print "problem with promoter at gene " . $loci[$_] . "!\n";
		}
		$skipNextLocus = 1;
	}

	#normal-cases
	elsif ( !$skipNextLocus ) {
		if ( $strands[$_] eq '+' ) {
			if ( ( $stops[ $_ - 1 ] < $starts[$_] - $from ) && ( $stops[$_] > $starts[$_] + $to ) ) {    #1
				push( @positions, $starts[$_] - $from );
				push( @positions, $starts[$_] + $to );
				push( @names,     "around_ts_of_" . $loci[$_] );
			}
			elsif ( ( $stops[ $_ - 1 ] >= $starts[$_] - $from ) && ( $stops[$_] > $starts[$_] + $to ) ) {   #2
				push( @positions, $stops[ $_ - 1 ] + 1 );
				push( @positions, $starts[$_] + $to );
				push( @names,     "around_ts_of_" . $loci[$_] );
			}
			elsif ( ( $stops[ $_ - 1 ] < $starts[$_] - $from ) && ( $starts[$_] + $to >= $stops[$_] ) ) {   #3
				push( @positions, $starts[$_] - $from );
				push( @positions, $stops[$_] );
				push( @names,     "around_ts_of_" . $loci[$_] );
			}
			elsif ( ( $stops[ $_ - 1 ] >= $starts[$_] - $from ) && ( $starts[$_] + $to >= $stops[$_] ) ) {  #7
				push( @positions, $stops[ $_ - 1 ] + 1 );
				push( @positions, $stops[$_] );
				push( @names,     "around_ts_of_" . $loci[$_] );
			}
			else {
				print "problem with promoter at gene " . $loci[$_] . "!\n";
			}
		}
		elsif ( $strands[$_] eq '-' ) {
			if ( ( $starts[$_] < $stops[$_] - $to ) && ( $starts[ $_ + 1 ] > $stops[$_] + $from ) ) {       #4
				push( @positions, $stops[$_] - $to );
				push( @positions, $stops[$_] + $from );
				push( @names,     "around_ts_of_" . $loci[$_] );
			}
			elsif ( ( $starts[$_] < $stops[$_] - $to ) && ( $starts[ $_ + 1 ] <= $stops[$_] + $from ) ) {   #5
				push( @positions, $stops[$_] - $to );
				push( @positions, $starts[ $_ + 1 ] - 1 );
				push( @names,     "around_ts_of_" . $loci[$_] );
			}
			elsif ( ( $starts[$_] >= $stops[$_] - $to ) && ( $starts[ $_ + 1 ] > $stops[$_] + $from ) ) {   #6
				push( @positions, $starts[$_] );
				push( @positions, $stops[$_] + $from );
				push( @names,     "around_ts_of_" . $loci[$_] );
			}
			elsif ( ( $starts[ $_ + 1 ] <= $stops[$_] + $from ) && ( $starts[$_] >= $stops[$_] - $to ) ) {  #8
				push( @positions, $starts[$_] );
				push( @positions, $starts[ $_ + 1 ] - 1 );
				push( @names,     "around_ts_of_" . $loci[$_] );
			}
			else {
				print "problem with promoter at gene " . $loci[$_] . "!\n";
			}
		}
		else {
			print "strand-problem at gene " . $loci[$_] . "!\n";
		}
	}
}

#open output files
open( my $outputSequences, ">", "all_promoterSequences_of_" . $contig . ".fasta" )
  || die "outputSequences error: " . $!;
open( my $outputPositions, ">", "all_promoterPositions_of_" . $contig . ".csv" )
  || die "outputSequences error: " . $!;
print $outputPositions "name\tstart\tend\tlength\n";    #print head to positions file
my $skipNextPosition = 0;
my $currentRegion    = 0;
for ( 0 .. $#positions ) {

	#skip every other position
	#@positions[start, end, start, end, start, ...]
	#> go through all start positions only
	if ($skipNextPosition) {
		$skipNextPosition = 0;
	}
	else {
		$length = $positions[ $_ + 1 ] - $positions[$_] + 1;
		print $outputPositions $names[$currentRegion] . "\t"
		  . $positions[$_] . "\t"
		  . $positions[ $_ + 1 ] . "\t"
		  . $length . "\n";
		my $comment = ">" . $names[$currentRegion] . "___" . $length . "bp\n";
		my $sequence = substr( $sequence, $positions[$_] - 1, $length ) . "\n";
		$sequence =~ s/(.{60})/$1\n/g;    #maximum of 60 characters per line
		print $outputSequences $comment;
		print $outputSequences $sequence;
		++$currentRegion;
		$skipNextPosition = 1;
	}
}
close $outputPositions;
close $outputSequences;
print $currentRegion. " promoters extracted.";

#############################
#error detection / debugging#
#############################

$skipNextPosition = 0;
for ( 0 .. $#positions ) {
	if ( $_ < 0 ) {
		print "\nError: Computed a negative position!";
	}
	if ( $_ > $nucleotides ) {
		print "\nError: Computed a position \"after\" the last nucleotide!";
	}

	#skip every other position
	#@positions[start, end, start, end, start, ...]
	#> look at all start positions only
	if ($skipNextPosition) {
		$skipNextPosition = 0;
	}
	else {
		$skipNextPosition = 1;

		#no promoter length should be negative
		if ( ( $positions[ $_ + 1 ] - $positions[$_] + 1 ) < 0 ) {
			print "\nError: Computed a promoter with negative length: "
			  . ( $positions[ $_ + 1 ] - $positions[$_] + 1 ) . " "
			  . $_;
		}

		#no promoter length should be longer than twice the
		#assumed start and end of a promoter related to the transcription-start
		if ( ( $positions[ $_ + 1 ] - $positions[$_] + 1 ) > ( $from + $to ) * 2 + 1 ) {
			print "\nError: Computed a promoter that is too long: "
			  . ( $positions[ $_ + 1 ] - $positions[$_] + 1 ) . " "
			  . $_;
		}

		#the start of a promoter should always be before the stop and not vice versa or equal
		if ( $positions[$_] >= $positions[ $_ + 1 ] ) {
			print "\nError: Computed twisted or equal positions: "
			  . $positions[$_] . ">="
			  . $positions[ $_ + 1 ] . " "
			  . $_;
		}

		#no overlapping of start and end should happen
		if ( $_ != 0 && ( $positions[$_] <= $positions[ $_ - 1 ] ) ) {
			print "\nError: Computed start of current promoter before end of last promtor: "
			  . $positions[$_] . "<="
			  . $positions[ $_ - 1 ] . " "
			  . $_;
		}
	}
}
