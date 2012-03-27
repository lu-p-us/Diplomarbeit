#!/usr/bin/perl

#################################################################
#Extract the sequences of possible TFBS regions from the genome.#
#################################################################

#start and end of motifsearch region related to transcription-start
$from = 1000;
$to   = 50;

#genes of interest (aflatoxin-cluster in aspergillus flavus)
@genes = qw(
  AFL2G_07209
  AFL2G_07210
  AFL2G_07211
  AFL2G_07212
  AFL2G_07213
  AFL2G_07214
  AFL2G_07215
  AFL2G_07216
  AFL2G_07217
  AFL2G_07218
  AFL2G_07219
  AFL2G_07220
  AFL2G_07221
  AFL2G_07222
  AFL2G_07223
  AFL2G_07224
  AFL2G_07225
  AFL2G_07226
  AFL2G_07227
  AFL2G_07228
  AFL2G_07229
  AFL2G_07230);
@genes            = sort(@genes);
$contig           = "supercont2.7";
$tabsBeforeLOCUS  = 0;
$tabsBeforeSTART  = 4;
$tabsBeforeSTOP   = 5;
$tabsBeforeSTRAND = 6;

#open annotation file
open( $annotation, "<", "../aspergillus_flavus_2_genome_summary_per_gene.txt" )
  || die "annotation input error";

#analyse each line of the annotation file
<$annotation>;    #skip first line
while (<$annotation>) {
	$line = $_;

	#extract LOCUS, START, STOP and STRAND
	$noOfTab = 1;
	$tabPos  = 0;
	do {
		$tabPos = index( $line, "\t", $tabPos + 1 );
		if ( $noOfTab == $tabsBeforeLOCUS + 1 ) {
			$LOCUS = substr( $line, 0, $tabPos );
		}
		elsif ( $noOfTab == $tabsBeforeSTART ) {
			$START_pos = $tabPos + 1;
		}
		elsif ( $noOfTab == $tabsBeforeSTOP ) {
			$STOP_pos = $tabPos + 1;
		}
		elsif ( $noOfTab == $tabsBeforeSTRAND ) {
			$STRAND_pos = $tabPos + 1;
		}
		++$noOfTab;
	} while ( $tabPos != -1 );
	$START  = substr( $line, $START_pos,  $STOP_pos - $START_pos - 1 );
	$STOP   = substr( $line, $STOP_pos,   $STRAND_pos - $STOP_pos - 1 );
	$STRAND = substr( $line, $STRAND_pos, 1 );
	push( @LOCUSES, $LOCUS );
	push( @STARTS,  $START );
	push( @STOPS,   $STOP );
	push( @STRANDS, $STRAND );
}
close $annotation;

#compute borders of regions for motifsearch
$currentGene = 0;
$skipNext    = 0;
for ( 0 .. $#LOCUSES ) {

	#select genes of interest
	if ( $currentGene <= $#genes && ( $genes[$currentGene] eq $LOCUSES[$_] ) ) {
		if ($skipNext) {
			$skipNext = 0;
		}
		elsif ( $_ == 0 ) {    #[ $_ -1 ] not present
			print "cluster at begin of genome/contig!\n";
		}
		elsif ( $_ == $#LOCUSES ) {    #[ $_ +1 ] not present
			print "cluster at end of genome/contig!\n";
		}

		#special-case 9
		elsif (    ( ( $STRANDS[$_] eq '-' ) && ( $STRANDS[ $_ + 1 ] eq '+' ) )
				&& ( $STOPS[$_] + $from >= $STARTS[ $_ + 1 ] - $from ) )
		{
			if ( ( $STOPS[$_] > $STARTS[$_] + $to ) && ( $STARTS[$_] < $STOPS[$_] - $to ) ) {    #9.1+4
				push( @positions, $STOPS[$_] - $to );
				push( @positions, $STARTS[ $_ + 1 ] + $to );
				push( @names,     "around_ts_of_" . $LOCUSES[$_] . "_and_" . $LOCUSES[ $_ + 1 ] );
			}
			elsif ( ( $STARTS[$_] < $STOPS[$_] - $to ) && ( $STARTS[ $_ + 1 ] + $to >= $STOPS[ $_ + 1 ] ) )
			{                                                                                    #9.3+4
				push( @positions, $STOPS[$_] - $to );
				push( @positions, $STOPS[ $_ + 1 ] );
				push( @names,     "around_ts_of_" . $LOCUSES[$_] . "_and_" . $LOCUSES[ $_ + 1 ] );
			}
			elsif ( ( $STARTS[$_] >= $STOPS[$_] - $to ) && ( $STOPS[ $_ + 1 ] > $STARTS[ $_ + 1 ] + $to ) )
			{                                                                                    #9.1+6
				push( @positions, $STARTS[$_] );
				push( @positions, $STARTS[ $_ + 1 ] + $to );
				push( @names,     "around_ts_of_" . $LOCUSES[$_] . "_and_" . $LOCUSES[ $_ + 1 ] );
			}
			elsif ( ( $STARTS[$_] >= $STOPS[$_] - $to ) && ( $STARTS[ $_ + 1 ] + $to >= $STOPS[ $_ + 1 ] ) )
			{                                                                                    #9.3+6
				push( @positions, $STARTS[$_] );
				push( @positions, $STOPS[ $_ + 1 ] );
				push( @names,     "around_ts_of_" . $LOCUSES[$_] . "_and_" . $LOCUSES[ $_ + 1 ] );
			}
			else {
				print "problem with motif-search-region at gene " . $LOCUSES[$_] . "!\n";
			}
			$skipNext = 1;
		}

		#normal-cases
		elsif ( !$skipNext ) {
			if ( $STRANDS[$_] eq '+' ) {
				if ( ( $STOPS[ $_ - 1 ] < $STARTS[$_] - $from ) && ( $STOPS[$_] > $STARTS[$_] + $to ) ) {   #1
					push( @positions, $STARTS[$_] - $from );
					push( @positions, $STARTS[$_] + $to );
					push( @names,     "around_ts_of_" . $LOCUSES[$_] );
				}
				elsif ( ( $STOPS[ $_ - 1 ] >= $STARTS[$_] - $from ) && ( $STOPS[$_] > $STARTS[$_] + $to ) )
				{                                                                                           #2
					push( @positions, $STOPS[ $_ - 1 ] + 1 );
					push( @positions, $STARTS[$_] + $to );
					push( @names,     "around_ts_of_" . $LOCUSES[$_] );
				}
				elsif ( ( $STOPS[ $_ - 1 ] < $STARTS[$_] - $from ) && ( $STARTS[$_] + $to >= $STOPS[$_] ) )
				{                                                                                           #3
					push( @positions, $STARTS[$_] - $from );
					push( @positions, $STOPS[$_] );
					push( @names,     "around_ts_of_" . $LOCUSES[$_] );
				}
				elsif ( ( $STOPS[ $_ - 1 ] >= $STARTS[$_] - $from ) && ( $STARTS[$_] + $to >= $STOPS[$_] ) )
				{                                                                                           #7
					push( @positions, $STOPS[ $_ - 1 ] + 1 );
					push( @positions, $STOPS[$_] );
					push( @names,     "around_ts_of_" . $LOCUSES[$_] );
				}
				else {
					print "problem with motif-search-region at gene " . $LOCUSES[$_] . "!\n";
				}
			}
			elsif ( $STRANDS[$_] eq '-' ) {
				if ( ( $STARTS[$_] < $STOPS[$_] - $to ) && ( $STARTS[ $_ + 1 ] > $STOPS[$_] + $from ) ) {   #4
					push( @positions, $STOPS[$_] - $to );
					push( @positions, $STOPS[$_] + $from );
					push( @names,     "around_ts_of_" . $LOCUSES[$_] );
				}
				elsif ( ( $STARTS[$_] < $STOPS[$_] - $to ) && ( $STARTS[ $_ + 1 ] <= $STOPS[$_] + $from ) )
				{                                                                                           #5
					push( @positions, $STOPS[$_] - $to );
					push( @positions, $STARTS[ $_ + 1 ] - 1 );
					push( @names,     "around_ts_of_" . $LOCUSES[$_] );
				}
				elsif ( ( $STARTS[$_] >= $STOPS[$_] - $to ) && ( $STARTS[ $_ + 1 ] > $STOPS[$_] + $from ) )
				{                                                                                           #6
					push( @positions, $STARTS[$_] );
					push( @positions, $STOPS[$_] + $from );
					push( @names,     "around_ts_of_" . $LOCUSES[$_] );
				}
				elsif ( ( $STARTS[ $_ + 1 ] <= $STOPS[$_] + $from ) && ( $STARTS[$_] >= $STOPS[$_] - $to ) )
				{                                                                                           #8
					push( @positions, $STARTS[$_] );
					push( @positions, $STARTS[ $_ + 1 ] - 1 );
					push( @names,     "around_ts_of_" . $LOCUSES[$_] );
				}
				else {
					print "problem with motif-search-region at gene " . $LOCUSES[$_] . "!\n";
				}
			}
			else {
				print "strand-problem at gene " . $LOCUSES[$_] . "!\n";
			}
		}
		++$currentGene;
	}
}

#open genome file
open( $genome, "<", "../aspergillus_flavus_2_supercontigs.fasta" ) || die "genome input error";

#read the genome line by line
$found = 0;
while ( !eof($genome) && !$found ) {
	$line = <$genome>;
	if ( index( $line, ">" ) != -1 && index( $line, $contig ) != -1 ) {
		$found = 1;
		$line  = <$genome>;
		while ( !eof && index( $line, ">" ) == -1 ) {
			$sequence .= $line;
			chomp($sequence);
			$line = <$genome>;
		}
	}
}
close $genome;

#open output file
open( $output, ">", "sequences_bs_regions_" . $from . "_" . $to . ".fasta" ) || die "output error";
$skipNext      = 0;
$currentRegion = 0;
for ( 0 .. $#positions ) {
	if ($skipNext) {
		$skipNext = 0;
	}
	else {
		$length         = $positions[ $_ + 1 ] - $positions[$_] + 1;
		$outputComment  = "";
		$outputSequence = "";
		$outputComment  = ">" . $names[ $currentRegion++ ] . "___" . $length . "bp\n";
		$outputSequence = substr( $sequence, $positions[$_] - 1, $length ) . "\n";
		$outputSequence =~ s/(.{60})/$1\n/g;    #maximum of 60 characters per line
		print $output $outputComment;
		print $output $outputSequence;
		$skipNext = 1;
	}
}
close $output;
print $currentRegion; #print number of extracted promoter regions
