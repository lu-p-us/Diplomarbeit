#!/usr/bin/perl

##########################
#Compute score per frame.#
##########################

use Chart::Clicker;
use Chart::Clicker::Renderer::Area;
use Graphics::Primitive::Brush;
use List::Util qw(max sum);

#delete multiple occurrences in array
sub del_double {
	my %all = ();
	@all{@_} = 1;
	return ( keys %all );
}

###############
#ALL promoters#
###############

#rebuild @positions and @loci from all_promotors_of.pl
my $contig = "supercont2.7";
my @positions;
my @lociAll;

#open file with promoter positions
open( my $positionsInput, "<", "../all_promotorPositions_of_" . $contig . ".csv" )
  || die "positions input error: " . $!;

#extract positions and loci
<$positionsInput>;
while (<$positionsInput>) {
	my $line          = $_;
	my $firstTab      = index( $line, "\t" );
	my $locus         = substr( $line, 0, index( $line, "\t" ) );
	my $promotorStart = substr( $line, $firstTab + 1 );
	my $secondTab     = index( $promotorStart, "\t" );
	my $promotorEnd   = substr( $promotorStart, $secondTab + 1 );
	my $thirdTab      = index( $promotorEnd, "\t" );
	$promotorStart = substr( $promotorStart, 0, $secondTab );
	$promotorEnd   = substr( $promotorEnd,   0, $thirdTab );
	push( @lociAll,   $locus );
	push( @positions, $promotorStart );
	push( @positions, $promotorEnd );
}
close $positionsInput;
print "Contig " . $contig . " has " . @lociAll . " promoters (binding site regions).\n";

##############################
#promoters with BINDING SITES#
##############################

#open file with binding site positions computed by FIMO
my @starts;
my $bindingSites;
open( my $fimo, "<", "fimo_promoters.txt" ) || die "FIMO input error: " . $!;
<$fimo>;

#extract start-numbers
while (<$fimo>) {
	++$bindingSites;
	my $start = $_;
	$start = substr( $start, index( $start, "\t" ) + 1 );
	my $sequenceName = $start;
	$sequenceName = substr( $sequenceName, 0, index( $sequenceName, "\t" ) );
	$sequenceName = substr( $sequenceName, 0, index( $sequenceName, "___" ) );
	$start = substr( $start, index( $start, "\t" ) + 1 );
	my $stop = $start;
	$start = substr( $start, 0, index( $start, "\t" ) );
	$stop = substr( $stop, index( $stop, "\t" ) + 1 );
	$stop = substr( $stop, 0, index( $stop, "\t" ) );
	push( @loci,   $sequenceName );
	push( @starts, $start );
	push( @stops,  $stop );
}
close $fimo;
print "FIMO (promoters only) found " . $bindingSites . " binding sites within these regions.\n";

#convert promoter-positions of binding sites to contig-positions of binding sites
open( my $converted, ">", "fimo_converted.csv" ) || die "converted output error: " . $!;
print $converted "sequence name\tstart\tstop\n";
for ( 0 .. $#loci ) {
	my $p            = $_;
	my $currentLocus = $loci[$p];
	for ( 0 .. $#lociAll ) {
		my $a = $_;
		if ( $currentLocus eq $lociAll[$a] ) {
			$starts[$p] += ( $positions[ $a * 2 ] - 1 );
			$stops[$p]  += ( $positions[ $a * 2 ] - 1 );
			print $converted $currentLocus . "\t" . $starts[$p] . "\t" . $stops[$p] . "\n";
			last;    #break out of lociAll loop
		}
	}
}

#########################################
#count binding sites per promoter region#
#########################################

@starts = sort { $a <=> $b } (@starts);    #sort numerically ascending
@stops  = sort { $a <=> $b } (@stops);     #sort numerically ascending
my $skipNextPosition = 0;
my @counts;
for ( 0 .. $#positions ) {
	$pos = $_;

	#skip every other position
	#@positions[start, end, start, end, start, ...]
	#> look at all start positions only
	if ($skipNextPosition) {
		$skipNextPosition = 0;
	}
	else {
		$skipNextPosition = 1;
		push( @counts, 0 );
		for ( 0 .. $#starts ) {
			if ( $starts[$_] >= $positions[$pos] && $stops[$_] <= $positions[ $pos + 1 ] ) {
				++$counts[-1];
			}
		}
	}
}

#write binding sites per promoter region to file
open( $hitsPerPro, ">", "bs_per_promoter.csv" ) || die "bs per promoter output error: " . $!;
print $hitsPerPro "promoter number\tlocus name\tfirst nucleotide\tlast nucleotide\tbinding sites\n";
for ( 0 .. $#counts ) {
	print $hitsPerPro ( $_ + 1 ) . "\t"
	  . $lociAll[$_] . "\t"
	  . $positions[ $_ * 2 ] . "\t"
	  . $positions[ $_ * 2 + 1 ] . "\t"
	  . $counts[$_] . "\n";
}
close $hitsPerPro;

################
#frames/scoring#
################

my $penalty = 3;
my $reward  = 1;

#scoring function
sub score {
	$x = shift;
	if ( $x > 0 ) {

		#promoter with one or more binding sites
		#increase score by reward * number of binding sites
		return $x * $reward;
	}
	else {

		#promoter without a binding site
		#decrease score by penalty value
		return -$penalty;
	}
}

my @frameLengths = ( 1 ... 25 );
#my @frameLengths = (5, 10, 15, 20);
open( my $heatmap, ">", "heatmap.csv" ) || die "heatmap output error: " . $!;
my $frameNumber;
my $na            = 0;
my $maxX          = 0;
my $maxY          = 0;
my $maxAppearance = 0;
foreach $frameLength (@frameLengths) {
	$frameNumber = 1;
	my $frameStart = 0;
	my $frameEnd   = $frameStart + $frameLength;
	my @frameScores;
	while ( $frameEnd <= @counts ) {
		my $score = 0;
		for ( 0 .. $frameLength - 1 ) {
			my $promoter = $frameStart + $_;
			$score += score( $counts[$promoter] );    #r=1, p=3

			#			if($counts[$promoter] > 0){	  #r=1, no penalty
			#				$score += 1;
			#			}
		}
		if ( $score <= 0 ) {
			$score = 0;                               #do not allow negative scores
		}
		push( @frameScores, $score );
		++$frameStart;
		++$frameEnd;
	
		###########
		#3D output#
		###########
		
		#heatmap-matrix for R
		print $heatmap $score . "\t";
	}
	my $i = $na;
	while ( $i-- > 0 ) {
		print $heatmap "NA\t";
	}
	++$na;
	print $heatmap "\n";
	write score per frame to file
	  open( $scorePerFra, ">", "scores_l=" . $frameLength . "_r=" . $reward . "_p=" . $penalty . ".csv" )
	  || die "scores output error: " . $!;
	#	open( $scorePerFra, ">", "scores_l=" . $frameLength . ".csv" )	#without penalty value
	#	  || die "scores output error: " . $!;
	print $scorePerFra "frame number\tfirst promoter\tlast promoter\tscore\n";
	for ( 0 .. $#frameScores ) {
		my $firstPro = $_;
		my $lastPro  = $_ + $frameLength;
		print $scorePerFra ( $_ + 1 ) . "\t"
		  . ( $firstPro + 1 ) . "\t"
		  . $lastPro . "\t"
		  . $frameScores[$_] . "\n";
	}
	close $scorePerFra;

	##################
	#graphical output#
	##################

	#build the series
	#array with frames > x-axis
	#array with counts > y-axis
	my @frames = ( 1 .. @frameScores );
	my $series = Chart::Clicker::Data::Series->new( keys => \@frames, values => \@frameScores );
	$series->name( $frameLength . " Promotoren pro Frame" );
	push( @allSeries, $series );

	#compute maxium x-axis and y-axis values
	if ( scalar(@frameScores) > $maxX ) {
		$maxX = scalar(@frameScores);
	}
	if ( max(@frameScores) > $maxY ) {
		$maxY = max(@frameScores);
	}
	
	################
	#console output#
	################
	
	#highest score of each length
	my $lengthsMaxScore = max(@frameScores);
	if ( $lengthsMaxScore > 0 )    #lengths with lower scores will not be printed to console
	{
		print "\nThe frame(s) with the best score of "
		  . $lengthsMaxScore
		  . " in all frames of length "
		  . $frameLength
		  . " is(are) frame(s) no. ";
		my $string;
		my $highscores = 0;
		for ( 0 .. $#frameScores ) {
			if ( $frameScores[$_] == $lengthsMaxScore ) {
				$string .= ( $_ + 1 ) . ", ";
				++$highscores;
			}
		}
		substr( $string, length($string) - 2, 2, "." );
		print $string. "\n";
		print $highscores. " highscore(s).\n";
	}

	#arithmetic mean of scores
	my $average = sum(@frameScores) / @frameScores;

	#normal distribution?
	my @differentScores = sort { $a <=> $b } ( del_double(@frameScores) );
	my @distribution;
	for ( 0 .. $#differentScores ) {
		push( @distribution, 0 );
	}
	for ( 0 .. $#differentScores ) {
		my $diffScorePos = $_;
		my $diffScore    = $differentScores[$_];
		for ( 0 .. $#frameScores ) {
			if ( $frameScores[$_] == $diffScore ) {
				++$distribution[$diffScorePos];
			}
		}
	}
	if ( $maxAppearance < max(@distribution) ) {
		$maxAppearance = max(@distribution);
	}
	print "different scores: @differentScores\n";
	print "distribution of scores: @distribution\n";

	#standardization of highscores
	my $variance = 0;
	for ( 0 .. $#frameScores ) {
		$variance += ( $frameScores[$_] - $average )**2;
	}
	$variance /= ( @frameScores - 1 );
	my $standardDeviation = sqrt($variance);
	my $standardScore     = ( $lengthsMaxScore - $average ) / $standardDeviation;
	print "average = $average\n";
	print "variance = $variance\n";
	print "standard deviation = $standardDeviation\n";
	print "standard score = $standardScore\n";
	print "highscore - (3 x standard deviation + average) = "
	  . ( $lengthsMaxScore - ( 3 * $standardDeviation + $average ) ) . "\n";

#	#highest score in all lengths
#	if ( $lengthsMaxScore > $maxX )
#	{
#		$maxX = $lengthsMaxScore;
#	}
#	if ( $frameLength == $frameLengths[-1] )
#	{    #last frame length in process
#		print "\n";
#		print "The highest score is " . $maxY . ".";
#	}

#	#highest standard score in all lengths
#	if ( $standardScore > $maxY )
#	{
#		$maxY = $standardScore;
#	}
#	if ( $frameLength == $frameLengths[-1] )
#	{    #last frame length in process
#		print "\n";
#		print "The highest standard score is " . $maxY . ".";
#	}

#	#highest difference between highscore and (3 x standard deviation + average)
#	if ( ($lengthsMaxScore - (3*$standardDeviation+$average)) > $maxY )
#	{
#		$maxY = ($lengthsMaxScore - (3*$standardDeviation+$average));
#	}
#	if ( $frameLength == $frameLengths[-1] )
#	{    #last frame length in process
#		print "\n";
#		print "The highest difference between \"highscore and (3 x standard deviation + average)\" is " . $maxY . ".";
#	}
}
close $heatmap;

#build the chart
my $chart = Chart::Clicker->new( format => "pdf", width => 4000, height => 500 );
$chart->border->width(0);
$chart->legend->font->size(40);

#set renderer (line is default)
my $area = Chart::Clicker::Renderer::Area->new( opacity => 0.5,
												brush   => Graphics::Primitive::Brush->new( width => 1 ) );
$chart->get_context("default")->renderer($area);

#compute x-axis ticks
my @xTicks        = (1);           #starting x tick
my $ending_x_tick = $maxX;
my $xStepwidth    = 90;
my $i             = $xStepwidth;
do {
	push( @xTicks, $i );
	$i += $xStepwidth;
} while ( $i < $ending_x_tick );
push( @xTicks, $ending_x_tick );

#compute y-axis ticks
my @yTicks        = (0);           #starting y tick
my $ending_y_tick = $maxY;
my $yStepwidth    = 4;
my $j             = $yStepwidth;
do {
	push( @yTicks, $j );
	$j += $yStepwidth;
} while ( $j < $ending_y_tick );
push( @yTicks, $ending_y_tick );

#format x-axis
my $xContext = $chart->get_context("default")->domain_axis;
$xContext->label("Frame-Nummer");
$xContext->format("%d");
$xContext->tick_values( [@xTicks] );
$xContext->label_font->size(45);
$xContext->tick_font->size(40);

#format y-axis
my $yContext = $chart->get_context("default")->range_axis;
$yContext->label("Frame-Score");
$yContext->format("%d");
$yContext->tick_values( [@yTicks] );
$yContext->label_font->size(45);
$yContext->tick_font->size(40);

#build the dataset
my $dataset = Chart::Clicker::Data::DataSet->new( series => \@allSeries );

#add the dataset to the chart
$chart->add_to_datasets($dataset);

#write the chart to a file
$chart->get_context("default")->renderer($area);
$chart->write_output( "scores_r=" . $reward . "_p=" . $penalty . "_area.pdf" );
#$chart->get_context("default")->renderer($area);
#$chart->write_output( "scores_area.pdf" );

#######
#zoom#
#######

#set range
my $min = 565;
my $max = 585;
$xContext->range( Chart::Clicker::Data::Range->new( min => $min, max => $max ) );
@zoomTicks = ( $min .. $max );
$xContext->tick_values( \@zoomTicks );

#$chart->width( ( $max - $min ) * 50 );
#write to file
#$chart->get_context("default")->renderer($area);
#$chart->write_output( "scores_r=" . $reward . "_p=" . $penalty . "_area_zoom.pdf" );

