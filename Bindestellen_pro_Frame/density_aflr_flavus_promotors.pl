#!/usr/bin/perl

##############################################
#Count the number of binding sites per frame.#
#Only use binding site from promoter regions.#
##############################################

use Chart::Clicker;
use Chart::Clicker::Renderer::Area;
use Graphics::Primitive::Brush;
use List::Util ("max");

#thousands seperator
sub thousands {
	local $_ = shift;
	1 while s/^([-+]?\d+)(\d{3})/$1.$2/;
	return $_;
}

###############
#ALL promoters#
###############

#rebuild @positions and @loci from all_promotors_of.pl
my $contig = "supercont2.7";
my @positions;
my @lociAll;

#open file with promotor positions
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
open( my $converted, ">", "fimo_converted.txt" ) || die "converted output error: " . $!;
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

######################
#count hits per frame#
######################

#sort numerically ascending
@starts = sort { $a <=> $b } (@starts);

#define frames
my $frameLength     = 50000;
#my $frameLength     = 10000;
my $noOfNucleotides = 2337902;
my $frameShift      = 25000;
#my $frameShift      = 5000;
my $maxX            = 0;
my $maxY            = 0;
my @allSeries;
my $frameStart = 0;
my $frameEnd   = $frameStart + $frameLength;
my @counts;

#count hits per frame
while ( $frameEnd <= $noOfNucleotides ) {
	push( @counts, 0 );
	for ( 0 .. $#starts ) {
		last
		  if ( $starts[$_] > $frameEnd );  #all upcoming starts downstream the current frame > frame completed
		if ( ( $starts[$_] > $frameStart ) && ( $starts[$_] <= $frameEnd ) ) {
			++$counts[-1];
		}
	}
	$frameStart += $frameShift;
	$frameEnd   += $frameShift;
}
print "Counted hits for "
  . thousands( scalar(@counts) )
  . " frames with a length of "
  . $frameLength
  . " nucleotides per frame.\n";

#############
#text output#
#############

#write counts to file
open( my $output, ">", "counts" . $frameLength . "_promoters.csv" ) || die "output error: $!";
print $output "hits per frame" . "\t" . "first frame-nucleotide" . "\t" . "last frame-nucleotide" . "\n";
my $lines;
$frameStart = 0;
$frameEnd   = $frameStart + $frameLength;
for ( 0 .. $#counts ) {

	#	if($counts[$_] > 0){
	++$lines;
	print $output $counts[$_] . "\t" . ( $frameStart + 1 ) . "\t" . ( $frameEnd + 1 ) . "\n";

	#	}
	$frameStart += $frameShift;
	$frameEnd   += $frameShift;
}
print thousands($lines) . " lines written to output file.\n";
close $output;

##################
#graphical output#
##################

#build the chart
my $chart = Chart::Clicker->new( format => "pdf", width => 4000, height => 500 );
$chart->legend->visible(0);
$chart->border->width(0);

#set renderer (line is default)
my $area = Chart::Clicker::Renderer::Area->new( opacity => 0.5,
												brush   => Graphics::Primitive::Brush->new( width => 1 ) );
$chart->get_context("default")->renderer($area);

#compute x-axis ticks
my @xTicks        = (1);               #starting x tick
my $ending_x_tick = scalar(@counts);
my $xStepwidth    = 20;
my $i             = $xStepwidth;
do {
	push( @xTicks, $i );
	$i += $xStepwidth;
} while ( $i < $ending_x_tick );
push( @xTicks, $ending_x_tick );

#compute y-axis ticks
my @yTicks        = (0);               #starting y tick
my $ending_y_tick = max(@counts);
my $yStepwidth    = 8;
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
$yContext->label("Bindestellen");
$yContext->format("%d");
$yContext->tick_values( [@yTicks] );
$yContext->label_font->size(45);
$yContext->tick_font->size(40);

#build the series
#array with frames > x-axis
#array with counts > y-axis
my @frames = ( 1 .. @counts );
my $series = Chart::Clicker::Data::Series->new( keys => \@frames, values => \@counts );
$series->name(   "Binding-sites per frame. Frame length: "
			   . $frameLength
			   . " nucleotides. Number of frames: "
			   . thousands( scalar(@counts) )
			   . "." );
push( @allSeries, $series );

#build the dataset
my $dataset = Chart::Clicker::Data::DataSet->new( series => \@allSeries );

#add the dataset to the chart
$chart->add_to_datasets($dataset);

#write the chart to a file
$chart->write_output( "counts" . $frameLength . ".pdf" );
