#!/usr/bin/perl

##############################################
#Count the number of binding sites per frame.#
##############################################

use Chart::Clicker;
use Chart::Clicker::Renderer::Area;
use Graphics::Primitive::Brush;
use List::Util ("max");

open( my $input, "<", "fimo_contig.txt" ) || die "input error: $!";
my @starts;

#skip first line
<$input>;

#extract start-numbers
while (<$input>) {
	my $start = $_;
	$start = substr( $start, index( $start, "\t" ) + 1 );
	$start = substr( $start, index( $start, "\t" ) + 1 );
	$start = substr( $start, 0, index( $start, "\t" ) );
	push( @starts, $start );
}
close $input;

#sort numerically ascending
@starts = sort { $a <=> $b } (@starts);

#define frames
my $frameLength     = 10000;
#my $frameLength     = 50000;
my $noOfNucleotides = 2337902;
my $noOfFrames      = int( $noOfNucleotides / ( $frameLength * 0.5 ) );    #overlapping
print $noOfFrames . " full-length frames\n";
my $frameStart   = 0;
my $frameEnd     = $frameStart + $frameLength;
my $currentFrame = 0;
my @counts;
my @frameStarts;
my @frameEnds;

#init counts
for ( 0 .. $noOfFrames ) {
	push( @counts, 0 );
}

#count hits per frame
while ( $currentFrame < $noOfFrames ) {
	for ( 0 .. $#starts ) {
		if ( ( $starts[$_] > $frameStart ) && ( $starts[$_] <= $frameEnd ) ) {
			++$counts[$currentFrame];
		}
	}
	++$currentFrame;
	push( @frameStarts, $frameStart );
	push( @frameEnds,   $frameEnd );
	$frameStart = $frameEnd - 0.5 * $frameLength;    #overlapping
	$frameEnd   = $frameStart + $frameLength;
}

#last frame
#the remaining nucleotides after to last full-length frame to the end of the contig
$frameEnd = $noOfNucleotides;
for ( 0 .. $#starts ) {
	if ( ( $starts[$_] > $frameStart ) && ( $starts[$_] <= $frameEnd ) ) {
		++$counts[$currentFrame];
	}
}
push( @frameStarts, $frameStart );
push( @frameEnds,   $frameEnd );

#############
#text output#
#############

#write counts to file
open( my $output, ">", "counts" . $frameLength . "_contig.csv" ) || die "output error: $!";
$currentFrame = 0;
print $output "framenumber" . "\t"
  . "hits per frame" . "\t"
  . "first frame-nucleotide" . "\t"
  . "last frame-nucleotide" . "\n";
for ( 0 .. $#counts ) {
	print $output ++$currentFrame . "\t"
	  . $counts[$_] . "\t"
	  . ( $frameStarts[$_] + 1 ) . "\t"
	  . $frameEnds[$_] . "\n";
}
print $output "\n";
close $output;

##################
#graphical output#
##################

#build the chart
my $chart = Chart::Clicker->new( format => "pdf", width => 4000, height => 500 );
$chart->border->width(0);
$chart->legend->visible(0);

#set renderer (line is default)
my $area = Chart::Clicker::Renderer::Area->new( opacity => 0.5,
												brush   => Graphics::Primitive::Brush->new( width => 1 ) );
$chart->get_context("default")->renderer($area);

#compute x-axis ticks
my @xTicks        = (1);               #starting x tick
my $ending_x_tick = scalar(@counts);
my $xStepwidth    = 100;
my $i             = $xStepwidth;
do {
	push( @xTicks, $i );
	$i += $xStepwidth;
} while ( $i < $ending_x_tick );
push( @xTicks, $ending_x_tick );

#compute y-axis ticks
my @yTicks        = (0);               #starting y tick
my $ending_y_tick = max(@counts);
my $yStepwidth    = 3;
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
my @frames = ( 1 .. $noOfFrames + 1 );
my $series = Chart::Clicker::Data::Series->new( keys => \@frames, values => \@counts );
$series->name("Binding-sites per frame");

#build the dataset
my $dataset = Chart::Clicker::Data::DataSet->new( series => [$series] );

#add the dataset to the chart
$chart->add_to_datasets($dataset);

#write the chart to a file
$chart->write_output( "counts" . $frameLength . "_contig.pdf" );
