#!/usr/bin/perl

#############################
#TFBS distribution in FunTF.#
#############################

use Chart::Clicker;
use Chart::Clicker::Data::Marker;
use List::Util 'max';

#min and max positions of binding sites
$min = 0;
$max = 0;

#open file
open( $fh, "<", "Location_Site.csv" ) || die "cannot open file";

#read file line by line
while (<$fh>) {
	$line = $_;

	#find first ";"
	$index1 = index( $line, ";" );

	#find second ";"
	$index2 = index( $line, ";", $index1 + 1 );

	#find third ";"
	$index3 = index( $line, ";", $index2 + 1 );

	#isolate position "from"
	$substring_from = substr( $line, $index1 + 2, $index2 - $index1 - 3 );

	#isolate position "to"
	$substring_to = substr( $line, $index2 + 2, $index3 - $index2 - 3 );

	#substitute "around" with 4 positions plus and minus the "around position"
	$index_around1 = index( $substring_from, "arnd" );
	if ( $index_around1 != -1 ) {
		$substring_from = substr( $substring_from, 5 );
		$substring_from -= 4;
		$substring_to = substr( $substring_to, 5 );
		$substring_to += 4;
	}

	#compute min and max positions
	if ( $substring_from < $min ) {
		$min = $substring_from;
	}
	if ( $substring_to > $max ) {
		$max = $substring_to;
	}

	#save positions and length in arrays
	if ( $substring_from <= $substring_to ) {
		if ( $substring_from <= 0 && $substring_to <= 0 ) {
			$l = -( $substring_from - $substring_to ) + 1;
		}
		elsif ( $substring_from >= 0 && $substring_to >= 0 ) {
			$l = ( $substring_to - $substring_from ) + 1;
		}
		elsif ( $substring_from <= 0 && $substring_to >= 0 ) {
			$l = -($substring_from) + $substring_to + 1;
		}
		push( @length, $l );
		push( @from,   $substring_from );
		push( @to,     $substring_to );
		++$sites;
	}
}

#close file
close $fh;

#initiate counts
for ( $min .. $max ) {
	push( @counts, 0 );
}

#fill counts
#increase counts for positions of @from to @to and shift per $min:
#range $min..$max to range 0..overall-length
for ( 0 .. $#from ) {
	$x = $from[$_];    #startposition of binding-site
	while ( $x <= $to[$_] ) {
		$counts[ $x - $min ]++;
		++$x;
	}
}

#compute min and max length
$min_length = -$min + $max;
$max_length = 0;
for ( 0 .. $#length ) {
	if ( $length[$_] < $min_length ) {
		$min_length = $length[$_];
	}
	if ( $length[$_] > $max_length ) {
		$max_length = $length[$_];
	}
}

#compute median length and average length
@slength    = sort(@length);
$med_length = @slength[ int( $#length / 2 ) ];    #int() > cutoff the decimal part > under-median
$sum_length = 0;
for ( 0 .. $#length ) {
	$sum_length += $length[$_];
}
$avg_length = $sum_length / $#length;

#compute positions with most binding-sites
$max_bs = max(@counts);
for ( 0 .. $#counts ) {
	if ( $counts[$_] == $max_bs ) {
		push( @max_bs_pos, ( $_ + $min ) );
	}
}

#open(create) output-file
open( $fh, ">", "distribution_fun.txt" ) || die "cannot create file";

#print length-statistics in file
print $fh "minimum length = " . $min_length . "\n";
print $fh "maximum length = " . $max_length . "\n";
print $fh "average length = " . $avg_length . "\n";
print $fh "median length = " . $med_length . "\n";
print $fh "most binding sites at one position = " . $max_bs . "\n";
print $fh "most binding sites at position(s) ";
print $fh "@max_bs_pos \n";
print $fh "number of binding sites = " . $sites . "\n" . "\n";

#print counts in file
for ( 0 .. $#counts ) {
	print $fh $_ + $min . ";" . $counts[$_] . "\n";
}

#close file
close $fh;

#############
#print chart#
#############

#build the chart
my $chart = Chart::Clicker->new( format => "pdf", width => 1000 );
$chart->legend->visible(0);
my @x = ( $min .. $max );
my $funtf = Chart::Clicker::Data::Series->new( keys => \@x, values => \@counts );
$funtf->name("FunTF");
my @series;
push( @series, $funtf );

#compute x-axis ticks
my @xTicks;
my $stepSize = 200;
my $i        = -2400;
do {
	push( @xTicks, $i );
	$i += $stepSize;
} while ( $i < $max );

#format x-axis
$chart->get_context("default")->domain_axis;
my $xContext = $chart->get_context("default")->domain_axis;
$xContext->label("Position der Bindestellen relativ zum Transkriptionsstart");
$xContext->format("%d");
$xContext->range( Chart::Clicker::Data::Range->new( min => -1000, max => 200 ) );

#format y-axis
my $yContext = $chart->get_context("default")->range_axis;
$yContext->label("Transkriptionsfaktorbindestellen");
$yContext->format("%d");

#add markers
my $tss = Chart::Clicker::Data::Marker->new( key => 0, key2 => 2 );
my $range = Chart::Clicker::Data::Marker->new(
					key   => -1000,
					key2  => 50,
					color => Graphics::Color::RGB->new( red => 0, green => 0.5, blue => 0, alpha => 1 ),
					inside_color => Graphics::Color::RGB->new( red => 0, green => 1, blue => 0, alpha => 0.1 )
);
$chart->get_context("default")->add_marker($tss);
$chart->get_context("default")->add_marker($range);

#build the dataset
my $dataset = Chart::Clicker::Data::DataSet->new( series => \@series );

#add the dataset to the chart
$chart->add_to_datasets($dataset);

#colors
my $ca = Chart::Clicker::Drawing::ColorAllocator->new(
													   colors => [
																   Graphics::Color::RGB->new(
																							  red   => 0,
																							  green => 0,
																							  blue  => 1,
																							  alpha => 1.0
																   )
													   ]
);
$chart->color_allocator($ca);

#write the chart to a file
$chart->write_output("distribution_fun.pdf");
