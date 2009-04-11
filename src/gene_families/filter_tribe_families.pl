#!/usr/bin/perl -w
use Getopt::Long;
my $min_size = 10;

GetOptions(
	  's|size:s' => \$min_size,
	  );

while(<>) {
    if( /^\d+/ ) {
	my @line = split(/\t/,$_);
	my $desc = pop @line;
	my $id = shift @line;
	my $sum;
	$sum += $_ for ( @line );
	if( $sum > $min_size ) { 
	    print;
	}
    } else {
	s/_(BRD|GLEAN|TIGR)//g;
	print;
    }
}
