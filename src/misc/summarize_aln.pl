#!/usr/bin/perl -w
use strict;

use Getopt::Long;
use Bio::AlignIO;

my $format = "clustalw";

GetOptions(
	   'f|format:s' => \$format);

for my $file ( @ARGV ) {
    my $in = Bio::AlignIO->new(-format => $format,
			       -file   => $file);
    my @alns;
    while( my $aln = $in->next_aln ) {
	push @alns, $aln;
    }
    my $len = @alns;
    if( $len > 1 ) {
	my $i = 1;
	for my $aln ( @alns ) {
	    print join("\t", "$file.$i",
		       $aln->average_percentage_identity,
		       $aln->length,
		       ), "\n";
	}
	    $i++;
    } else {
	print join("\t", $file,
		   $alns[0]->average_percentage_identity,
		   $alns[0]->length,
		   ), "\n";
    }
}
