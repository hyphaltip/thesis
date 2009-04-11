#!/usr/bin/perl -w
use strict;
use Getopt::Long;
my $type = 'target';
my  $ptag = 'HSP';
my  $stag = 'BLASTN';

GetOptions(
	   't|type:s' => \$type,
	   'p|ptag:s' => \$ptag,
	   's|stag:s' => \$stag,);

while(<>) {
    my ($q,$h,$pid,$hsp_len,$gap,$mm,$qstart,$qend,$hstart,$hend,
	$e_value,$score) = split;
    my $st = 1;
    if( $hstart > $hend) {
	($hstart,$hend,$st) = ($hend,$hstart,-1);
    }
    elsif( $qstart > $qend ) {
	($qstart,$qend,$st) = ($qend,$qstart,-1);
    }
    if( $type =~ /^q/i ) {
	print join("\t", ( $q,
			   $stag,
			   $ptag,
			   $qstart,
			   $qend,
			   $score,
			   $st,
			   '.',
			   "Target=$h,$hstart,$hend")), "\n";		   
    } elsif( $type =~ /^(h|t)/i ) {
	
	print join("\t",
		   ( $h,
		     $stag,
		     $ptag,
		     $hstart,
		     $hend,
		     $score,
		     $st,
		     '.',
		     "Target=$q,$qstart,$qend")), "\n";
    }
}
