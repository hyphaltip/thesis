#!/usr/bin/perl -w
use strict;

while(<>) {
    if( /^>/ ) {
    } else { 
	my ($exon,$start,$end,$strand,undef,undef,undef,undef,$grp) = split;
	if( $strand eq '-' ) {
	    ($start,$end) = ( $end,$start);
	}
	$_ = join("\t", $exon,$start,$end,$grp)."\n";
    }
    print $_;
}
