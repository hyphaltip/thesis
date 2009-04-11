#!/usr/bin/perl -w
use strict;
my %ac;

# run this on Pfam_ls or Pfam_fs database to get a description database of the form
# >NAME|ACC DESCRIPTION
# which is used in the map2genome.pl script to map Pfam hits (from proteins) back onto the
# genome suitable for loading into Gbrowse
while(<>) {
    chomp;
    if(/^(NAME|ACC|DESC)\s+(.+)/ ) {
	my ($type,$d) = ($1,$2);
	$ac{$type} = $d;
	if( $type eq 'DESC' ) {
	    printf ">%s|%s %s\n",map { $ac{$_} } qw(NAME ACC DESC);
	    %ac = ();
	}
    }
}
