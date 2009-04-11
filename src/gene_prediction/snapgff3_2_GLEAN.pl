#!/usr/bin/perl -w
while(<>) {
    chomp;
    next if /^\#/;
    my ($seqid,$src,$type,
	$start,$end,$score,
	$strand,$phase,
	$group) = split(/\s+/,$_,9);
    my @groups;
    if( $type eq 'gene' ) { $type = 'mRNA' }
    for my $field ( split(/;/,$group) ) {
	my ($key,$vals) = split(/=/,$field,2);
	my @vals = split(/,/,$vals);
	if( $key =~ /Parent|ID/i) { 
	    $key = 'GenePrediction';
	}
	if( $key eq 'GenePrediction' ) {
	    unshift @groups, sprintf("%s %s",$key,join(' ', @vals));
	} else {
	    push @groups, sprintf("%s %s",$key,join(' ', @vals));
	}
    }	
    print join("\t",
	       $seqid,$src,$type,$start,$end,$score,
	       $strand,$phase,
	       join("; ", @groups)), "\n";
}
