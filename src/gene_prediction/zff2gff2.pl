#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
my $tag = 'Genomewise';
GetOptions(
	   't|tag:s' => \$tag,
	   );
my $seqid;
my (%exons,%loc,@order);
while(<>) {
    if(/>(\S+)/ ) {
	$seqid = $1;
    } else {
	my ($exon,$start,$end,$group) = split;

	my $strand;
	if( $start < $end ) {
	    $strand = '+';
	} else {
	    ($start,$end,$strand) = ($end,$start,'-');
	}
	my ($min,$max);
	if( ! defined $loc{$group} ) {
	    $loc{$group} = { 'min'    => $start, 
			     'max'    => $end, 
			     'strand' => $strand,
			     'seqid'  => $seqid,
			 };
	    push @order, $group;
	} else {
	    if( $loc{$group}->{min} > $start ) {
		$loc{$group}->{min} = $start;
	    }
	    if( $loc{$group}->{max} < $end ) {
		$loc{$group}->{max} = $end;
	    }
	}
	push @{$exons{$group}}, join("\t", ($seqid, $tag,$exon,
					    $start,$end,'.',$strand,'.',
					    "GenePrediction $group"));
    }
}

for my $gene ( @order ) {
    print join("\t", $loc{$gene}->{'seqid'}, $tag, 'mRNA', 
	       $loc{$gene}->{'min'}, $loc{$gene}->{'max'},
	       '.', $loc{$gene}->{'strand'}, "GenePrediction $gene"), "\n";
    print join("\n", @{$exons{$gene}}), "\n";
}
