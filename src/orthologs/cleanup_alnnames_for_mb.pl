#!/usr/bin/perl -w
use strict;
use Bio::AlignIO;
use Getopt::Long;
use List::Util qw(sum);
my $file = shift || die;

my $min_fraction = 0.25;

my ($iformat,$oformat) = qw(fasta nexus);
my $trim = 0;
GetOptions(
	   'if|iformat|informat|format:s' => \$iformat,
	   'of|outformat|oformat:s'       => \$oformat,
	   't|trim'                       => \$trim,
	   );


my $in = Bio::AlignIO->new(-format => $iformat,
			   -file   => $file);
my $out = Bio::AlignIO->new(-format => $oformat,
			    -show_endblock => 0,
			    -show_symbols  => 0);


while( my $aln  = $in->next_aln ) {
    my $newaln = $aln->new;
    my @cols;
    if( $trim ) {
	my $matrix = $aln->gap_col_matrix;
	my %colstoremove;
	my $i = 0;
	for my $cl ( @{$matrix} ) {
	    my $num = scalar keys %$cl;
	    my $sum = sum values %$cl;
	    if( ($sum / $num) > $min_fraction ) {
		$colstoremove{$i}++;
	    }
	    $i++;
	}
	@cols = sort { $b <=> $a} keys %colstoremove;
    }

    my %seen;    
    for my $seq ( $aln->each_seq ) {
	my $id= $seq->display_id;
	my $newid = $seq->display_id;

	if( $trim) { 
	    my $str = $seq->seq;
	    for my $col ( @cols  ) {
		substr($str,$col,1,''); # remove the column
	    }
	    $str =~ s/X/-/g;
	    $seq->seq($str);
	}
	$newid =~ s/:GLEAN_(gz\d*_)?/_/;	
	$newid =~ s/_[^:]+:/_/;
	$newid =~ tr/-/_/;
	$newid =~ s/\.p\d+$//;
	$newid =~ s/\.\d+$//;
	$newid =~ s/_GLEAN//g;
	$newid =~ s/[.:]/_/g;
	next if $seen{$newid}++;
	$seq->display_id($newid);
	
	$newaln->add_seq($seq);
    }
    $out->write_aln($newaln);
}
