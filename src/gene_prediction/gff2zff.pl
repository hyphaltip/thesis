#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Bio::Tools::GFF;
use Bio::SeqFeature::Generic;
use Getopt::Long;
my $field = 'CDS';
my $species;
GetOptions(
	   'sp|species:s' => \$species,
	   'f|field:s'    => \$field,
	   );
die("Must provide species\n") unless defined $species;
my $zff;
open($zff => ">$species.ann") || die $!;
my $in = Bio::Tools::GFF->new(-gff_version => 3, -file => shift);
my %groups;
while( my $f = $in->next_feature ) {
    next unless $f->primary_tag eq $field;
    my $group;
    for my $tag ( qw(ID Parent Gene Transcript Protein Target) ) {
	if( $f->has_tag($tag) ) {
	    ($group) = $f->get_tag_values($tag);
	    last;
	}
    }
    if( ! $group ) { 
	warn("cannot find group for ",$f->gff_string, "\n");
	exit;
    }
    push @{$groups{$f->seq_id}->{$group}}, $f;
}

while( my ($seqid,$features) = each %groups ) {
    my @line;
    while(  my ($gene,$d) = each %{$features} ) {
	my $d = $features->{$gene};
	my $i = 1;
	my $count = @$d;
	for my $df ( map { $_->[0] }
		     sort { $a->[1] <=> $b->[1]} 
		     map { [$_, $_->start * $_->strand ] }
		     @$d ) {

	    my $type = 'Exon';
	    if( $count == 1 ) {
		$type = 'Esngl';
	    } elsif( $i == 1 ) {
		$type = 'Einit';
	    } elsif( $i == $count ) {
		$type = 'Eterm';
	    }

	    my ($l,$r) = ($df->start, $df->end);
	    if( $df->strand < 0 ) { 
		($l,$r) = ( $r, $l);
	    }
	    push @line, [$df->start,
			 sprintf("%s\t%d\t%d\t%s\n",$type,$l,$r, $gene)];
	    $i++;
	}	
    }
    print $zff ">$seqid\n";
    for my $l ( sort { $a->[0] <=> $b->[0] } @line ) {
	print $zff $l->[1];
    }
}
