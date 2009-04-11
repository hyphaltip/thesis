#!/usr/bin/perl -w
use strict;
use Bio::Tools::GFF;
use Bio::SeqFeature::Generic;
use Data::Dumper;
use Bio::Location::Split;
use constant CodonSize => 3;
use Getopt::Long;
my $genetag = 'gene';
my $srctag  = 'SNAP';
my $munge = '';
GetOptions(
	   'm|munge:s'    => \$munge,
	   's|src:s'      => \$srctag);

my $in = Bio::Tools::GFF->new(-gff_version => 1,
			      -fh          => \*ARGV );
my $out = Bio::Tools::GFF->new(-gff_version => 3);

my %genes;

while( my $f = $in->next_feature ) {
    my ($group) = $f->get_tag_values('group');
    warn($f->frame,"\n");
    push @{$genes{$group}}, $f;
}

while ( my ($gene,$exons) = each %genes ) {
    my $splitloc = Bio::Location::Split->new();
    my $strand;
    $gene = $munge . $gene if( $munge );
	
    my ($runninglength,$phase) = (0,0);
    for my $exon ( @$exons ) {
	my $elen = $exon->length;
	if( defined $strand && $strand != $exon->strand) {
	    die("assertion: mixed strand exons!!! $gene\n") ;
	}
	$strand = $exon->strand;
	if( $exon->start !~ /\d+/ || $exon->end !~ /\d+/ ) {
	    die("$gene\n");
	}
	$splitloc->add_sub_Location($exon->location);
	$exon->remove_tag('group');
	$exon->add_tag_value('Parent', $gene);
	$exon->frame();
	$phase = ($runninglength += $elen) % CodonSize;
	$exon->source_tag($srctag);
    }
    my $genef = Bio::SeqFeature::Generic->new
	(-start => $splitloc->start,
	 -end   => $splitloc->end,
	 -strand=> $exons->[0]->strand,
	 -seq_id=> $exons->[0]->seq_id,
	 -source_tag => $exons->[0]->source_tag,
	 -primary_tag=> $genetag,
	 -tag => { 'ID' => $gene }
	 );
    $out->write_feature($genef,@$exons);
}

    
