#!/usr/bin/perl -w
use strict;
use Bio::Tools::GFF;
use Bio::SeqFeature::Generic;
my $out= Bio::Tools::GFF->new(-gff_version => 3);
my $header = <>;
chomp($header);
my @header = split(/\t/,$header);
my %seen;

while(<>) {
    chomp;
    my ($chrom,$strand,$geneid,$desc,$pepid,$exonid,
	$cds_start,$cds_end) = split(/\t/,$_);
    
    next unless $cds_start && $cds_end; # UTR if no CDS positions
    if( $cds_start > $cds_end ) { 
	if(abs($cds_start-$cds_end) <= 3) {
	    $cds_start = $cds_end;
	}
    } 
    my $exon = Bio::SeqFeature::Generic->new
	(-seq_id      => $chrom,
	 -primary_tag => 'CDS',
	 -source_tag  => 'EnsEMBL',
	 -strand      => $strand,
	 -start       => $cds_start,
	 -end         => $cds_end,
	 -tag         => { 
	     'Parent'  => $pepid,
	 });
    $out->write_feature($exon);
    if( ! defined $seen{$pepid} ) {
	$seen{$pepid} = Bio::SeqFeature::Generic->new
	    (-start => $cds_start,
	     -end   => $cds_end,
	     -strand=> $strand,
	     -seq_id=> $chrom,
	     -primary_tag => 'gene',
	     -source_tag  => 'EnsEMBL',
	     -tag   => { 'ID' => $pepid,
			 'Note'=> $desc,
			 'GeneName' => $geneid,
		     }
	     );
    } else {
	my ($s,$e) = ($seen{$pepid}->start,
		      $seen{$pepid}->end);
	$seen{$pepid}->start($cds_start) if( $s > $cds_start);
	$seen{$pepid}->end($cds_end) if( $e < $cds_end);
    }
}

for my $gene ( values %seen ) {
    $out->write_feature($gene);
}
