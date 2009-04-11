#!/usr/bin/perl -w
use strict;
use Bio::SeqFeature::Generic;
use Bio::Location::Split;
use Bio::DB::Fasta;
use Bio::SeqIO;
use Bio::Location::Split;
use Bio::PrimarySeq;

use constant CodonSize => 3;
use Getopt::Long;
my $genetag = 'gene';
my $srctag  = 'SNAP';
my $munge = '';
my ($cdsout,$dbname);
GetOptions(
	   'm|munge:s'    => \$munge,
	   'db:s'         => \$dbname,
	   'cds:s'        => \$cdsout,
	   's|src:s'      => \$srctag);

my $db;
if( $dbname && $cdsout ) {
    $db = Bio::DB::Fasta->new($dbname);
}
my $out = Bio::Tools::GFF->new(-gff_version => 3);
my $outseq;
if( $cdsout ) {
    $cdsout = Bio::SeqIO->new(-format => 'fasta', -file => ">$cdsout");
}
my %genes;
my $seqid;
while(<>) {
    if (/^>(\S+)/){
	$seqid =  $1;
    } else {
	my ($type,$start,$end,$strand,$score,$fiveprimeO,
	    $threeprimeO,$frame,$group) = split;
	for ( $fiveprimeO, $threeprimeO,$frame ) {
	    $_ = 0 if $_ eq '.';
	}
	$strand = '-1' if $strand eq '-';
	$strand = '1' if $strand eq '+';
	push @{$genes{$group}}, Bio::SeqFeature::Generic->new
	    (-start  => $start,
	     -end    => $end,
	     -frame  => $frame,
	     -strand => $strand,
	     -score  => $score,
	     -primary_tag => $type,
	     -source_tag  => $srctag,
	     -seq_id  => $seqid,
	     -tag => {
		 'five_prime_overhang' => $fiveprimeO,
		 'three_prime_overhang' => $threeprimeO}
	     );
    }
}

while ( my ($gene,$exons) = each %genes ) {
    my $splitloc = Bio::Location::Split->new();
    my $strand;
    $gene = $munge . $gene if( $munge );
    my $seqname = $exons->[0]->seq_id;
    my ($phase,$next_phase) = (0,0);
    my ($cds);
    for my $exon ( map { pop @$_ }
		   sort {$a->[0] <=> $b->[0] } 
		   map { [$_->start * $_->strand, $_] }
		   @$exons ) {
	my $elen = $exon->length;
	if( defined $strand && $strand != $exon->strand) {
	    die("assertion: mixed strand exons!!! $gene\n") ;
	}
	$strand = $exon->strand;
	if( $exon->start !~ /\d+/ || $exon->end !~ /\d+/ ) {
	    die("$gene\n");
	}
	my ($five_overhang) = $exon->get_tag_values('five_prime_overhang');
	my ($three_overhang) = $exon->get_tag_values('three_prime_overhang');
	$splitloc->add_sub_Location($exon->location);
	$exon->remove_tag('group');
	$exon->add_tag_value('Parent', "Gene:$gene");
	$exon->frame($phase);
	$phase = ($elen % CodonSize) + $three_overhang - $five_overhang;
	$phase *= -1 if($phase < 0 );
	$phase -= 3 if $phase > 3;
	if( $db ) {
	    if( $exon->strand > 0 ) {
		$cds .= $db->seq($seqname, $exon->start => $exon->end);
	    } else {
		$cds .= $db->seq($seqname, $exon->end => $exon->start);
	    }
	}
    }
    if( $db ) {
	$cdsout->write_seq(Bio::PrimarySeq->new(-seq => $cds,
						-id  => $gene,
						-description => 
						sprintf("%s:%s",
							$seqname,
							$splitloc->to_FTstring())));
    }
    my $genef = Bio::SeqFeature::Generic->new
	(-start => $splitloc->start,
	 -end   => $splitloc->end,
	 -strand=> $exons->[0]->strand,
	 -seq_id=> $exons->[0]->seq_id,
	 -source_tag => $exons->[0]->source_tag,
	 -primary_tag=> $genetag,
	 -tag => { 'ID' => "Gene:$gene" }
	 );
    $out->write_feature($genef,@$exons);
}

    
