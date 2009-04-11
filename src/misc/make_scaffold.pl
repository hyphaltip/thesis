#!/usr/bin/perl -w

=head2 NAME

makE_scaffold - make a scaffold file for Gbrowse

=cut

use strict;

use Bio::SeqIO;
use Bio::Tools::GFF;
use Bio::SeqFeature::Generic;

my $in = new Bio::SeqIO(-file => shift);
my $version = shift ||3;
my $out = new Bio::Tools::GFF(-gff_version=>$version);

while( my $seq = $in->next_seq ) {
    my $feature = new Bio::SeqFeature::Generic(-start => 1,
					       -end   => $seq->length,
					       -strand=>1,
					       -seq_id=> $seq->display_id,
					       -source_tag => 'chromosome',
					       -primary_tag  => 'scaffold',
					       -tag => { 
						   'ID'       => $seq->display_id,
						   'Sequence' => $seq->display_id 
					       }
					       );
    $out->write_feature($feature);
}
