#!/usr/bin/perl -w
use strict;
use Bio::Tools::RepeatMasker;
use Bio::Tools::GFF;
my $parser = Bio::Tools::RepeatMasker->new(-file => shift);
my $out = Bio::Tools::GFF->new(-gff_version => 3);
while( my $result = $parser->next_result ) {
    $result->feature2->primary_tag('similarity');
    $out->write_feature($result);
}
