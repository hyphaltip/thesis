#!/usr/bin/perl -w
# Author Jason Stajich <jason@cgt.mc.duke.edu>

use strict;
use Bio::SeqFeature::Generic;
use Bio::Tools::GFF;
use Getopt::Long;
my $version = 3;
GetOptions('v|version:s'   => \$version);
my $out = new Bio::Tools::GFF(-gff_version => $version);

my %genes;
my $header = <>;
chomp($header);
my $i = 0;
$header =~ s/^\#//;
my @h = split(',',$header);
my %header_fields = map { $_ => $i++ } @h;
$i=0;
$,=';';
while(<>) {
    my $cp = $_;
    my $orig = $_;
    # remap so we can deal with commas in quotes, etc
    # some boundary cases still exist I suspect
    while( $cp =~ /(\,\"[^\"]+)\,([^\"]+\"\,)/ ) {
	$cp =~ s/(\,\"[^\"]+)\,([^\"]+\"\,)/$1:::$2/;
    }
    chomp;
    my @fields = map { s/:::/,/; s/\"//g; $_;} split(',',$cp);
    my $i = 0;
    my $contig = $fields[ $header_fields{CONTIG_NAME} ];
    my @parts = split(/\s+/,$contig);
    unless( $parts[3] ) {
        print STDERR "You probably need to edit the HEADER of the whitehead file.  Default whitehead column names are missing some columns ($contig)\n$orig\n$_\n";
        last;
    }
    $contig = sprintf("Contig_%s",$parts[3]);
     
    my $group = $fields[$header_fields{GENE_LOCUS}];
    my $note = $fields[$header_fields{FULL_GENE_NAME}];
    $note =~ s/\s+$//;

    unless( defined $genes{$group} ) {
        push @{ $genes{$group} }, Bio::SeqFeature::Generic->new
            (
             -start  => $fields[$header_fields{GENE_START}],
             -end    => $fields[$header_fields{GENE_STOP}],
             -strand => $fields[$header_fields{GENE_STRAND}],
             -seq_id => $contig,
             -primary_tag => 'mRNA',
             -source_tag  => 'broad',
             -tag => {
		 'ID'   => "Gene:$group",
		 'Note' => $note ,
	     },
             );
    }
    push @{ $genes{$group} },
    Bio::SeqFeature::Generic->new
        ( -start  => $fields[$header_fields{EXON_START}],
          -end    => $fields[$header_fields{EXON_STOP}],
          -strand => $fields[$header_fields{EXON_STRAND}],
          -seq_id => $contig,
          -primary_tag => 'CDS',
          -source_tag  => 'broad',
          -tag   => { 
	      'Parent' => "Gene:$group" }
          );
}

foreach my $gene  ( values %genes ) {
    foreach my $exon ( @$gene ) {
        $out->write_feature($exon);
    }
}

