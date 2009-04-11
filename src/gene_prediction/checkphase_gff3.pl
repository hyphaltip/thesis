#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Bio::Tools::GFF;
use List::Util qw(min max);
my $inf = 3;
my $outf = 3;
GetOptions('if:s'     => \$inf,
	   'of:s'     => \$outf);


my $in = Bio::Tools::GFF->new(-fh => \*ARGV, -gff_version => $inf);
my $out = Bio::Tools::GFF->new(-fh => \*STDOUT, -gff_version => $outf);

my @gene;
my $last;
while (my $f = $in->next_feature) {
    if( $f->primary_tag =~ /(gene)|mRNA|CDS|(Exon|Einit|Eterm|Esngl|internal|terminal|initial|single-exon|internal-exon|terminal-exon|final-exon|initial-exon)/ ) {
	my $name;
	for my $p ( qw(ID Parent Gene Target Sequence Protein EST Group)){
	    if( $f->has_tag($p) ) {
		($name) = $f->get_tag_values($p);
		last;
	    }
	}
	if( $1 ) {
	    $f->primary_tag('mRNA');
	} elsif( $2 ) {
	    $f->primary_tag('CDS');
	    $f->add_tag_value('ExonType', $2);
	}
	if(defined($last) && $last ne $name) {
	    calc_phase(\@gene,$name);
	    $out->write_feature(@gene);
	    @gene = ();
	}
	$last = $name;
	push @gene, $f;
    } else {
	$out->write_feature($f);
    }
}

calc_phase(\@gene,$last) if @gene;
$out->write_feature(@gene) if @gene;

sub calc_phase {
    my $gene = shift @_;
    my $name = shift @_;
    my ($mRNA) = grep { $_->primary_tag eq "mRNA" } @{$gene};
    my ($vulgar) = $mRNA->get_tag_values("vulgar") if $mRNA;
    if ($vulgar && $vulgar =~ m/ F \d+ \d+/) {
      warn $name, " has frameshifts, cannot (yet) be processed\n";
      @$gene = ();
      return;
    }

    my @gene = sort { $a->start <=> $b->start } grep { $_->primary_tag eq "CDS" } @{$gene};
    unless (@gene) {
      warn $gene->[0]->gff_string . " doesn't have any CDS subfeatures\n";
      return;
    }

    @gene = reverse @gene if $gene[0]->strand() == -1;

    my $cdslen = $gene[0]->frame;

    if($gene[0]->has_tag("ExonType")) {
	my ($type) = $gene[0]->get_tag_values("ExonType");
	if ( ($type =~ /Internal|Exon|Esngl/ ) ||
	     ($type =~ /Terminal|Esngl|Eterm/ && @gene == 1)) {
	    my $len = 0;
	    for my $f (@gene) {
		$len += $f->length;
	    }
	    my $phase = $len % 3;
	    $cdslen = $phase == 1 ? 2 : $phase == 2 ? 1 : 0;
	} else {
	    $cdslen = 0;
	}
	for (@gene) { $_->remove_tag("ExonType"); }
    } elsif ($cdslen eq ".") {
	$cdslen = 0;
    } else {
	$cdslen = (3 - $cdslen) % 3;
    }

    for my $f (@gene) {
	my $phase = $cdslen % 3;
	$f->frame($phase == 1 ? 2 : $phase == 2 ? 1 : 0);
	$cdslen += $f->length;
    }

    if (@gene == @$gene) {
	my $mRNA = $gene[0]->new;
	$mRNA->seq_id($gene[0]->seq_id);
	$mRNA->start(min(map { $_->start } @gene));
	$mRNA->end(max(map { $_->end } @gene));
	$mRNA->primary_tag("mRNA");
	$mRNA->source_tag($gene[0]->source_tag);
	$mRNA->strand($gene[0]->strand);
	$mRNA->add_tag_value('ID', $gene[0]->get_tag_values("Parent") );
	unshift @$gene, $mRNA;
    }

}
