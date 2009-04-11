#!/usr/bin/perl -w
use strict;
use Bio::DB::Fasta;
use Bio::SeqIO;
use File::Spec;

my $pepdir = 'pep_seqs';
my $outdir = 'ortholog_sequence_sets';

my $datfile = shift || die("need a dat file");
my $pref = shift || 'multi';

open my $in => $datfile || die("$datfile: $!");

mkdir($outdir) unless -d $outdir;

my $db = Bio::DB::Fasta->new($pepdir, -glob => "*.{fa,fas,fasta,seq,pep,cds,aa,nt,dna}");

my $paircount = 0;
open my $outfh => ">$pref\_match.tab" || die $!;
while(<$in>) {
    my @genes = split;

    my $out = new Bio::SeqIO(-file => ">$outdir/$pref\_$paircount.pep.fa",
			     -format=> 'fasta');
    for my $id ( @genes ) {
	my $seq = $db->get_Seq_by_acc($id);
	if( ! defined $seq ) {
	    my(undef,$nid) = split(/:/,$id);
	    if( ! ($seq = $db->get_Seq_by_acc($nid)) ) { 
		warn("No seq found for $id\n");
		next;
	    } 
	    $seq->display_id($id);
	}
	$out->write_seq($seq);
    }
    print $outfh join("\t", $paircount, @genes), "\n";
    $paircount++;
}
