#!/usr/bin/perl -w
use strict;
use Bio::DB::Fasta;
use Bio::SeqIO;
use File::Spec;

my $pepdir = 'pep_seqs';
my $cdsdir = 'cds_seqs';
my $outdir = 'ortholog_sequence_sets';

my $datfile = shift || die("need a dat file");
my $pref = shift || 'multi';

open my $in => $datfile || die("$datfile: $!");

mkdir($outdir) unless -d $outdir;

my $cdsdb = Bio::DB::Fasta->new($cdsdir, -glob => "*.{fa,fas,fasta,seq,pep,cds,aa,nt,dna}");
my $pepdb = Bio::DB::Fasta->new($pepdir, -glob => "*.{fa,fas,fasta,seq,pep,cds,aa,nt,dna}");

my $paircount = 0;
open my $outfh => ">$pref\_match.tab" || die $!;
while(<$in>) {
    my @genes = split;

    my $out = new Bio::SeqIO(-file => ">$outdir/$pref\_$paircount.cds.fa",
			     -format=> 'fasta');
    for my $id ( @genes ) {
	my $cdsseq = $cdsdb->get_Seq_by_acc($id);
	
	if( ! defined $cdsseq ) {
	    my($prefix,$nid) = split(/:/,$id);
	    if( ! ($cdsseq = $cdsdb->get_Seq_by_acc($nid)) ) { 
		# now we have to do some remapping
		my ($pepid,$pepdesc) = split(/\s+/,$pepdb->header($id),2);
		my ($transcript_id) = ( $pepdesc =~ /transcript:(\S+)/);
		if( ! $transcript_id ) {
		    warn("cannot find id $id\n");
		    next;
		}
		$cdsseq = $cdsdb->get_Seq_by_acc("$prefix:$transcript_id");
		next unless $cdsseq;
		$cdsseq->description($pepdesc);
		$cdsseq->display_id($pepid);
	    } else {
		$cdsseq->description((split(/\s+/,$cdsdb->header($nid),2))[1]);
		$cdsseq->display_id($id);
	    }
	} else {
	    $cdsseq->description((split(/\s+/,$cdsdb->header($id),2))[1]);
	}
	$out->write_seq($cdsseq);
    }
    print $outfh join("\t", $paircount, @genes), "\n";
    $paircount++;
}
