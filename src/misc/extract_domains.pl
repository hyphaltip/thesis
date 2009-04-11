#!/usr/bin/perl -w
use strict;
use Bio::SearchIO;
use Bio::DB::Fasta;
use Bio::SeqIO;
use Bio::PrimarySeq;
my $file = shift || die "need a file";
my $dbname = shift || "peps";

my $db = Bio::DB::Fasta->new($dbname, -glob => "*.{fa,fasta,fast,FA,FASTA,FAST,dna,pep,aa,peps,cds}");

my $in = Bio::SearchIO->new(-format => 'hmmer',
			    -file   => $file);

my $out = Bio::SeqIO->new(-format => 'fasta');
my %domains;
while( my $r = $in->next_result ) {
    while( my $h = $r->next_hit ) {
	my $hname = $h->name;
#	next if( $hname =~ /^(scer_yjm789|scer_rm11|spar|skud|smik|sklu|sbay|scas|cbri|cele|dmel|ddis)/);
	while( my $hsp = $h->next_hsp ) {
	    next if $hsp->evalue > 1e-10;
	    my ($hstart,$hend) = ($hsp->hit->start,
				  $hsp->hit->end);
	    push @{$domains{$hname}}, [ $hstart,$hend, $hsp->evalue,
			     $db->seq($hname, $hstart => $hend)];
	}
    }
}
my %seen;
for my $hname ( keys %domains ) {
    my @domains = @{$domains{$hname}};
    $hname =~ s/:/_/;
    $hname =~ s/\.p\d+$//;
    $hname =~ s/[\.\-]/_/g;
    if($hname =~ s/^(\w+)\_(\1)/$1/ ) {
	
    }
    $hname =~ s/GLEAN_(gz2?_)?//;
    next if $seen{$hname}++;
    my $i = 1;
    for my $s ( sort { $a->[0] <=> $b->[0] } @domains ) {
	$out->write_seq(Bio::PrimarySeq->new
			(-seq => pop @$s,
			 -id  => "$hname\_A$i",
			 -desc => sprintf("protein_location:%d..%d evalue:%s",@$s)));
	$i++;
    }
}
