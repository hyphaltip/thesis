#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Bio::DB::Fasta;
use Bio::PrimarySeq;
use Getopt::Long;
use File::Spec;
use File::Basename;
my ($db,$Max_members,$Min_members) = (undef,50,5);
my $printseq = 1;
GetOptions('d|db:s'  => \$db,
	   'max:i'   => \$Max_members,
	   'min:i'   => \$Min_members,
           's|seq!'  => \$printseq,
	   );
my $dbh;
if( $printseq) {
    unless( $db ) { 
	die( "did not provide a valid db name, please provide a -db opption\n");
    }
    $dbh = Bio::DB::Fasta->new($db, -glob => "*.{fa,aa,cds,pep,fasta,fas,fast,FA,FAST,FASTA,FAS,nt,dna}");
}

for my $file ( @ARGV ) {
    my $countfh;
    my $spcountfh;
    my $fileold = $file;
    my ($prefix) = fileparse($fileold, qw(.families .fam));
    $prefix .= ".d";
    open($countfh, ">$file.count") || die $!;
    open($spcountfh, ">$file.spcount") || die $!;
    mkdir($prefix) unless -d $prefix;

    my $in;
    open($in, File::Spec->catfile($file)) || die $!;
    my @fams;
    my @spfam;
    my %spnames;
    while(<$in>) {
	my ($fam,$id) = split;
	my ($sp,$gn) = split(/:/,$id);
	push @{$fams[$fam]}, $id;
	$spfam[$fam]->{$sp}++;
	$spnames{$sp} = 1 unless $spnames{$sp};
    }
    my @spec = sort keys %spnames;
    my $i = 0;
    print $spcountfh join("\t", 'FAM', @spec),"\n";
    for my $fam ( @fams ) {
	my $count = scalar @$fam;
	print $countfh join("\t", $i, $count), "\n";
	print $spcountfh join("\t",$i, map {$spfam[$i]->{$_} || 0} 
			      @spec ),"\n";
	if( $printseq )  {
	    if( ($Max_members && $count < $Max_members) && ($Min_members && $count > $Min_members)) {
		my $out = Bio::SeqIO->new(-format => 'fasta',
					  -file   => ">$prefix/$i.pep.fa");
		for my $id ( @$fam ) {
		    my $seq = $dbh->seq($id);
		    my ($idpart,$desc) = split(/\s+/,$dbh->header($id),2);
		    if( $seq ) { 
			$out->write_seq(Bio::PrimarySeq->new
					(-seq         => $seq,
					 -display_id  => $id,
					 -description => $desc));
		    } else {
			warn("cannot find $id in db $db\n");
		    }
		}	    
	    }
        }
	$i++;
    }
}

