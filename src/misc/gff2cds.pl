#!/usr/bin/perl -w
use strict;
use Bio::DB::Fasta;
use Getopt::Long;
use Bio::SeqIO;
use Bio::PrimarySeq;
use Bio::Location::Split;
use Bio::Location::Simple;
my $db_dir = '/scratch/jes12/fa';
my $gff_version = 3;
GetOptions(
	   'gffv:i'            => $gff_version,
	   'd|db|database:s'   => \$db_dir,
	   );
my $db = Bio::DB::Fasta->new($db_dir);
my %transcripts;
while(<>) {
    chomp;
    my $record = $_;
    my ($seqid,$src,$type,$start,$end,$score,
	$strand,$frame,$group) = split(/\s+/,$record,9);
    next unless $type =~ /CDS|ORF|Exon|Einit|Eterm|Esngl/i; 
    my ($primary_id,%grp);
    
    if( $gff_version ==3 ) {
	my @groups = split(/\s*;\s*/,$group);
	for ( @groups ) {
	    my ($key,$vals) = split(/=/,$_,2);
	    my @vals= split(/,/,$vals);
	    if( exists $grp{$key} )  {
		warn("What! double defined group tag name $key for $record\n");
	    }
	    $grp{uc $key} = [@vals];
	}	
    } else {
	die("unsupported gff version");
    }
    for my $id ( qw(ID PARENT GENE TRANSCRIPT TARGET PROTEIN EST ) ) {
	if( $grp{$id} ) {
	    ($primary_id) = @{$grp{$id}};
	    last;
	}
    }
    $strand = 1  if $strand eq '+';
    $strand = -1 if $strand eq '-';

    if( ! $primary_id ) {
	warn("no good ids found in ", join(',',keys %grp),"\n");
    }
    push @{$transcripts{$primary_id}}, [$seqid,$start,$end,$strand];
}
for my $v(  values %transcripts ) {
    @$v = sort { $a->[1] * $a->[3] <=> $b->[1] * $b->[3] } @$v;
}

my @order = ( map { $_->[2] } 
	      sort { $a->[0] cmp $b->[0] || $a->[1] <=> $b->[1] } 
	      map{
		  [$transcripts{$_}->[0]->[0], 
		    $transcripts{$_}->[0]->[1] * $transcripts{$_}->[0]->[3],
		    $_ ] }
	      keys %transcripts);
my $out = Bio::SeqIO->new(-format => 'fasta');
for my $transcript ( @order ) {
    my $cdsseq = '';
    my $splitloc = Bio::Location::Split->new();
    my $seqid;
    for my $cds ( @{ $transcripts{$transcript} } ) {
	$seqid = $cds->[0] unless $seqid;
	if( $cds->[3] < 0 ) {
	    $cdsseq .= $db->seq($cds->[0], $cds->[2] => $cds->[1]);
	} else {
	    $cdsseq .= $db->seq($cds->[0], $cds->[1] => $cds->[2]);
	}
	$splitloc->add_sub_Location(Bio::Location::Simple->new
				    (-start => $cds->[1],
				     -end   => $cds->[2],
				     -strand=> $cds->[3]));
    }
    if( ! $cdsseq ) {
	die("No sequence found for $seqid, $transcript\n");
    } else {
	$out->write_seq(Bio::PrimarySeq->new
			(-display_id => $transcript,
			 -description => 
			 sprintf("%s:%s cdslen=%d",
				 $seqid, $splitloc->to_FTstring(),
				 length($cdsseq)),
			 -seq => $cdsseq));
    }
					 
}
