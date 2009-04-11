#!/usr/bin/perl -w
use strict;
use File::Spec;
use Bio::AlignIO;
use Bio::DB::Fasta;
use Bio::SimpleAlign;
use Bio::LocatableSeq;
use Bio::Factory::FTLocationFactory;
use IO::String;
use Bio::TreeIO;
use Bio::SeqIO;
use Bio::PrimarySeq;

sub Bio::LocatableSeq::validate_seq {
    my ($self,$seqstr) = @_;
    if( ! defined $seqstr ){ $seqstr = $self->seq(); }
    return 0 unless( defined $seqstr); 
    if((CORE::length($seqstr) > 0) && ($seqstr !~ /^([A-Za-z\-\.\*\?\d]+)$/)) {
	$self->warn("seq doesn't validate, mismatch is " .
		    ($seqstr =~ /([^A-Za-z\-\.\*\?\d]+)/g));
	return 0;
    }
    return 1;
}
use constant CODONLEN => 3;
my %COLOR = ('A' => 'green',
	     'C' => 'blue',
	     'G' => 'orange',
	     'T' => 'red');
# amino color table from 
#http://info.bio.cmu.edu/Courses/BiochemMols/RasFrames/SHAPELY.HTM
my %AMINOPEPCOLOR = ('D'    => 'E60A0A',  # bright red
		'E'    => 'E60A0A',  # bright red
		'C'    => 'E6E600',  # yellow
		'M'    => 'E6E600',  # yellow
		'K'    => '145AFF',  # blue
		'R'    => '145AFF',  # blue
		'S'    => 'FA9600',  # orange
		'T'    => 'FA9600',  # orange
		'F'    => '3232AA',  # mid blue
		'Y'    => '3232AA',  # mid blue
		'N'    => '00DCDC',  # cyan
		'Q'    => '00DCDC',  # cyan
		'G'    => 'brown',  # light grey
		'L'    => '0F820F',  # green
		'V'    => '0F820F',  # green
		'I'    => '0F820F',  # green
		'A'    => 'C8C8C8',  # dark grey
		'W'    => 'B45AB4',  # ping
		'H'    => '8282D2',  # pale blue
		'P'    => 'DC9682',  # peach
		);

# shapely color table from 
#http://info.bio.cmu.edu/Courses/BiochemMols/RasFrames/SHAPELY.HTM
my %SHAPELYPEPCOLOR = ('D'    => 'A00042',  # dark red
		       'T'    => 'A00042',  # dark red
		       'E'    => '660000',  # red-brown
		       'C'    => 'FFFF70',  # bright yellow
		       'M'    => 'B8A042',  # dark yellow
		       'Y'    => 'B8A042',  # dark yellow
		       'K'    => '4747B8',  # blue
		       'R'    => '00007C',  # dark blue
		       'S'    => 'FF4C4C',  # orange
		       'Q'    => 'FF4C4C',  # orange
		       'F'    => '534C42',  # dark grey
		       'P'    => '534C42',  # dark grey
		       'W'    => '534C42',  # dark grey
		       'N'    => 'FF7C70',  # peach
		       'G'    => 'C8C8C8',  # grey
		       'V'    => 'C8C8C8',  # grey
		       'I'    => '004C00',  # dark green
		       'L'    => '455E45',  # grey-green
		       'A'    => '8CFF8C',  # light green
		       'H'    => '7070FF',  # pale blue
		       );		
my %PEPCOLOR = %AMINOPEPCOLOR;
my $locfactory = Bio::Factory::FTLocationFactory->new(-verbose => 0);
my $force = 0;
my $idx = Bio::DB::Fasta->new('pep_seqs', -glob => "*.{fa,fasta,fast,FA,FASTA,FAST,dna,pep,aa,peps,cds}");

my $file = shift || die;
my $alnformat = 'fasta';

my $db;# = Bio::DB::Fasta->new('intronstats');

my $outfhrpt;

my $in = Bio::AlignIO->new(-format => $alnformat,
			   -file   => $file);

my ($outfh,$outfh_map);
my ($stem) = $file;
$file =~ s/\.(fasaln|clustal|aln)$/;
open($outfh, ">$stem.intronmap.html") || die $!;
open($outfh_map, ">$stem.intronmap.fasaln") || die $!;
print $outfh "<html>\n";
if( my $aln = $in->next_aln ) {
    $aln->set_displayname_flat(1);
    my $final_seqs = Bio::SimpleAlign->new;    
    my (%add_columnsp,%add_columns,%colmap);

    for my $g ( $aln->each_seq ) {
	my $id = $g->display_name;
	my $header = $idx->header($id);
	my @d = split(/\s+/,$header,3);
	shift @d;		# get rid of ID/name
	my ($chrom,$loc);
	for ( @d ) {
	    if(/(\S+):(\S+)/ ) {
		($chrom,$loc) = ($1,$2);
		last;
	    }
	}
	if( ! $loc ) {
	    warn("cannot find loc in ",$idx->header($id),"\n");
	}
	$loc = $locfactory->from_string($loc);
	$loc->seq_id($chrom);
	my $runninglen = 0;
	my @intron_positions;
	my $offset = 0;
	my $poffset = 0;
	my $ct = 1;
	for my $subloc ( sort { $a->start * $a->strand <=> 
				    $b->start * $b->strand } 
			 $loc->each_Location ) {
	    if( $runninglen && $subloc->length > 3 ) { 
		my $intronphase = $runninglen % CODONLEN;
		# warn("$runninglen ",$g->length, "\n");
		next if( $runninglen > $g->length );		    
		my $col;
		eval { 
		    $col= $g->column_from_residue_number($runninglen);
		};
		if( $@ || ! defined $col ) { 
		    warn("no column for $runninglen $id\n");
		} else { 
		    $colmap{$id}->{$col} = $ct++;
		    push @intron_positions, [$col, $intronphase];
		}
	    }
	    $runninglen += $subloc->length;
	}
	my $pseq = $g->seq();

	for my $pos ( @intron_positions ) {
	    next unless defined $pos;
	    $add_columns{$pos->[0]}->{$id} = $pos->[1];
	    my $ppos = int((($pos->[0] -1)/ CODONLEN)) +1;
	    $add_columnsp{$ppos}->{$id} = $pos->[1];
	    substr($pseq, $ppos+$poffset,0,$pos->[1]);
	    $poffset++;
	    $offset++;
	}

	my $ps = Bio::LocatableSeq->new(-verbose => -1,
					-display_id => $id);
	#$ps->no_validate(1);
	$ps->alphabet('protein');
	$ps->seq($pseq);
	$final_seqs->add_seq($ps);
    }	
    $final_seqs->set_displayname_flat(1);
    my %nmset;	
    $final_seqs->set_displayname_flat(1);
    for my $seq ( $final_seqs->each_seq ) {
	my $pseq = $seq->seq();
	my $offset = 0;
	for my $col ( sort {$a <=> $b } keys %add_columnsp ) {
	    if( ! defined $add_columnsp{$col}->{$seq->display_id} ) {
		substr($pseq, $col+$offset,0,'-');
		$add_columnsp{$col}->{$seq->display_id} = '-';
	    }
	    $offset++;
	}
	$seq->seq($pseq);
    }

    # pep seqs
    my $str = IO::String->new();
    my $out = Bio::AlignIO->new(-format => 'clustalw',
				-linelength => 100,
				-fh     => $str);
    $out->write_aln($final_seqs);
    $out->close();
    my $str2 = IO::String->new(${$str->string_ref});
    print $outfh "<pre>\n";
    $out = Bio::AlignIO->new(-format => 'fasta',
			     -fh     => $outfh_map);
    $out->write_aln($final_seqs);
    while(<$str2>) {
	if( /^CLUSTAL|MUSCLE/ ) {
	} elsif( /^(\s*\S+\s+)(\S+)(\s*\d*\s*)$/ox ) {
	    my ($pre,$seq,$post) = ($1,$2,$3);
	    # add 'bold' tags
	    $seq =~ s/(\d+)/<b>$1<\/b>/g;
	    $seq =~ s/([DTECMYKRSQFPWNGVILAH])/<font color="$PEPCOLOR{$1}">$1<\/font>/g;
	    $_ = join('',$pre,$seq,$post);
	}
	print $outfh $_;
    }
    print $outfh "<\/pre>\n<\/html>\n";
} else { 
    warn("no aln for $file\n");
}
