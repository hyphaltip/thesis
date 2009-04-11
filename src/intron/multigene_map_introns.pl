#!/usr/bin/perl -w
use strict;
use File::Spec;
use Getopt::Long;
use Bio::AlignIO;
use Bio::DB::Fasta;
use Bio::SimpleAlign;
use Bio::LocatableSeq;
use Bio::Factory::FTLocationFactory;
use IO::String;
use Bio::TreeIO;
use Bio::SeqIO;
use Bio::PrimarySeq;

my $html = 0;

GetOptions(
	   'html!' => \$html,
	   );
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
my $idx = Bio::DB::Fasta->new('cds_seqs', 
			      -glob => "*.{fa,fasta,fast,FA,FASTA,FAST,dna,pep,aa,peps,cds}");
my $pepidx = Bio::DB::Fasta->new('pep_seqs', 
				 -glob => "*.{fa,fasta,fast,FA,FASTA,FAST,dna,pep,aa,peps,cds}");

my $dir = shift || 'ortholog_sequence_sets-aln';
my $match = shift || 'multi_match.tab';

my %gffidx; # index gff locations

my $db;# = Bio::DB::Fasta->new('intronstats');

open(MATCH, $match) || die "$match: $!";
my @num2seqs;
while(<MATCH>) {
    my ($num,@seqnames) = split;
    for my $g ( @seqnames ) { 
	my ($sp,$gn) = split(/:/,$g,2);
	$num2seqs[$num]->{$sp} = $gn;
    }
}

opendir(DIR, $dir) || die $!;
for my $file ( readdir(DIR) ) { 
    next unless $file =~ /multi\_(\d+)\.cds\.nex$/;
    my ($num) = $1;
    next unless $num2seqs[$num];
    my $in = Bio::AlignIO->new(-format => 'nexus',
			       -file   => File::Spec->catfile($dir,$file));
    
    my ($outfhp,$outfh,$outfhp_map,$outfh_map);
    if( $html ) {
	next if -e "$dir/multi_$num.cds.intronmap.html" && ! $force;
	open($outfh, ">$dir/multi_$num.cds.intronmap.html") || die $!;
	open($outfhp, ">$dir/multi_$num.pep.intronmap.html") || die $!;
	print $outfh "<html>\n";
	print $outfhp "<html>\n";
    }

    open($outfhp_map, ">$dir/multi_$num.pep.intronmap.fasaln") || die $!;
    open($outfh_map, ">$dir/multi_$num.cds.intronmap.fasaln") || die $!;
    if( my $aln = $in->next_aln ) {
	$aln->set_displayname_flat(1);
	my $final_seqs = Bio::SimpleAlign->new;
	my $final_pseqs = Bio::SimpleAlign->new;
	my (%add_columnsp,%add_columns,%colmap);
	
	for my $g ( $aln->each_seq ) {
	    my $sp = $g->display_name;
	    if( my $mapname = $num2seqs[$num]->{$sp} ) {
		my $id = "$sp:$mapname";		
		if( ! $idx->offset($id) ) { 
		    if( ! $pepidx->offset($id) ) {
			die("Cannot find $id\n");
		    } else {
			my $header = $pepidx->header($id);
			my ($mapname2) = ($header =~ /transcript:(\S+)/);
			if( $mapname2 ) {
			    my $cds_id = "$sp:$mapname2";
			    if( ! $idx->offset($cds_id) ) {
				die("Cannot find '$id' or '$cds_id'\n");
			    }
			    $id = $cds_id;
			    $mapname = $mapname2;
			} else {
			    die("no transcript mapping in protein header ($header)\n");
			}
		    }
		}
		my $header = $idx->header($id);
		my @d = split(/\s+/,$header,3);
		shift @d; # get rid of ID/name
		my ($chrom,$loc);
		for ( @d ) {
		    if(/(\S+):(\S+)/ ) {
			($chrom,$loc) = ($1,$2);
			last;
		    }
		}
		if( ! $loc ) {
		    # get location from GFF 
		    $loc = $gffidx{$id};
		    if( ! $loc ) {
			warn("cannot find loc for $id in ",$idx->header($id),"\n");
			exit;
		    }
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
			    warn("no column for $runninglen $sp\n");
			} else { 
			    $colmap{$sp}->{$col} = $ct++;
			    push @intron_positions, [$col, $intronphase];
			}
		    }
		    $runninglen += $subloc->length;
		}
		my $gseq = $g->seq();
		my $pseq = $g->translate()->seq();
		
		for my $pos ( @intron_positions ) {
		    next unless defined $pos;
		    $add_columns{$pos->[0]}->{$id} = $pos->[1];
		    my $ppos = int((($pos->[0] -1)/ CODONLEN)) +1;
		    $add_columnsp{$ppos}->{$id} = $pos->[1];
		    
# 		    print $outfhrpt join("\t", 
# 					 $num, 
# 					 $ppos,				      
# 					 $pos->[1],
# 					 $pos->[0],
# 					 $id),"\n";
		    
		    substr($gseq, $pos->[0]+$offset,0,$pos->[1]);
		    substr($pseq, $ppos+$poffset,0,$pos->[1]);
		    $poffset++;
		    $offset++;
		}

		my $ls = Bio::LocatableSeq->new(-verbose => -1,
						-display_id => $id);
		$ls->alphabet('dna');
		$ls->seq($gseq);
		$final_seqs->add_seq($ls);
		my $ps = Bio::LocatableSeq->new(-verbose => -1,
						-display_id => $id);
		#$ps->no_validate(1);
		$ps->alphabet('protein');
		$ps->seq($pseq);
		$final_pseqs->add_seq($ps);
	    }
	}
	$final_seqs->set_displayname_flat(1);
	my %nmset;
	for my $seq ( $final_seqs->each_seq ) {
	    my $gseq = $seq->seq();
	    my $id = $seq->display_id;
	    my ($sp) = split(/_/,$id);
	    $nmset{$id} = $sp;
	    my $offset = 0;
	    for my $col ( sort {$a <=> $b } keys %add_columns ) {
		if( ! defined $add_columns{$col}->{$id} ) {
		    substr($gseq, $col+$offset,0,'-');
		    $add_columns{$col}->{$id} = '-';
		}
		$offset++;
	    }
	    $seq->seq($gseq);
	}
	$final_pseqs->set_displayname_flat(1);
	for my $seq ( $final_pseqs->each_seq ) {
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

	# fasta writing
	my $out = Bio::AlignIO->new(-format => 'fasta',
				 -fh     => $outfh_map);
	$out->write_aln($final_seqs);
	$out = Bio::AlignIO->new(-format => 'fasta',
				 -fh     => $outfhp_map);
	$out->write_aln($final_pseqs);
	$out->close;
	if( $html ) {	    
	    my $str = IO::String->new();	    
	    $out = Bio::AlignIO->new(-format => 'clustalw',
				     -linelength => 100,
				     -fh     => $str);
	    $out->write_aln($final_seqs);	
	    $out->close();
	    my $str2 = IO::String->new(${$str->string_ref});
	    print $outfh "<pre>\n";
	    while(<$str2>) {
		if( /^CLUSTAL|MUSCLE/ ) {
		} elsif( /^(\s*\S+\s+)(\S+)(\s*\d*\s*)$/ox ) {
		    my ($pre,$seq,$post) = ($1,$2,$3);
		    # add 'bold' tags
		    $seq =~ s/(\d+)/<b>$1<\/b>/g;
		    $seq =~ s/(A+)/<font color="$COLOR{'A'}">$1<\/font>/g;
		    $seq =~ s/(C+)/<font color="$COLOR{'C'}">$1<\/font>/g;
		    $seq =~ s/(G+)/<font color="$COLOR{'G'}">$1<\/font>/g;
		    $seq =~ s/(T+)/<font color="$COLOR{'T'}">$1<\/font>/g;
		    $_ = join('',$pre,$seq,$post);
		}
		print $outfh $_;
	    }
	    print $outfh "<\/pre>\n<\/html>\n";
	
	    # pep seqs	
	    $str = IO::String->new();
	    $out = Bio::AlignIO->new(-format => 'clustalw',
				     -linelength => 100,
				     -fh     => $str);
	    $out->write_aln($final_pseqs);
	    $out->close();
	    
	    $str2 = IO::String->new(${$str->string_ref});	
	    print $outfhp "<pre>\n";
	
	    while(<$str2>) {
		if( /^CLUSTAL|MUSCLE/ ) {
		} elsif( /^(\s*\S+\s+)(\S+)(\s*\d*\s*)$/ox ) {
		    my ($pre,$seq,$post) = ($1,$2,$3);
		    # add 'bold' tags
		    $seq =~ s/(\d+)/<b>$1<\/b>/g;
		    $seq =~ s/([DTECMYKRSQFPWNGVILAH])/<font color="$PEPCOLOR{$1}">$1<\/font>/g;
		    $_ = join('',$pre,$seq,$post);
		}
		print $outfhp $_;
	    }
	    print $outfhp "<\/pre>\n<\/html>\n";
	}
    } else { 
	warn("no aln for $num\n");
    }
}
