
#!/usr/bin/perl -w
use strict;
use Bio::DB::GFF;
use Bio::SeqIO;
use Bio::Perl;
use constant DEBUG     => 0;
use constant CodonSize => 3;
my %dbnames = (
#	       'cneo' => 'JEC21_20040730',
	       #'seroA' => 'seroA',
#	       'ccin' => 'coprinus',
	       #'h99' =>'H99_WI',
	       #'umay' => 'ustilago',
#	       'pchr' => 'pchrysosporium',
	       #'anid' => 'anidulans',
	       #'calb' => 'calbicans19',
#	       'fgra' => 'fusarium',
	       #'fgram'=> 'fusarium',
	       #'ncra' => 'ncrassa',
	       #'mgri' => 'mgrisea',
	       'spom' => 'spombe',
	       #'scer' => 'yeast'
	       #'afum' => 'afumigatus',
	       );
my %agg = ( 
	    'coprinus'  => ['transcript:Genewise' ],
	    'anidulans' => ['transcript:whitehead'  ],
	    'mgrisea'   => ['transcript:whitehead'  ],
	    'ncrassa'   => ['transcript:whitehead'  ],
	    'fusarium'  => ['transcript:whitehead'  ],
	    'calbicans19' => ['transcript:Stanford' ],
	    'yeast'     => ['transcript:SGD'] ,
	    'afumigatus' => ['transcript:Genewise' ],
	    'H99_WI' => ['processed_transcript:Genomewise' ],
	    'JEC21_20040730' => ['transcript:TIGR' ],
	    'seroA'     => ['alignments:SIM4' ],
	    'pchrysosporium'  => ['transcript:DOE'],
	    'spombe'    => ['transcript:gbk'],
	    );

my $defaultagg = ['transcript'];
my $alignments = Bio::DB::GFF::Aggregator->new(-method       => 'alignments',
					       -main_method  => 'match',
					       -sub_parts    => [qw(similarity)]);


for my $species ( values %dbnames ) {
    my $aggs = $agg{$species} || $defaultagg;
    mkdir($species) unless -d $species;
    my $CDS = ($species =~ /afum|cneo|JEC21|H99|coprinus|ustilago|whiterot|calbicans/i) ? 'CDS' : 'exon';
    $CDS = 'similarity' if $species eq 'seroA';
    warn("$species\n");
    my $dbh = Bio::DB::GFF->new(-dsn => "dbi:mysql:database=$species;host=neptune.duhs.duke.edu",
				-user=> 'gbUser', -pass => 'Gbr0wse',
				-aggregators => ['transcript',$alignments,
						 'processed_transcript'],
				) || die($!);
    my $iterator = $dbh->get_seq_stream(-type => $aggs);
    my %contig;
    while( my $f = $iterator->next_seq ) {
	my @exons = map { $_->[0] }
	sort { $a->[1] <=> $b->[1] }
	map { [$_, $_->start*$_->strand] }
	( grep { $_->primary_tag =~ /$CDS/i }
	  $f->sub_SeqFeature );
	my $phase = 0;
	my $runninglength= 0;
	my ($len,$i,$last,$lastf,$intron) = (scalar @exons,1);
	my $transcript;
	my @exon_features;
	for my $exon ( @exons ) {
	    my ($start,$end) = sort { $a <=> $b } ($exon->start, $exon->end);
	    ($start,$end) = ($end,$start) if $f->strand < 0;
	    my $elen = $exon->length;
            next if $elen == 0;
	    my $exonseq = $dbh->get_dna($f->seq_id,
					$start, $end); 

	    warn("phase $phase for ",$i-1," ",$f->group, "\n") 
		if DEBUG && $i > 1;
	    
	    $phase = ($runninglength += $elen) % CodonSize;
	    my $ftype;
	    if( $i == 1
		#&& uc(substr($exonseq,0,3)) eq 'ATG'
		) {
		if( $i == $len ) { # only one exon 
		    $ftype = 'Esngl';
		} else { $ftype = 'Einit' }
	    } elsif( $i == $len ) { $ftype = 'Eterm' }
	    else { $ftype = 'Exon' }


	    push @exon_features, join("\t",  ($ftype, 
					      $start,$end,
					      $f->group));
	    $i++;
	    $transcript .= $exonseq;
	}

	if( $transcript && 
	    &translate_as_string($transcript) =~ /^([^\*]+)\*?$/ ) { 
	    push @{$contig{$f->seq_id}}, @exon_features;
	}
    }
    open my $zff => ">$species/$species.ann" || die $!;
    my $seqio = Bio::SeqIO->new(-file   => ">$species/$species.dna",
				-format => 'fasta');
    for my $c ( keys %contig ) {
	my $seq = $dbh->segment($c);
	$seq->display_id($c);
	$seq->description("");
	$seqio->write_seq(Bio::PrimarySeq->new(-seq => $seq->seq,
					       -id  => $c,
					       ));
	printf $zff ">%s\n",$c;
	for my $v ( @{$contig{$c}} ) {
	    print $zff $v, "\n";
	}
    }
}
