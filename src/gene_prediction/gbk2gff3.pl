#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Bio::Tools::GFF;

my $pref = shift || die "need a prefix";

my $gffout = Bio::Tools::GFF->new(-gff_version => 3, 
				  -file        => ">$pref.ncbi.gff3");
my $seqout = Bio::SeqIO->new(-format => 'fasta',
			     -file   => ">$pref.ncbi.fa");

my $pepout = Bio::SeqIO->new(-format => 'fasta',
			     -file   => ">$pref-gbk_pep.fa");

my $cdsout = Bio::SeqIO->new(-format => 'fasta',
			     -file   => ">$pref-gbk_cds.fa");

my $in= Bio::SeqIO->new(-format => 'genbank',
			-fh => \*ARGV);

my %extras_to_grab = ( 'protein_id' => 'genpept_id',
		       'db_xref'    => 'db_xref',
		       'transl_table' => 'translation_table',
		       'codon_start' => 'codon_start',
		       'evidence'   => 'evidence',
		       );
my %mappings = ( 'tRNA' => 'tRNA',
		 'CDS'  => 'gene',
		 'rRNA' => 'rRNA',
		 'misc_feature' => 'misc_feature',
		 'mRNA' => 'mRNA',
                 'LTR' => 'LTR',
                 'conflict' => 'conflict',
		 'unsure' => 'low_quality');
my $source = 'gbk';
my %count;
while( my $seq = $in->next_seq ) {
    my %seen; # insure uniqueness per chromosome
    my @features = $seq->get_SeqFeatures;
    my ($src) = grep { $_->primary_tag eq 'source' } @features;
    my $seqid = $seq->display_id;
    for my $f ( @features ) {
	my $note = '';
	if( $f->has_tag('note') ) {
	    ($note) = $f->get_tag_values('note');
	}
	my $ptag = $f->primary_tag;
	if($ptag eq 'source') {
	    if( $f->length == $seq->length ) {
		$gffout->write_feature(Bio::SeqFeature::Generic->new
				       ( -start => $f->start,
					 -end   => $f->end,
					 -strand=> 1,
					 -seq_id=> $seqid,
					 -primary_tag => 'Component',
					 -source_tag  => 'chromosome',
					 -tag   => { 
					     'ID' => "Sequence:".$seqid,
					 }));
	    } elsif( $f->has_tag('clone') ) { 
		# BAC clones
		my ($clone) = $f->get_tag_values('clone');
		if( ! ($clone =~ s/BAC\s+clone\s+//) ) {
		    warn("clone $clone did not match\n");
		}
		$gffout->write_feature(Bio::SeqFeature::Generic->new
				       ( -start => $f->start,
					 -end   => $f->end,
					 -strand=> 1,
					 -seq_id=> $seqid,
					 -primary_tag => 'assembly',
					 -source_tag  => 'BAC',
					 -tag   => { 
					     'ID' => "BAC:".$clone,
					 }));
	    }
	} elsif( $ptag eq 'CDS' || $ptag eq 'tRNA' || $ptag eq 'mRNA' ) {
	    my $gn;
	    if( $ptag eq 'tRNA' ) {
		if( $f->has_tag('product') ) {
		    ($gn) = $f->get_tag_values('product');
		} else { 
		    warn("no product tag for ", $f->gff_string, "\n");
		    ($gn) = 'unk-tRNA';
		}
		# insure that tRNAs gene names are unique
		$gn = "$gn.".$seen{$gn}++;
	    } else { 
		if( $f->has_tag('gene') ) {
		    ($gn) = $f->get_tag_values('gene');
		} elsif( $f->has_tag('locus_tag') ) { 
                    ($gn) = $f->get_tag_values('locus_tag');
                } elsif( $note =~ /mRNA from (\S+)/) {
		    $gn = $1;
		} else {
		    warn("no locus for ", $f->gff_string, "\n");
		    next;
		}
	    }

	    if( $note =~ /^(\S+)\s*\,/ ) {
		my ($nm) = $1;
		if( $gn ne $nm && $nm !~ /anticodon/ && $nm !~ /^len/ ) {
		    $note = "$gn $note";
		    $gn = $nm;
		}
	    }
	  if( $seen{$gn}++ > 1) {
		$gn .= "-".$seen{$gn};
	    }
	    my $mapptag = $mappings{$ptag};
	    my $genef = Bio::SeqFeature::Generic->new
				   ( -start => $f->start,
				     -end   => $f->end,
				     -seq_id => $seqid,
				     -strand=> $f->strand,
				     -primary_tag => $mapptag,
				     -source_tag  => $source,
				     -tag   => { 
					 'ID'   => "Gene:$gn",
				     });
	    while( my ($tag,$newtag) = each %extras_to_grab) {
		next unless $f->has_tag($tag) ;
		$genef->add_tag_value($newtag, $f->get_tag_values($tag));
	    }
	    $genef->add_tag_value('Note', $note) if defined $note;
	    my ($tbl) = $f->has_tag('transl_table') ? 
		$f->get_tag_values('transl_table') : 1;
	    my ($frame) = $f->has_tag('codon_start') ? $f->get_tag_values('codon_start') : 1;
	    $frame--;
	    
	    if( $ptag ne 'mRNA' ) {
		my $exon = ( $ptag eq 'tRNA' ) ? 'tRNA_exon' : 'exon';
	        my $first = 1;
		for my $loc ( $f->location->each_Location ) {
		    my ($s,$e) = ($loc->start,$loc->end);
	 	    if( $first && $frame ) {
			if( $loc->strand < 0 ) { $e -= $frame;
			} else { $s += $frame }
		    }
		    $first = 0;
		    my $exon = Bio::SeqFeature::Generic->new
			( -start => $s,
			  -end   => $e,
			  -strand=> $loc->strand,
			  -primary_tag => $exon, # for tRNA what do we do
			  -source_tag  => $source,
			  -seq_id => $seqid,
			  -tag => {
			      'Parent' => "Gene:$gn",
			  });		    
		    $gffout->write_feature($exon);
		}
	    }
	    if( $ptag eq 'CDS' ) {
		my $cds = $f->spliced_seq(undef,1);
		$cds->display_id($gn);
		$cds->description("$seqid $note");
	        if( $frame ) {
		    if( $genef->strand < 0 ) {
			$genef->end($genef->end - $frame);
		    } else { 
			$genef->start( $genef->start + $frame);
		    }
		    $cds->seq( substr($cds->seq(),$frame) );
		}
		$cdsout->write_seq($cds);
		$pepout->write_seq($cds->translate(undef,undef,undef,$tbl));
	    }
	    $gffout->write_feature($genef);
	} elsif( $ptag eq 'unsure' || $ptag eq 'misc_feature' || 
		 $ptag eq 'conflict' || 
	         $ptag eq 'LTR' || $ptag eq 'promoter' ||
		 $ptag eq 'repeat_region' ) {
	    my $mapptag = $mappings{$ptag} || $ptag;
	    if( $note =~ /centromere/ ) {
		my $unsure = Bio::SeqFeature::Generic->new
		    ( -start => $f->start,
		      -end   => $f->end,
		      -strand=> $f->strand,
		      -seq_id=> $seqid,
		      -primary_tag => 'CEN',
		      -source_tag => $source,
		      -tag => { 
			  'ID' => "CEN".$seqid,
		  });
		while( my ($tag,$newtag) = each %extras_to_grab) {
		    next unless $f->has_tag($tag);
		    $unsure->add_tag_value($newtag, $f->get_tag_values($tag));
		}
		$unsure->add_tag_value('Note', $note) if defined $note;
		$gffout->write_feature($unsure);
	    } else { 
		my $unsure = Bio::SeqFeature::Generic->new
		    ( -start => $f->start,
		      -end   => $f->end,
		      -strand=> $f->strand,
		      -seq_id=> $seqid,
		      -primary_tag => $mapptag,
		      -source_tag => $source,
		      -tag => { 
			  'ID' => "$mapptag.".$count{$mapptag}++,
		  });
		while( my ($tag,$newtag) = each %extras_to_grab) {
		    next unless $f->has_tag($tag);
		    $unsure->add_tag_value($newtag, $f->get_tag_values($tag));
		}
		$unsure->add_tag_value('Note', $note) if defined $note;
		$gffout->write_feature($unsure);
	    }
	} elsif( $ptag =~ /RNA/ ) {
	    my $mapptag = $mappings{$ptag} || $ptag;
            my $product;
            if( ! $f->has_tag('product') ) {
             # warn("No product tag for ", $f->gff_string, "\n");
             $product = "$mapptag.".$count{$mapptag}++;
            } else {
	      ($product) = $f->get_tag_values('product');
            }
	    my $rna = Bio::SeqFeature::Generic->new
		( -start => $f->start,
		  -end   => $f->end,
		  -strand=> $f->strand,
		  -seq_id=> $seqid,
		  -primary_tag => $mapptag,
		  -source_tag => $source,
		  -tag => { 
		      'ID' => $product,
		  });
	    while( my ($tag,$newtag) = each %extras_to_grab) {
		next unless $f->has_tag($tag) ;
		$rna->add_tag_value($newtag, $f->get_tag_values($tag));
	    }
	    $rna->add_tag_value('Note', $note) if defined $note;
	    $gffout->write_feature($rna);
	} elsif( $ptag eq 'gene' || $ptag eq 'intron' || $ptag eq 'exon') { 
	} else { 
	    warn("ptag $ptag\n");
	}
    }
    $seq->display_id($seqid);
    $seqout->write_seq($seq);
}
