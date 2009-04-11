#!/usr/bin/perl -w
use strict;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::LocatableSeq;
use List::Util qw(shuffle sum);

use Getopt::Long;

my ($iformat,$oformat,$ext) = ('nexus','phylip','nex');
my $numgenes = 50;
my $min_fraction = 0.25;
my $gapchar = '-';
GetOptions(
	   'mf:f'      => \$min_fraction,
	   'if:s'      => \$iformat,
	   'of:s'      => \$oformat,
	   'ext:s'     => \$ext,
	   'g|gapchar:s'=> \$gapchar,
	   'n|num:s'   => \$numgenes,
	   );

my $dir = shift || '.';

opendir(DIR, $dir ) || die $!;

my ($c,%seqs,@f);
for my $f ( readdir(DIR) ) {
    next unless $f =~ /\.$ext$/;
    push @f, $f;
}
$c = 0;
shuffle(\@f);
my @partitions;
for my $f ( @f ) {
    last if $c >= $numgenes;
    my $in = Bio::AlignIO->new (-format => $iformat,
				-file   => "$dir/$f");
    my ($fname) = split(/\./,$f);
    my $first = 1;
    while( my $aln = $in->next_aln ){
	my $matrix = $aln->gap_col_matrix;
	my %colstoremove;
	my $i = 0;
	for my $cl ( @{$matrix} ) {
	    my $num = scalar keys %$cl;
	    my $sum = sum values %$cl;
	    if( ($sum / $num) > $min_fraction ) {
		$colstoremove{$i}++;
	    }
	    $i++;
	}

	my @cols = sort { $b <=> $a} keys %colstoremove;
	for my $seq ( $aln->each_seq ) {
	    my ($id) = split(/:/,$seq->display_id);
	    my ($sp,$src) = split(/\_/,$id);
	    $src =~ s/rm11-1a/rm11/;
	    my $str = $seq->seq;
	    for my $col ( @cols  ) {
		substr($str,$col,1,''); # remove the column
	    }
	    $str =~ s/X/-/g;
	    if( $first ) { 
		push @partitions, [ $fname, 
				    length($seqs{"$sp\_$src"} || '')+1, 
				    length($seqs{"$sp\_$src"} || '')+
				    length($str)];
		$first = 0;
	    }
	    $seqs{"$sp\_$src"} .= $str;

	}
    }
    $c++;
}

warn("$c genes used\n");
my $newaln = Bio::SimpleAlign->new();
while( my ($id,$seq) = each %seqs ) {
    $newaln->add_seq(Bio::LocatableSeq->new(-id   => $id,
					    -seq  => $seq,
					    -start=>1,
					    -end  => length($seq)));
}
my $out = Bio::AlignIO->new(-show_symbols => 0,
			    -show_endblock => 0,
			    -format => $oformat);
$out->write_aln($newaln);

if( $oformat eq 'nexus' ) {
    print "begin mrbayes;\n";
    
    for my $p ( @partitions ) {
	printf "  charset %s = %d-%d;\n", @$p;
    }
    
    printf " partition all = %d: %s;\n", scalar @partitions, 
    join(", ", map { $_->[0] } @partitions);
    
    print "set partition = all;\n";
    print "prset applyto=(all) ratepr=variable;\n";
    print "end;\n";
}
