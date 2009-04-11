#!/usr/bin/perl -w
use strict;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::LocatableSeq;
use List::Util qw(shuffle sum);

use Getopt::Long;

my ($iformat,$oformat,$ext) = ('nexus','phylip');
my $min_fraction = 0.25;
my $gapchar = '-';
my ($file,$output);
GetOptions(
	   'mf:f'      => \$min_fraction,
	   'if:s'      => \$iformat,
	   'of:s'      => \$oformat,
	   'i|in:s'    => \$file,
	   'o|out:s'   => \$output,
	   'g|gapchar:s'=> \$gapchar,
	   );

$file ||= shift @ARGV || die;

my $in = Bio::AlignIO->new (-format => $iformat,
			    -file   => $file);
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

while( my ($id,$seq) = each %seqs ) {
    $newaln->add_seq(Bio::LocatableSeq->new(-id   => $id,
					    -seq  => $seq,
					    -start=>1,
					    -end  => length($seq)));
}
my $out = Bio::AlignIO->new(-show_symbols  => 0,
			    -show_endblock => 0,
			    -format        => $oformat);
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


sub strnum2range {
    return unless @_;
    local $_ = shift;
    s/(?<!\d)(\d+)(?:,((??{$++1})))+(?!\d)/$1-$+/g;
    return $_;
}

sub list2ranges {
    my @list = sort {$a<=>$b} @_;
    my $current = $list[0];
    my @ranges = ($current) ;
    for (my $i=1;$i <= @list; $i++){
        if ($list[$i] && $list[$i] - $current == 1 || $list[$i] && $list[$i] - $current < 1  && 
            substr($list[$i],0,-1) eq substr($current,0,-1) ){
            $current = $list[$i];
        }
        else{
              $ranges[-1] .= "-$current" if $ranges[-1] != $current;
              $list[$i] && push @ranges, $current = "".$list[$i];
        }
    }

    return wantarray ? @ranges : join(",",@ranges);
}

sub num2range {
    return unless @_;
    local $_ = join ',' => sort { $a <=> $b } @_;
    s/(?<!\d)(\d+)(?:,((??{$++1})))+(?!\d)/$1-$+/g;
    return $_;
}
