#!/usr/bin/perl -w
use strict;
use Bio::AlignIO;
use Env;
use File::Spec;
use File::Copy;
use Cwd;

my $min_percent = 0.75;

my $indir = shift || 'ortholog_sequence_sets-aln';
$indir= File::Spec->rel2abs($indir);
my (undef,undef,$cwd) = File::Spec->splitpath(cwd);
my $outdir = "/scratch/jes12/$cwd\_masked";
mkdir $outdir unless -d $outdir;
my $statsfh;
open($statsfh => ">alignment.stats.tab") || die $!;

print $statsfh join("\t", qw(ALN LENGTH_BEFORE LENGTH_AFTER 
			    5PRIME_REMOVED 3PRIME_REMOVED)), "\n";
 
opendir(DIR, $indir) || die $!;
my $i =0;
for my $file ( readdir(DIR) ) {
    next unless $file =~ /(multi_(\d+)\.(cds|pep)).(fasaln|nex)$/;
    my $fmt = $4 eq 'nex' ? 'nexus' : 'fasta';
    my $stem = $1;
    my ($num,$type) = ( $2,$3);
    mkdir(File::Spec->catfile($outdir,$type)) unless -d File::Spec->catfile($outdir,$type);
    my $alnstats_fh;
    my $alnin = Bio::AlignIO->new(-format => $fmt,
				  -file   => File::Spec->catfile($indir,
								 $file));
    if( my $aln = $alnin->next_aln ) {
	$aln->missing_char('?');
	my $out = Bio::AlignIO->new(-format => 'nexus',
				    -file   => ">".
				    File::Spec->catfile($outdir,$type,
							"$stem.nex"));
	
	open($alnstats_fh => ">".File::Spec->catfile
	     ($outdir,"multi_$num.$type.stats.tab")) || die $!;
	print $alnstats_fh join("\t", qw(ID UNALIGNED_LEN ALIGNED_LEN 
					 5PRIME_GAPS 3PRIME_GAPS)), "\n";

	# now walk down the columns, first left to right
	my $gap_line = $aln->gap_line;
	my ($p3_gaps, $p5_gaps) = (0,0);
	if( $gap_line =~ /^(\-+)/ ) {
	    $p5_gaps = length($1);
	}
	if( $gap_line =~ /(\-+)$/ ) {
	    $p3_gaps = length($1);
	}
	my $alen = $aln->length;
	print $statsfh join("\t", $num, $alen, 
			    $alen - ($p5_gaps + $p3_gaps),
			    $p5_gaps, $p3_gaps), "\n";

	my @cols;
	for my $seq ( $aln->each_seq ) {
	    my $str = $seq->seq();
	    for my $col ( 0..$alen ) {
		if( substr($str,$col,1) eq '-' ) {
		    $cols[$col]++;
		}
	    }
	    my ($trailing,$leading) = (0,0); 
	    if( $str =~ /^(\-+)/ ) { 
		$leading = length($1);
		my $v = $1;
		my $n = '?'x$leading;
		$str =~ s/^$v/$n/;
	    }

	    if( $str =~ /(\-+)$/ ) {
		$trailing = length($1);
		my $v = $1;
		my $n = '?'x$trailing;
		$str =~ s/$v$/$n/;
	    }
	    $seq->seq($str);

	    print $alnstats_fh join("\t",
				    $seq->display_id,
				    $seq->end - $seq->start + 1, 
				    $seq->length,
				    $leading,
				    $trailing), "\n";
	}
	my $num_seqs = $aln->no_sequences;
	my @gaps;
	for my $col (0..$alen ) {
	    $cols[$col] ||= 0;
	    my $keep = (1 - ($cols[$col] /  $num_seqs)) > $min_percent ? 1 : 0;
	    if( $keep ) {
		if( $col == $alen ) {
		    push @gaps, $col;
		} else { 
		    push @gaps, $col+1;
		}
	    }
	}
	my @ranges= &list2ranges(@gaps);

	$out->write_aln($aln);
	my $fh = $out->_fh;
	print $fh "BEGIN sets;\n";
	if( $p5_gaps && $p5_gaps >= $alen ) {
	    printf $fh "\t\tcharset fiveprimegap  = %d-%d;\n",1,$alen;
	} elsif ( $p5_gaps ) {
	    printf $fh "\t\tcharset fiveprimegap  = %d-%d;\n",1,$p5_gaps;
	} else {
	    print $fh "\t\tcharset fiveprimegap  =; \n";
	}
	if( !$p3_gaps || $p3_gaps >= $alen ) {
	    printf $fh "\t\tcharset threeprimegap = ;\n";
	} elsif ( $p3_gaps ) {
	    printf $fh "\t\tcharset threeprimegap = %d-%d;\n", 
	    1+$alen - $p3_gaps, $alen;
	}
	
	printf $fh "\t\tcharset nogap = %s;\n",join(" ", @ranges);
	if( $p5_gaps == $p3_gaps && $p5_gaps >= $alen ) {
	    printf $fh "\t\tcharset alncore =;\nend;\n";
	} elsif ($p5_gaps == 0 && $p3_gaps == 0 ) {
	    printf $fh "\t\tcharset alncore =%d-%d;\nend;\n", 1, $alen;
	} elsif ($p3_gaps == 0) {
	    printf $fh "\t\tcharset alncore =%d-%d;\nend;\n", 1+$p5_gaps, 
	    $alen;
	} elsif ($p5_gaps == 0) {
	   printf $fh "\t\tcharset alncore =%d-%d;\nend;\n", 1, 
	   $alen - $p3_gaps; 
        } else {
	    printf $fh "\t\tcharset alncore = %d-%d;\nend;\n",
	    1+$p5_gaps,$alen - $p3_gaps;
	}

	
	$out->close();
    }
    $i++;
    warn("$i processed\n") unless ($i % 100);
}

sub num2range {
    return unless @_;
    local $_ = join ',' => sort { $a <=> $b } @_;
    s/(?<!\d)(\d+)(?:,((??{$++1})))+(?!\d)/$1-$+/g;
    return $_;
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

