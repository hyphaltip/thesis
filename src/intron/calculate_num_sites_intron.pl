#!/usr/bin/perl -w

# this is to calculate the number of effective intron positions

use strict;
use Getopt::Long;
use Bio::AlignIO;
use Bio::LocatableSeq;
use Bio::TreeIO;
use Math::BigInt;
use File::Spec;
use Bio::AlignIO::fasta;

$Bio::AlignIO::fasta::MATCHPATTERN .= '\d';
$Bio::PrimarySeq::MATCHPATTERN .= '\d';
my %Groups = ( 'H' => [ qw(scas sbay smik skud sklu spar
			   scer_sgd scer_yjm789_SNAP scer_rm11_SNAP
			   cgla agos klac dhan ylip calb ctro cgui clus)],
	       'E' => [ qw(mgri fgra ncra anid cimm afum hcap pans cglo sscl
			   ater bcin snod aory uree fver) ],
	       'A' => [ qw(spom) ],
	       'B' => [ qw(umay pchr ccin
			   cneo_H99 cneo
			   cneo_JEC21 cneo_R265
			   cneo_WM276)],
	       'Z' => [qw (rory) ],
               'O' => [qw(ddis)],
	       'M' => [ qw(hsap cbri cele mmus rnor xtro ggal cfam
			   ptro ggal drer tnig frub agam dmel amel cint)],
	       'P' => [ qw(atha crei)]  );

use constant {
    CODONLEN    => 3,
    MATCH_STATE => 'M',
    INTRON_STATE=> 'I',
    GAP_STATE   => 'G',
};

my $gap = GAP_STATE;
my $windowup = 3;
my $windowdn = 3;
my $debug = 0;
my $ct = 0;
my $other = 0;
my $fullname = 0;
my $phylofile;
my $pepalndir;
my $min_pid = 40;
GetOptions(
	   'p|pepaln:s'            => \$pepalndir,
	   'f|full|fullname'       => \$fullname,
	   't|tree:s'              => \$phylofile,
	   'wd|windowdn:i'         => \$windowdn,
	   'wu|windowup:i'         => \$windowup,
	   'v|verbose'             => \$debug,
#	   'min_pid:i'             => \$min_pid,
	   );

die("need a treefile with -tree option") unless $phylofile && -f $phylofile;
my $tformat ='newick';
open(T, $phylofile) || die $!;
if(<T> =~ /\#NEXUS/) {
    $tformat = 'nexus';
}
close(T);
my $tree = Bio::TreeIO->new(-format => $tformat,
			    -file   => $phylofile)->next_tree;
my %phylogeny;
my (%reverse_phy,@tips);
open(MAPPING, ">sp_id_mapping.numsites");
for my $tip ( map { $_->[0] }   # schwartzian transformation
	      sort {$a->[1] cmp $b->[1] } 
	      map { [$_,$_->id] } $tree->get_leaf_nodes ) {
    my $n = Math::BigInt->new(2)->bpow($ct++);
    $phylogeny{$tip->id} = $n;
    print MAPPING join("\t",$tip->id,$n),"\n";
    $reverse_phy{$n} = $tip->id;
    $tip->id($n);
    push @tips, $tip;
}
close(MAPPING);
$pepalndir ||= 'ortholog_sequence_sets-aln';

my ($debugaln,$debugfh);
if( $debug ) {
    open($debugfh => ">aln.debug2") || die $!;
    $debugaln = Bio::AlignIO->new(-format => 'clustalw',
				  -fh     => $debugfh);
}
warn("processing dir $pepalndir\n");
opendir(SETS, $pepalndir) || die("$pepalndir: $!");
for my $file ( map { $_->[0] } 
	       sort {$a->[1] <=> $b->[1] } 
	       map { [$_, /multi_(\d+)\.pep.fasaln$/] } grep {/multi_\d+\.pep.fasaln$/ } readdir(SETS) ) {
    next unless $file =~ /multi_(\d+)\.pep.fasaln$/;
    my $num = $1;
    my $aln = Bio::AlignIO->new(-format => 'fasta', 
				-file => File::Spec->catfile($pepalndir,
							     $file))->next_aln;
    $aln->verbose(-1);
    my $len = $aln->length;
    my @seqs;
    my @names;
    my $spct;    
    for my $g ( $aln->each_seq ) {
	push @seqs, $g->seq;
	if( $fullname ) {
	    push @names, (split(/:/,$g->display_id))[0]
	} else {
	    push @names, (split(/_/,$g->display_id))[0];
	}
	$spct++;
    }
    my $match_line = $aln->match_line;
    my $state_str;
    my @all_cols;
    print $debugfh "len is $len\n" if $debug;
    for( my $i = 0; $i < $len; $i++) {
	my $in_gap = 0;
	my @col;
	my $snum = 0;
	my $phase;
	for my $s ( @seqs ) {	    
	    my $char = substr($s,$i,1);
	    if( $char eq '-' || $char eq '.' || $char eq '~' ) {
		$in_gap++;
	    }
	    push @col, $char;
	    $snum++;
	}
	if ( $in_gap ) {
	    $state_str .= GAP_STATE;
	} else {
	    $state_str .= MATCH_STATE;
	}
	push @all_cols, [$i+1,$in_gap,\@col];
    }
    
    for my $column ( @all_cols ) {
	my ($coln,$state,$col_val) = @$column;
	
	my $l = $coln;
	my $r = $coln;
	# warn("$l $r phase=$phase coln=$coln intron_ct=$intron_ct\n");
	
	$l -= $windowup;
	$r += $windowdn;

	if( $l <= 0 ) { $l = 1 }
	if( $r > $len ) { $r = $len }
	print $debugfh ("l is $l, r is $r\n") if $debug;
	my $slice = $aln->slice($l,$r,1);
	my $skip = 'keep';
	if( $slice->gap_line =~ /\-/ ) {
	    $skip = 'has_gap';
	} 
	if( $debug ) {
	    $debugaln->write_aln($slice);
	    printf $debugfh "site is @$col_val \n";
	    printf $debugfh "num=$num col is $coln pid is %.2f  pairwise avg pid is %.2f for $coln status=$skip\n",
	    $slice->average_percentage_identity,
	    &paired_percent_identity($slice);
#
#	    $debugaln->write_aln($aln->slice($coln - 
#					     $intron_ct - (3*($windowup+1)),
#					     $coln - $intron_ct+(3*$windowdn)));
	    printf $debugfh "----\n\n";
	}
	my $ppid = &paired_percent_identity($slice);
	$ppid = $slice->average_percentage_identity unless defined $ppid;
	if( $skip eq 'keep' && $ppid < $min_pid) {
	    $skip = "low_id";
	}
	print join("\t",
		   sprintf("%04d:",$num),
		   $coln,
		   $len,
		   $skip,
		   sprintf("%.2f",$slice->average_percentage_identity),
		   sprintf("%.2f",$ppid),
		   ), "\n";
#	$intron_ct++;
    }
}

sub paired_percent_identity {
    my $aln = shift;
    my %groups;
    my @seqs = $aln->each_seq();
    my $len = $aln->length();
    my @pid;
    my @grps = keys %Groups;
    for( my $i =0; $i < scalar @grps; $i++ ) {
	my @grp1_seqs;
	for my $id ( @{$Groups{$grps[$i]}} ) {
	    for my $s ( @seqs ) {
		if($s->display_id =~ /^$id/ ) {
		    push @grp1_seqs, $s;
		}
	    }
	}
	for( my $j = $i+1; $j < scalar @grps; $j++ ) {
	    my @grp2_seqs;
	    for my $id ( @{$Groups{$grps[$j]}} ) {
		for my $s ( @seqs ) {
		    if($s->display_id =~ /^$id/ ) {
			push @grp2_seqs, $s;
		    }
		}			
	    }
#	    next unless( @grp1_seqs && @grp2_seqs ); 
#	    warn( join(',', map { ref($_) ? $_->display_id : $_} @grp1_seqs, "\n"));
#	    warn( join(',', map { ref($_) ? $_->display_id : $_ } @grp2_seqs, "\n"));
	    my $avg_pid = &average_pairwise(\@grp1_seqs, \@grp2_seqs);
	    push @pid, $avg_pid if defined $avg_pid;
	}
    }
    my $sum;
    for ( @pid ) { $sum += $_ }
    if( ! @pid ) {
	return undef;
    } else { 
	return ($sum / (scalar @pid));
    }
}

sub average_pairwise {
    my ($grp1,$grp2,$len) = @_;
    my @pid;
    for my $s1 ( @$grp1 ) {
	my $s1seq = $s1->seq;
	for my $s2 ( @$grp2 ) {
	    my $s2seq = $s2->seq;
	    my $id_bases = 0;
	    for( my $i = 0; $i < $s1->length; $i++ ) {
		my ($ch1,$ch2) = (substr($s1seq,$i,1),
				  substr($s2seq,$i,1));
		next if( $ch1 eq '-' || $ch2 eq '-');
		if( $ch1 eq $ch2 ) {
		    $id_bases++;
		}
	    }
#	    warn($id_bases, " for \n", $s1->display_id, ":\t",
#		 $s1seq, "\n",$s2->display_id, ":\t", $s2seq," $id_bases\n");
	    push @pid, 100 * ($id_bases / $s1->length);
	}
    }
    my $sum;
    for ( @pid ) { $sum += $_; }
    #warn(sprintf("g1=%s g2=%s pid=%.2f\n", 
    #	 join(",",map { $_->display_id} @$grp1),
    #	 join(",",map { $_->display_id} @$grp2),
    #	 @pid ? $sum / scalar @pid : 0));
    
    if( ! @pid ) {
	return undef;
    } else { 
	return $sum / scalar @pid;
    }
}

# 	my ($l,$r) = ( $coln - $windowup - $phase,
# 		       $coln + $windowdown + $phase);
# 	# boundary cases
# 	if( $l < 0 ) {$l = 0 }
# 	if( $r >= $len ) { $r = $len-1 }
# 	my $state_window = substr($state_str, $l, abs($r-$l)+1);
# 	my $skip = 'keep';
# 	if( $state_window =~ /$gap/ ) {
# 	    $skip = 'skip_gap';
# 	} else { 
# 	    my $slice = $aln->slice($l+1,$r+1);
# 	    my $newslice = Bio::SimpleAlign->new();
# 	    for my $seq ( $slice->each_seq ) {
# 		$newslice->add_seq($seq->translate);
# 	    }
# 	    if( $newslice->average_percentage_identity < $min_pid ) {
# 		$skip = 'skip_lowpid';
# 		if( $debug ) {
# 		    printf $debugfh "pid is %.2f for $coln\n",
# 		    $newslice->average_percentage_identity;
# 		    $debugaln->write_aln($newslice);
# 		}
# 	    }
# 	}
