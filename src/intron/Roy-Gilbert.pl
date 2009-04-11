#!/usr/bin/perl -w
# This is the script which does everything!

use strict;
use Bio::TreeIO;
use Getopt::Long;
use POSIX;
my $ct = 0;
my $debug;
my $min_pid = 0;
my $filter =0;
GetOptions(
	   'mp|min_pid:s' => \$min_pid,
	   'v|verbose'    => \$debug,
	   'f|filter'     => \$filter,
	   );

my @states = ("", qw(n1 n2 n12 n3 n13 n23 n123
		     n4 n14 n24 n124 n34 n134 n234 n1234));

my $phylofile = shift || die "Need a treefile\n";
my $matrix    = shift || die "Need a matrix\n";

my ($nm) = ($matrix =~ /(\S+)\.\S+$/);
my $outtre = shift || "$nm.exloss.tre";

my $tformat ='newick';
open(T, $phylofile) || die $!;
if(<T> =~ /\#NEXUS/) {
    $tformat = 'nexus';
}
close(T);

my $tree = Bio::TreeIO->new(-format => $tformat,
			    -file   => $phylofile)->next_tree;

my (%phylogeny,%reverse_phy,@tips);
for my $tip ( map { $_->[0] }   # schwartzian transformation
	      sort {$a->[1] cmp $b->[1] } 
	      map { [$_,$_->id] } $tree->get_leaf_nodes ) {
    my $n = 2**$ct++;
    $phylogeny{$tip->id} = $n;
    # warn("$n ,",$tip->id,"\n");
    $reverse_phy{$n} = $tip->id;
    $tip->id($n);
    push @tips, $tip;
}

my (%seen,@node,%exnodes,%lookup2);
for my $node ( grep { ! $_->is_Leaf } $tree->get_nodes ) {
    my @tips = grep { $_->is_Leaf } $node->get_all_Descendents;
    my $sum = 0;
    for ( @tips ) { $sum += $_->id }
    my $childstr = $sum. " ".join(",", map { $_->id } @tips)."\n";
    $lookup2{$sum} = $node->id;
    $node->id($sum) if $sum;
}

for my $node ( $tree->get_leaf_nodes ) {
    next unless ( $node->ancestor );
    #next if $seen{$node->id};
    my $left = $node->id;
    my $right = $node->ancestor->id - $left;
    my @xtips = grep { $_->internal_id != $node->internal_id } grep { $_->is_Leaf } $node->ancestor->get_all_Descendents;
    my $childstr = 'external='.$reverse_phy{$left}  .'; '.join(",", map { $reverse_phy{$_->id} } @xtips)."\n";
    push @node, [$left,$right,$childstr];
    $exnodes{$left}++;
    $exnodes{$right}++;
    $seen{$left}++; $seen{$right}++;
}

@node = sort {$a->[0] <=> $b->[0] } @node;

my $tree2 = Bio::TreeIO->new(-format => $tformat,
			     -file   => $phylofile)->next_tree;

my @alltips = keys %phylogeny;
my @vals;
NODE: 
    for my $node ( grep { ! $_->is_Leaf } $tree2->get_nodes ) {
    next unless( $node->ancestor && $node->ancestor->ancestor );
    my @tips = grep { $_->is_Leaf } $node->get_all_Descendents;
    my $sum = 0;
    my $int_nodeid;
    for ( @tips ) { $sum += $phylogeny{$_->id} }
    my $y = $node->ancestor;
    my (@g1,@g2,@g3,@g4);

    # G3
    for my $child ( $y->each_Descendent ) {
	next if $child->internal_id == $node->internal_id;
	push @g3, grep { $_->is_Leaf } ($child,$child->get_all_Descendents);
    }
    # G1 & G2
    if ( $node->each_Descendent > 2 ) { 
	my $childstr = $sum. " ".join(",", map { $_->id } @tips)."\n";
	warn("Not a bifurcating tree at node $childstr\n");
    } else {
	my ($l,$r) = $node->each_Descendent;
	@g2 = grep { $_->is_Leaf } ($l,$l->get_all_Descendents);
	@g1 = grep { $_->is_Leaf } ($r,$r->get_all_Descendents);
    }
    # do G4
    for my $t ( @alltips ) {
	if( ! grep { /$t/i } map { $_->id } (@g1,@g2,@g3) ) {
	    push @g4, $t;
	}
    }
#    print "g1 is ", join(',',map { ref($_) ? $_->id : $_ } @g1), " $childstr";
#    print "g2 is ", join(',',map { ref($_) ? $_->id : $_ } @g2), " $childstr";
#    print "g3 is ", join(',',map { ref($_) ? $_->id : $_ } @g3), " $childstr";
#    print "g4 is ", join(',',map { ref($_) ? $_->id : $_ } @g4), " $childstr";
#    print "\n";
    my ($g1,$g2,$g3,$g4);
    for ( @g1 ) {
	$g1 += $phylogeny{$_->id};
    }
    for ( @g2 ) {
	$g2 += $phylogeny{$_->id};
    }
    for ( @g3 ) {
	$g3 += $phylogeny{$_->id};
    }
    for ( @g4 ) {
	$g4 += $phylogeny{$_};
    }
    #next if @g1 <= 1 || @g2 <= 1;
    my $childstr = $lookup2{$g1+$g2+$g3} . " sum=$sum g1=". join(",", map { $_->id } @g1).' g2='.
	join(",", map { $_->id } @g2)." g3=".
	join(",", (map {$_->id} @g3)) . " g4=".join(",",@g4);
    
    push @vals, [$sum,$g1,$g2,$g3,$g4,$childstr];
#    $childstr = $sum. " ".join(",", map { $_->id } @g2).';'.
#	join(",", map { $_->id } @g1)."-vs-".
#	join(",", (map {$_->id} @g3));
#
#    push @vals, [$sum,$g2,$g1,$g3,$g4,$childstr];
}

@vals = sort {$a->[0] <=> $b->[0] } @vals;


open(IN,$matrix) || die "$matrix: $!";

my (%tally);
while(<IN>) {
    next if /^\s*$/;
    my ($id,$type,@line) = split;
    my ($status,$avgpid,$paired_pid) = ($line[-3],$line[-2],$line[-1]);
    next if $paired_pid < $min_pid;
    next if $filter && $status ne 'keep';
    warn($_);;
    $tally{$type}++
}
close(IN);
my %counts;
my $max = $ct;
foreach my $pow (0..($max-1)) {
    my $n = 2**$pow;
    my $x = 0;
    for my $code (keys %tally) {
	if ($code & $n) {
	    $x += $tally{$code};
	}
    }
    print join("\t",$n,$reverse_phy{$n},$x)," ($pow)\n";
    $counts{$n} = $x;
    $reverse_phy{$n} .= "_$x";
}
my $i =0;
%seen = ();
my %rename;
print "external branch stuff + ancestral state calculation (MLE)\n";
# external node calculations + ancestral state calculation
for my $n ( @node ) {
    my ($g1,$g2,$childstr) = @$n;
    my $g3 = (2**$ct -1) - $g1 - $g2;
#    warn("$g1 $g2 $g3 ",(2**$ct),"\n");
    my @tally2;
    for my $code ( keys %tally ) {
	my $t =  ( (($g1&$code)>0) *1 + 
		   (($g2&$code)>0) *2 +
		   (($g3&$code)>0) *4);
	$tally2[$t] += $tally{$code};
    }
    my $sum = $g1+$g2;
    print "$lookup2{$sum} sum=$sum $g1:$g2 $childstr\n";
    
    for my $j ( 1..7 ) {
	$tally2[$j] ||= 0;
	print join("\t", '',$states[$j],$tally2[$j]),"\n";
    }
    $tally2[7] || next; #skip when no introns found in common in all groups

    my $MLE = $tally2[7]*($tally2[3]+$tally2[7])*($tally2[5]+$tally2[7])*($tally2[6]+$tally2[7])/($tally2[7]**3);
    
    my ($node) = $tree->find_node(-id => $g1);
    if( ! $node ) {
	warn("no node for $g1\n");
    }
    if( ! $seen{$node->ancestor->internal_id}++ ) {
	$rename{$node->ancestor->id} = sprintf("%s_%d",$node->ancestor->id,
					       POSIX::ceil($MLE));
    }
    if( $exnodes{$node->id} ) {
	my $MLL = $MLE*($tally2[6]/($tally2[6]+$tally2[7]));
	my $MLG = $tally2[1]+$tally2[3]+$tally2[5]+$tally2[7]-$MLE+$MLL;
	$rename{$node->id} = sprintf("%s_%d/%d",
				     $node->id,
				     POSIX::ceil($MLL),
				     POSIX::ceil($MLG));
	printf "MLE(N) = %d\nMLE(L) = %d\nMLE(G) = %d\n",
	POSIX::ceil($MLE),
	POSIX::ceil($MLL),
	POSIX::ceil($MLG);
#	printf "MLL/site/blen=%.2f MLG/site/blen=%.2f count=%d\nMLL/blen=%.2f MLG/blen=%.2f (blen=%s,id=%s)\n",
#	($MLE / $node->ancestor->branch_length / $tally2[7]),
#	($MLG / $node->ancestor->branch_length / $tally2[7]),
#	$tally2[7],
#	($MLE / $node->ancestor->branch_length),
#	($MLG / $node->ancestor->branch_length),
#	$node->ancestor->branch_length,
#	$node->ancestor->id;
	print "\n";
	
    }
    $i++;
}
print "internode/branch calculations\n";
# internode calculations
for my $i ( 0..$#vals ) {
    my @tally2;
    my ($internal_node_id,$g1,$g2,$g3,$g4,$childstr) = @{$vals[$i]};
    for my $code ( keys %tally ) {
	my $t = ( 
		  (( $g1 & $code) > 0)*1 + 
		  (( $g2 & $code) > 0)*2 + 
		  (( $g3 & $code) > 0)*4 + 
		  (( $g4 & $code) > 0)*8 
		  );
	$tally2[$t] += $tally{$code};
    }
    print "$childstr $g4\n";
    for my $j ( 1,2,4,8, 3,5,6,9,10,12, 7,11,13,14,15) {
	$tally2[$j] ||=0;
	$states[$j] ||=0;
	print join("\t",'',$states[$j],$tally2[$j]),"\n";;
    }
    next unless( $tally2[6] + $tally2[7] + $tally2[10] + $tally2[11] + 
		 $tally2[14] + $tally2[15]) &&
		 ( $tally2[5] + $tally2[7] + $tally2[9] + $tally2[11] + 
		   $tally2[13] + $tally2[15]);
    next unless ( $tally2[12] + $tally2[13] + $tally2[14] + $tally2[15] );
    my $o1 = ($tally2[7] + $tally2[11] + $tally2[15]) / 
	( $tally2[6] + $tally2[7] + $tally2[10] + $tally2[11] + 
	  $tally2[14] + $tally2[15] );
    
    my $o2 = ($tally2[7] + $tally2[11] + $tally2[15]) / 
	( $tally2[5] + $tally2[7] + $tally2[9] + $tally2[11] + 
	  $tally2[13] + $tally2[15]);
    my $o12 = ($tally2[13] + $tally2[14] + $tally2[15]) / 
	( $tally2[12] + $tally2[13] + $tally2[14] + $tally2[15] );
    my $r = $o12 / ( $o1 + $o2 - ( $o1 * $o2 ) );
    my $o3 = ($tally2[13] + $tally2[14] + $tally2[15] )/
	( $tally2[9] + $tally2[10] + $tally2[11] + 
	  $tally2[13] + $tally2[14] + $tally2[15]);

    my $o4 = ($tally2[13] + $tally2[14] + $tally2[15] )/
	( $tally2[5] + $tally2[6] + $tally2[7] + 
	  $tally2[13] + $tally2[14] + $tally2[15]);
    my $nx = ($tally2[3] / ( $o1 * $o2) );
    my $N = ( $tally2[13] + $tally2[14] + $tally2[15] ) / 
	( $o12 * $o3 * $o4 );
    my $MLL = ($N * ( 1 - $r));
    my $MLG = ($nx - ( $N * ( $r * ( 1 - $o3) * ( 1 - $o4))));
    printf "\tMLE(L)=%d\n\tMLE(G)=%d\n\tr=$r\n\tMLE(N)=%d\n\n",
    POSIX::ceil($MLL), POSIX::ceil($MLG), POSIX::ceil($N);
    $rename{$internal_node_id} ||= $internal_node_id;
    $rename{$internal_node_id} .= sprintf("_%d/%d/%.2f",
					  POSIX::ceil($MLL),
					  POSIX::ceil($MLG),$r);
}

while (my ($r,$v) = each %rename ) {
    my ($n) = $tree->find_node(-id => $r);
    $n->id("'$v'");
}
for my $node ( $tree->get_leaf_nodes ) {
    my $id = $node->id;
    $id =~ s/\'//g;
    my ($x) = split(/_/,$id);
    $node->id(sprintf("'%s_%s'",$reverse_phy{$x},$id));
}


my $out = Bio::TreeIO->new(-format =>'newick', -file => ">$outtre");
$out->write_tree($tree);
