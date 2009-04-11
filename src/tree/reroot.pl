#!/usr/bin/perl -w
use strict;
use Bio::TreeIO;
use Getopt::Long;

my ($oformat,$format,@outgroup) = qw(newick newick);
GetOptions(
	   'of|oformat:s' => \$oformat,
	   'if|iformat:s' => \$format,
	   'o|outgroup:s' => \@outgroup,
	   );
my $in = Bio::TreeIO->new(-fh => \*ARGV, -format => $format);
my $out= Bio::TreeIO->new(-format => $oformat);
while( my $t = $in->next_tree ) { 
    if( @outgroup ) {
	my @tips;
	for ( @outgroup ) {
	    my $n = $t->find_node($_);
	    if( ! $n ) {
		warn("could not find node $_\n");
	    } else {
		push @tips, $n;
	    }
	}
	while( @tips > 1 ) {
	    my $lca = $t->get_lca(-nodes => [shift @tips, shift @tips]);
	    unshift @tips, $lca;
	}
	my ($new_root) = @tips;
	if( $new_root->internal_id == $t->get_root_node->internal_id ) {
	    warn("rerooting with existing root\n");
	} else {
	    $t->reroot($new_root);
	}
    }
    $out->write_tree($t);
}

