#!/usr/bin/perl -w
use strict;
use Bio::TreeIO;
use Getopt::Long;

my ($oformat,$format) = qw(newick nexus);
GetOptions(
	   'of|oformat:s' => \$oformat,
	   'if|iformat:s' => \$format,
	   );
my $in = Bio::TreeIO->new(-fh => \*ARGV, -format => $format);
my $out= Bio::TreeIO->new(-format => $oformat);
while( my $t = $in->next_tree ) { $out->write_tree($t) }
