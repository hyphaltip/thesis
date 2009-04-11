#!/usr/bin/perl -w

=head1 NAME
 tribe2count

=head1 DESCRIPTION

count the number of members of each tribe family.
Expects tribe 2 column data
=cut

use warnings;
use strict;
use List::Util qw(sum);

my @grp;
my %spn;

while(<>) {
    my ($id,$gene) = split;
    my ($sp,$gn) = split(/:/,$gene);
    $spn{$sp}++;
    $grp[$id]->{$sp}++;
}

my @nm = sort keys %spn;
print join("\t", 'cluster', @nm), "\n";
my $i=0;
for my $x ( @grp ) {    
    my $total = sum values %{$x};    
    next if $total == 1;
    print join("\t", $i++, map { $x->{$_} || 0 } @nm ), "\n";
}
