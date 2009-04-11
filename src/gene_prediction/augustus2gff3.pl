#!/usr/bin/perl -w
use strict;
use List::Util qw(min max);
print "##gff-version 3\n";
my $CDS = 'CDS';
my $prefix = 'Gene:';
my @groups_order;
my %groups;
while(<>) {    
    next if (/^\#/);    
    chomp;
    my @line = split(/\t/,$_);
    next unless $line[2] =~ /^$CDS$/;
    my $group = pop @line;
    $group =~ s/gene_id\s+\S+;//;
    $group =~ s/transcript_id\s+\"(\S+)\";/$line[0]-$1/;
    if( ! defined $groups{$group} ) {	
	push @groups_order, $group;
    }

    push @{$groups{$group}->{'CDS'}}, join("\t", @line,"Parent=$prefix$group")."\n";
    if( ! defined $groups{$group}->{'min'} ||
	$groups{$group}->{'min'} > $line[3] ) {
	$groups{$group}->{'min'} = $line[3];
    } 
    if( ! defined $groups{$group}->{'max'} ||
	$groups{$group}->{'max'} < $line[4] ) {
	$groups{$group}->{'max'} = $line[4];
    }
    if( ! defined $groups{$group}->{'strand'} ) {
	$groups{$group}->{'seq_id'} = $line[0];
	$groups{$group}->{'strand'} = $line[6];
    }
}

for my $x ( @groups_order ) {
    print join("\t", ( $groups{$x}->{'seq_id'},
		       'AUGUSTUS', 'mRNA', 
		       $groups{$x}->{'min'},
		       $groups{$x}->{'max'},
		       '.',
		       $groups{$x}->{'strand'},
		       '.',
		       "ID=$prefix$x")), "\n";
    for my $cds ( @{$groups{$x}->{'CDS'}} ) {
	print $cds;
    }
}
