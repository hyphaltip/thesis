#!/usr/bin/perl -w
use strict;
use File::Spec;
use Env;

use strict;

use GO::TermFinder;
use GO::AnnotationProvider::AnnotationParser;
use GO::OntologyProvider::OntologyParser;

use GO::TermFinderReport::Text;

use GO::Utils::File    qw (GenesFromFile);
use GO::Utils::General qw (CategorizeGenes);

$|=1;

my $godir = "$HOME/lib/GO";

my $process   = GO::OntologyProvider::OntologyParser->new(ontologyFile => "$godir/process.ontology");
my $component = GO::OntologyProvider::OntologyParser->new(ontologyFile => "$godir/component.ontology");
my $function  = GO::OntologyProvider::OntologyParser->new(ontologyFile => "$godir/function.ontology");

my $assocfile = "$HOME/projects/H99/domains/all.associations";

my @files = @ARGV;

# now get our annotation file and number of genes
# now set up the objects we need



my $annotation = GO::AnnotationProvider::AnnotationParser->new(annotationFile => $assocfile);

my $cutoff = 0.01;
my $report = GO::TermFinderReport::Text->new();

# now go through each file
foreach my $file (@ARGV){

    warn "Analyzing $file\n";
    my @genes;
    open(IN,$file) || die $!;
    my $totalNum = 0;
    my %spnm;
    while(<IN>) {
	my ($cid,$g) = split;
	push @{$genes[$cid]}, $g;
	$totalNum++;
	my ($spid) = split(/:/,$g);
	$spnm{$spid} = 1;
    }
    my $termFinderP = GO::TermFinder->new(annotationProvider=> $annotation,
					  ontologyProvider  => $process,
					  totalNumGenes     => $totalNum,
					  aspect            => 'P');
	
    
    my $termFinderC = GO::TermFinder->new(annotationProvider=> $annotation,
					  ontologyProvider  => $component,
					  totalNumGenes     => $totalNum,
					  aspect            => 'C');
	
    my $termFinderF = GO::TermFinder->new(annotationProvider=> $annotation,
					  ontologyProvider  => $function,
					  totalNumGenes     => $totalNum,
					  aspect            => 'F');

    my @sp = sort keys %spnm;
    my $outfile = $file.".GO_terms";
    open (OUT, ">$outfile") || die "Cannot make $outfile : $!"; 
    open (WARN, ">$outfile.warn") || die "Cannot make $outfile.warn : $!"; 
    print OUT join("\t", qw(CLUSTER), @sp, qw(GO_P GO_F)), "\n";
    my $cid = 0;
    for my $cluster ( @genes ) {
	$cid++, next unless defined $cluster;

	my %gene_by_sp;
	for my $genename ( @$cluster ) {
	    my ($sp) = split(/:/,$genename);
	    $gene_by_sp{$sp}++;
	}
	my (@list, @notFound, @ambiguous);    
	CategorizeGenes(annotation  => $annotation,
			genes       => $cluster,
			ambiguous   => \@ambiguous,
			unambiguous => \@list,
			notFound    => \@notFound);

	if (@list) {
	    #print WARN (scalar @list), " gene(s) will be considered:\n";
#	    foreach my $gene (@list){
#		print OUT $gene, "\t", $annotation->standardNameByName($gene), "\n";		
#	    }
#
#	print OUT "\n";
	} else{
	    print WARN "None of the gene names were recognized in $cid\n";
#	    print WARN "They were:\n\n";
#	    print WARN join("\n", @notFound), "\n";
	    # close WARN;
	    next;
	}

	if (@ambiguous) {
	    print WARN scalar @ambiguous, " gene(s) are ambiguously named and unused in $cid\n";
	    
#	    print WARN join("\n", @ambiguous), "\n\n";
	}
	if (@notFound) {
	    # print WARN "The following gene(s) were not recognized, and will not be considered:\n\n";
	    # print WARN join("\n", @notFound), "\n\n";
	}
	my %go_seen;
	for my $termFinder ($termFinderP, $termFinderC, $termFinderF){
	    #print OUT "Finding terms for ", $termFinder->aspect, "\n\n";
	    my @pvalues = $termFinder->findTerms(genes        => \@list,
						 calculateFDR => 1);
	    
	    my $numHypotheses = $report->print(pvalues  => \@pvalues,
					       numGenes => scalar(@list),
					       totalNum => $totalNum,
					       cutoff   => $cutoff,
					       fh       => \*WARN);
	    
	    # if they had no significant P-values
	    if ($numHypotheses == 0){
		print WARN "No terms were found for this aspect with a corrected P-value <= $cutoff.\n";
	    } else {
		my @lastresort;
		for my $pv ( @pvalues ) {
		    if( $pv->{CORRECTED_PVALUE} <= $cutoff) {
			my @paths = $pv->{NODE}->pathsToRoot;
			for my $path ( sort { scalar @$a <=> scalar @$b } @paths ) {
			    if( @$path < 3 ) {
				push @{$lastresort[@$path]}, [ $pv->{NODE}->term, sprintf("p=%.3f",$pv->{CORRECTED_PVALUE})];
			    } else {
				$go_seen{$termFinder->aspect}->{$pv->{NODE}->term} = sprintf("p=%.3f",$pv->{CORRECTED_PVALUE});
			    }
			    last;
			}
		    }
		}
		if( ! exists $go_seen{$termFinder->aspect} && @lastresort) {
		    # just report the deepest node(s)
		    $go_seen{$termFinder->aspect}->{$_->[0]} = $_->[1] for @{$lastresort[-1]};
		}
	    }
#	    print WARN "\n\n";
	}	
	print OUT join("\t", $cid, (map { $gene_by_sp{$_} } @sp), 
		       (map { my $x = $_;
			      join(";", 
				   map { sprintf("%s:%s",
						 $_, $go_seen{$x}->{$_} ) }
				   keys %{ $go_seen{$x} }) } qw(P F) ),
		       ), "\n";
	$cid++;
    }
    close (OUT);
    close (WARN);
}

=pod

=head1 NAME

fam2GO.pl - gene family GO analysis

=head1 SYNOPSIS


=head1 AUTHORS

Jason Stajich, jason.stajich@duke.edu

=cut


# This script is based on the following

# Date   : 16th October 2003
# Author : Gavin Sherlock

# License information (the MIT license)

# Copyright (c) 2003 Gavin Sherlock; Stanford University

# Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated documentation files
# (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge,
# publish, distribute, sublicense, and/or sell copies of the Software,
# and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:

# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT.  IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
# BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
# ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
