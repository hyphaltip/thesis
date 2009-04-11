#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $allowed_missing = 0;
my $DEBUG = 0;
GetOptions(
	   'a|am|allowed_missing:i'  => \$allowed_missing,
	   'v|verbose|d|debug'       => \$DEBUG,
         );
my %pairs;
my %total_sp;
while(<>) {
    my ($a,$b) = split;
    my ($spa) = split(/:/,$a);
    my ($spb) = split(/:/,$b);
    $total_sp{$spa}++; $total_sp{$spb}++;
    $pairs{$a}->{$b} = 1;
    $pairs{$b}->{$a} = 1;
}

my @genes = keys %pairs;
my (%lookup,%done);
my @toprint;
my $sp_count = scalar keys %total_sp;
warn("$sp_count species\n");
GENE: for my $g ( @genes ) {
    next if $done{$g};
    my @orthologs = keys %{ $pairs{$g} };
    my @extra;
    my %s;

    for ( @orthologs ) {
	push @extra, keys %{ $pairs{$_} };
    }    
    @orthologs = uniq ( $g, @extra, @orthologs );
    my %t;
    for my $gn ( @orthologs ) { 
	my ($sp,$gene) = split(/:/,$gn);
	next GENE if( $t{$sp}++ );  # if a species is present more
        # than once or gene already used bail	
	if( $done{$gn} ) {
	    if( $lookup{$gn} ) {
		# warn("lookup was ", $lookup{$gn}, "\n");
		$toprint[$lookup{$gn}] = '';
	    }
	    next GENE;
	}
    }
    my $ct = keys %t;
    if( $ct <= $sp_count && 
	$ct >= ($sp_count - $allowed_missing) ) {
	push @toprint, join("\t", sort (@orthologs)). "\n";
	for (@orthologs) { 
	    $lookup{$_} = $#toprint;
	    $done{$_} = 1;
	}
    } elsif( $ct < $sp_count ) {
	my %missing = map { $_ => 1} keys %total_sp;
	for my $g ( (@orthologs) ) {
	    my ($sp,$gn) = split(/:/,$g,2);
	    $missing{$sp} = 0;
	}
#	warn("missing: ", join(",",sort grep { $missing{$_} } keys %missing)," out of $g matches\n");
	# 
    } elsif( $ct > $sp_count ) {
	my %seen;
	for my $g ( @orthologs ) {
	    my ($sp,$gn) = split(/:/,$g,2);
	    push @{$seen{$sp}}, $g;
	} 
	my @too_many;
	for my $x ( sort keys %seen ) {
	    if( @{$seen{$x}} > 1) {
		push @too_many, sprint("%d;%s",scalar @{$seen{$x}},
				       join(",",@{$seen{$x}}));
	    }
	}
	warn("too_many ",join(" ", @too_many),"\n");
    } else {
	warn("ct is $ct sp_count is $sp_count\n");
    }
}
print for @toprint;

sub uniq {
    my %seen;    
    grep { ! $seen{$_}++ } @_;
}
