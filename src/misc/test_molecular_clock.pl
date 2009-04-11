#!/usr/bin/perl -w
use strict;
my @dat;
use Statistics::Distributions;

my $signif_high = 0.01;
my $signif_mod  = 0.05;
# $type would be either 'clock' or 'noclock' depending on which file you are reading.
# so you'd want this type of loop to run over all the result files from PAML
# in my case I've numbered all the gene clusters and that number is stored in $num
# designed to parse PAML mlc files from codeml
for my $file ( @ARGV ) {
    my $type;
    if( $file =~ /no_clock/ ) {
	$type = 'noclock';
    } elsif( $file =~ /clock/ ) {
	$type = 'clock';
    } else { 
	die "unknown type $file";
    }
    my $num = 0;
    open(IN, $file) || die $!;
    while(<IN>) {
	if( /lnL\([^\)]+\):\s*(\-?\d+\.\d+)/ ) {
	    $dat[$num]->{$type} = $1 * -1; # this is just a lnL it is not -lnL
	} elsif( /^ns\s+=\s+(\d+)/ ) {
	    $dat[$num]->{'ns'} = $1;
	}
    }
    close(IN);
}
# to then compare the clock vs no-clock
my $i = 0;
print join("\t", qw(CLUSTER NOCLOCK CLOCK LRT SIGNIF)),"\n";
for my $d ( @dat ) {    
    if( defined $d ) {
        # H1,H0 = - log (likelihood)
         
        # LRT = 2*(lnL1 - lnL2);
        my $lrt = 2 * ($d->{'clock'} - $d->{'noclock'} );
        if( ! defined $d->{clock} || ! defined $d->{'noclock'} ) {
            warn("no dat for $i\n");
        }
        my $p = Statistics::Distributions::chisqrprob ($d->{'ns'} - 2,$lrt);
        printf "%s\t%9.3f\t%9.3f\t%-7.2f\t%-5.2g%s\n",$i,$d->{'noclock'},
        $d->{'clock'}, $lrt, $p, $p < $signif_high ? '**' : $p < $signif_mod ? '*' : '';
    }
    $i++;
}
