#!/usr/bin/perl -w
use strict;
use List::Util qw(min sum);
use Getopt::Long;
use File::Spec;
use Env;
my $cutoff = '0.01';
my $Min_rate = '0.05'; # must be at least 5%
my %info;
my $outgroup;
my $descfile = File::Spec->catfile($HOME, qw(projects coprinus genome_analysis 
					     pfam_18.desc));
GetOptions(
	   'd|desc:s' => \$descfile,
	   'cutoff'   => \$cutoff,
	   'min'      => \$Min_rate,
	   'outgroup:s' => \$outgroup,
	   );

my $famfile = shift @ARGV;
if( $famfile =~ /\.gz$/ ) {
    open(FAM, "zcat $famfile |") || die "$famfile: $!";
} else {
    open(FAM, $famfile) || die "$famfile: $!";
}

for my $file (@ARGV) {
    if( $file =~ /\.gz/ ) {
	open(IN, "zcat $file |") || die "$file: $!";
    } elsif( $file =~ /\.bz2/ ) {
	open(IN, "bunzip2 -c $file |") || die "$file: $!";
    } else {
	open(IN,$file) || die "$file: $!";
    }
    while(<IN>) {
	my ($gene,$st,$en,$domain,$dst,$dend,$score,$evalue) = split;
	if( $info{$gene}->{$domain} ) {
	    $info{$gene}->{$domain} = $evalue if $evalue < $info{$gene}->{$domain};
	} else {
	    $info{$gene}->{$domain} = $evalue;
	}
    }
    close(IN);
}

my @grp;
my %spn;
my @cluster;
while(<FAM>) {
    my ($id,$gene) = split;
    my ($sp,$gn) = split(/:/,$gene);
    if( ! $sp ) {
	warn("$gene is $gene for $id\n");
	next;
    }
    if( my $g = $info{$gene} ) {
	
	# sort by evalue
	my @domains = sort { $g->{$b} <=> $g->{$a} } keys %{$g};
#	if( @domains > 1 ) {
#	    @domains = grep { $g->{$_} < $cutoff } @domains;
#	}
	$cluster[$id]->{$_}++ for @domains;
    }

    $spn{$sp}++;
    $grp[$id]->{$sp}++;
}
close(FAM);
open(DESC, $descfile) || die "$descfile: $!";
my %domaindesc;
while(<DESC>) {
    chomp;
    if(/^>(\S+)\|(\S+)\s+(.+)/ ) {
	$domaindesc{$1} = $3;
    }
}


my @nm = sort keys %spn;
print join("\t", 'cluster', @nm, "description"), "\n";
my $i=0;
my (@creations,@ex);
for my $x ( @grp ) { 
    my $total = sum values %{$x}; 
    $total |=0;
    $i++, next if $total < @nm;

    if( $outgroup ) {
	if(! $x->{$outgroup} ) {
	    push @creations, [$i,$x] if $total > 1;
	    $i++;
	    next; # skip where the outgroup is 0 (creation)
	} elsif( $x->{$outgroup} == $total  ) {
	    # skip where only the outgroup has a value (extinction);	    
	    push @ex, [$i,$x] if $total > 1;
	    $i++;
	    next;	    
	}
    }
    my @desc;
    for my $c ( sort { $cluster[$i]->{$b} <=> $cluster[$i]->{$a} } 
		keys %{$cluster[$i] || {}} ) {
	# must be in at least 1/10 of the seqs?
	next if( ($cluster[$i]->{$c} / $total) < $Min_rate);
	
	if( $domaindesc{$c} ) {
	    push @desc, $domaindesc{$c} . sprintf("(%d)",$cluster[$i]->{$c});
	}
    }
    print join("\t", $i++, (map { $x->{$_} || 0 } @nm), 
	       join("; ", @desc)), "\n";
}
$famfile =~ s/\.gz$//;
open(REPORT, ">$famfile.sumreport" ) || die $!;
printf REPORT "creations: %d\n",scalar @creations;
printf REPORT "extinctions: %d\n",scalar @ex;
close(REPORT);
open(REPORT, ">$famfile.creations" ) || die $!;
print REPORT join("\t", 'cluster', @nm, "description"), "\n";
for my $xx ( @creations ) {
    my ($i,$x) = @$xx;
    my @desc;
    my $total = sum values %{$x};     
    for my $c ( sort { $cluster[$i]->{$b} <=> $cluster[$i]->{$a} } 
		keys %{$cluster[$i] || {}} ) {
	# must be in at least 1/10 of the seqs?
	next if( ($cluster[$i]->{$c} / $total) < $Min_rate);
	
	if( $domaindesc{$c} ) {
	    push @desc, $domaindesc{$c} . sprintf("(%d)",$cluster[$i]->{$c});
	}
    }
    print REPORT join("\t", $i++, (map { $x->{$_} || 0 } @nm), 
	       join("; ", @desc)), "\n";
}
close(REPORT);
open(REPORT, ">$famfile.exinctions" ) || die $!;
print REPORT join("\t", 'cluster', @nm, "description"), "\n";
for my $xx ( @ex ) {
    my ($i,$x) = @$xx;
    my $total = sum values %{$x}; 
    my @desc;
    for my $c ( sort { $cluster[$i]->{$b} <=> $cluster[$i]->{$a} } 
		keys %{$cluster[$i] || {}} ) {
	# must be in at least 1/10 of the seqs?
	next if( ($cluster[$i]->{$c} / $total) < $Min_rate);
	
	if( $domaindesc{$c} ) {
	    push @desc, $domaindesc{$c} . sprintf("(%d)",$cluster[$i]->{$c});
	}
    }
    print REPORT join("\t", $i++, (map { $x->{$_} || 0 } @nm), 
	       join("; ", @desc)), "\n";
}
close(REPORT);
