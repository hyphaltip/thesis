#!/usr/bin/perl -w
# This is based on blast2table.pl code 
# which was written by
# Korf, Yandell, Bedell as part of the BLAST Book O'Reilly and Associates 2003
use strict;
use Getopt::Std;
my @FIELDS = qw(percent q_align mismatch s_align q_begin q_end 
		s_begin s_end expect bits similar);
use vars qw($opt_p $opt_b $opt_e $opt_m $opt_n);
getopts('p:b:e:m:n:');
my $PERCENT = $opt_p ? $opt_p : 0;
my $BITS    = $opt_b ? $opt_b : 0;
my $EXPECT  = $opt_e ? $opt_e : 1e30;
my $START   = $opt_m ? $opt_m : 0;
my $END     = $opt_n ? $opt_n : 1e30;

my ($Query, $Qlen,$Sbjct,$Slen);
my $HSP = "";
while (<>) {
    if    (/^Query=\s+(\S+)/) {outputHSP(); $Query = $1}
    elsif( /\((\-?[\d,]+)\s+letters.*\)/ ) { $Qlen = $1; $Qlen =~ s/,//g; }
    elsif( /Length\s*=\s*(\-?[\d,]+)/ ) { $Slen =$1; $Slen =~ s/\,//g; }
    elsif (/^>(\S+)/)         {outputHSP(); $Sbjct = $1}
    elsif (/^ Score = /) {
        outputHSP();	
        my $stats = $_;
        while (<>) {
            last unless /\S/;
	    if(/Length\s*=\s*(\-?[\d,]+)/ ) {
		$Slen = $1;
		$Slen =~ s/\,//g;
	    }
            $stats.= $_;
        }
	my ($score) = ($stats =~ /Score\s+=\s+(\-?\d+(\.\d+)?)/);
        my ($bits) = $stats =~ /(\d\S+) bits/;
        my ($expect) = $stats =~ /Expect\S* = ([\d\.\+\-e]+)/;
		$expect = "1$expect" if $expect =~ /^e/;
        my ($match, $total, $percent)
            = $stats =~ /Identities = (\d+)\/(\d+)\s+\((\d+)%\)/;
        my $mismatch = $total - $match;
        my ($similar) = ($stats =~ /Positives = (\d+)\/(\d+)\s+\((\d+)%\)/);
	
        $HSP = {score => $score,
		bits => $bits, expect => $expect, mismatch => $mismatch,
		percent => $percent, q_begin => 0, q_end => 0, q_align => 0,
		s_begin => 0, s_end => 0, s_align => "",
		similar => $similar };
    }
    elsif( /^(Query|Sbjct):\s+(\d+)\s+\-\s+(\d+)/) {

	my $type = $1;
	my $fchar = lc substr($type,0,1);
        $HSP->{$fchar.'_begin'}  = $2 unless $HSP->{$fchar.'_begin'};
        $HSP->{$fchar.'_end'}    = $3;
        $HSP->{$fchar.'_align'} .= '';
    } elsif (/^(Query|Sbjct):\s+(\d+)\s+(\S+)\s+(\d+)/) {
	my $type = $1;
	my $fchar = lc substr($type,0,1);
        $HSP->{$fchar.'_begin'}  = $2 unless $HSP->{$fchar.'_begin'};
        $HSP->{$fchar.'_end'}    = $4;
        $HSP->{$fchar.'_align'} .= $3;
    }
}
outputHSP();

sub outputHSP {
    return unless $HSP;
    return if $HSP->{percent}  < $PERCENT;
    return if $HSP->{bits}     < $BITS;
    return if $HSP->{expect}   > $EXPECT;
    return if ($HSP->{q_begin} < $START or $HSP->{q_end} < $START);
    return if ($HSP->{q_begin} > $END   or $HSP->{q_end} > $END);
    foreach my $field ( @FIELDS ) {		      
    	print "$field not defined\n" if not defined $HSP->{$field};
    }
    my $qgaps = countGaps($HSP->{q_align});
    my $sgaps = countGaps($HSP->{s_align});
    print join("\t", $Query, $Sbjct, $HSP->{percent},
	       length($HSP->{q_align}), $HSP->{mismatch},
	       $qgaps + $sgaps,
	       $HSP->{q_begin}, $HSP->{q_end}, $HSP->{s_begin}, $HSP->{s_end},
	       $HSP->{expect}, $HSP->{bits},
	       $HSP->{score},'',
	       $HSP->{similar},
	       $Qlen, $Slen,
	       $qgaps,$sgaps), "\n";
    $HSP = "";
}

sub countGaps {
    my $count = 0;
    my $pos =0;
    while(($pos = index($_[0],'-',$pos+1)) >= 0) {
	$count++;
    }
    return $count;
}
