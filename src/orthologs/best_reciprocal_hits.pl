#!/usr/bin/perl -w

# This requires FASTA output turned into TAB format with the
# fastam9_to_table script

use strict;
use Getopt::Long;
use DB_File;
use File::Spec;
use Getopt::Long;
use constant SEP => ';';
my $min_q_aligned = 10;
my $min_h_aligned = 10;
my $min_pid = 15;
$DB_BTREE->{'flags'} = R_DUP ;
my $idxdir = 'idx';
my $force = 0;
GetOptions(
	   'qaln|q:i' => \$min_q_aligned,
	   'haln|h:i' => \$min_h_aligned,
	   'p|pid:i' => \$min_pid,
	   'f|force' => \$force,
	   'idx:s'   => \$idxdir,
	   );
my $count = 0;
my $species_included = 1;
my $outdir = 'pairwise_orthologs_all';
mkdir($outdir) unless -d $outdir;
my ($lastq,%rank, %rankkept);

my (%pairs,%sp2num);
my $unique_spnum = 0;
mkdir($idxdir) unless -d $idxdir;
my $spfile   = File::Spec->catfile($idxdir,'species2num.idx');

if( $force ) { 
    unlink($spfile);
}

my $update = 1;
if( -e $spfile  && ! -z $spfile ) {
    $update = 0;
}

my $x2 = tie(%sp2num,'DB_File',$spfile); 

my (%dbhs);
if( $update ) {
    while( <> ) {
	next if( /^\#/ || /^\s+$/ );
	my ($qnamein,$hnamein, $percent_id, $hsp_len, $mismatches,$gapsm,
	    $qstart,$qend,$hstart,$hend,$evalue,$bits,$fasta_score,
	    $sw_score,undef,$percent_similar,$query_len,$hit_len,$query_gaps,
	    $hit_gaps) = split;
	
	my $q_aligned= (abs($qend - $qstart) / $query_len ) * 100;
	my $h_aligned= (abs($hend - $hstart) / $hit_len ) * 100;
	next unless ($q_aligned >= $min_q_aligned && 
		     $h_aligned >= $min_h_aligned && 
		     $percent_similar >= $min_pid);
	
	my ($qspecies,$qname,$hspecies,$hname);
	if( $species_included ) { 
	    ($qspecies,$qname) = split(/:/,$qnamein);
	    ($hspecies,$hname) = split(/:/,$hnamein);
	} else { 
	    ($qspecies,$hspecies) = (1,2);
	    ($qname,$hname) = ($qnamein,$hnamein);
	}
	next if ($qnamein eq $hnamein );
	if ( ($qspecies =~ /^(afum|sscl)$/ && $qname =~ /^GLEAN/ ) ||
	     ($hspecies =~ /^(afum|sscl)$/ && $hname =~ /^GLEAN/ ) ||
	     ($hspecies eq 'anid' && $hname !~ /\.2\.p\d+/ ) || 
	     ($qspecies eq 'anid' && $qname !~ /\.2\.p\d+/ ) 
	     ) {
	    #warn("skipping $qnamein $hnamein\n");
	    #next;
	}
	my $sp_id_nm = "$qspecies-$hspecies";
	$sp2num{$qspecies} = $unique_spnum++ unless defined $sp2num{$qspecies};
	$sp2num{$hspecies} = $unique_spnum++ unless defined $sp2num{$hspecies};
	my $dbh = $dbhs{$sp_id_nm};
	unless( defined $dbh ) {
	    $dbh = $dbhs{$sp_id_nm} = 
		tie(%{$pairs{$sp_id_nm}},
		    'DB_File', 
		    File::Spec->catfile($idxdir,$sp_id_nm),
		    O_RDWR|O_CREAT, 0666, $DB_BTREE);
	    if( ! $dbh ) {
		for my $v ( values %dbhs ) {
		    if( defined $v ) {
			$v->sync;
		    }
		}
		%dbhs = ();
		for my $v ( values %pairs ) {		    
		    untie %$v;
		}
		%pairs = ();
		$dbh = $dbhs{$sp_id_nm} = 
		    tie(%{$pairs{$sp_id_nm}},
			'DB_File', 
			File::Spec->catfile($idxdir,$sp_id_nm),
			O_RDWR|O_CREAT, 0666, $DB_BTREE) || 
			die("cannot initialize a dbh for $sp_id_nm\n");
	    }
	}
	if( ! defined $dbh ) { 
	    warn("no dbh for $sp_id_nm \n");
	    next;
	}
	if( defined $lastq && $lastq ne $qnamein ) { 
	    %rank = %rankkept = ();
	}
	$lastq = $qnamein;

	my $rank = $rankkept{$hname};
	unless ( defined $rank ) {
	    $rank = $rank{$hspecies}++;
	    $rankkept{$hname} = $rank;
	}
	if(  my @vals = $dbh->get_dup($qname) ) {
	    for my $h ( @vals ) {
		my ($hrank,$hn,$e,$len) = split(SEP,$h);
		if( $hn eq $hname ) {
		    if( $evalue < $e ) { 
			$e = $evalue; # save BEST hsp evalue for now, sum stats later on maybe
		    }		    
		    $dbh->del_dup($qname,$h);
		    $dbh->put($qname,join(SEP,($hrank,$hn,$e,
					       $len + $hsp_len)));
		    last;
		}
	    }
	} else { 
	    $dbh->put($qname,join(SEP,($rank,$hname,$evalue,$hsp_len)));
	}
	warn("$count\n") if ++$count % 250000 == 0;
    } 
}

my @species = keys %sp2num;
warn("species are ",scalar @species, ": @species \n");

for my $v ( values %dbhs ) {
    if( defined $v ) {
	$v->sync;
    }
}
%dbhs = ();
for my $v ( values %pairs ) {		    
    untie %$v;
}
%pairs = ();
for( my $i = 0; $i < scalar @species; $i++ ) {
    my $qspecies = $species[$i];
    for ( my $j = $i+1;$j < scalar @species; $j++ ) {
	my $hspecies = $species[$j];
	my $sp_id_nm = "$qspecies-$hspecies";
	
	my (%seen,%pairsf,%pairsr); 
	my $pairout;
	open($pairout, ">$outdir/$sp_id_nm.orthologs") || die($!);
	my $fwd_dbh = tie(%pairsf,
			  'DB_File', 
			  File::Spec->catfile($idxdir,$sp_id_nm),
			  O_RDWR, 0666, $DB_BTREE) ||
			      die("cannot initialize a dbh for $sp_id_nm\n");
	my $rev_dbh = tie(%pairsr,
			  'DB_File', 
			  File::Spec->catfile($idxdir,"$hspecies-$qspecies"),
			  O_RDWR, 0666, $DB_BTREE) ||
		 die("cannot initialize a dbh for $hspecies-$qspecies\n");
	Q: for my $qseq ( keys %pairsf ) {
	    my @pair;
	    my @hits = $fwd_dbh->get_dup($qseq);
	    next unless @hits && defined $hits[0];
	    H: for my $h ( @hits ) {
		my ($rank,$hit,$signif,$len) = split(SEP,$h);
		if( $rank == 0 ) {
		    next Q if ( $seen{$hit."||".$qseq}++ ||
				$seen{$qseq."||".$hit}++);
		    
		    my @rev_hits = $rev_dbh->get_dup($hit);
		    next unless ( @rev_hits && defined $rev_hits[0] );
		    for my $r ( @rev_hits ) {
			my ($revrank,$rev_name,$rev_signif,
			    $rev_hsp) = split(SEP,$rev_hits[0]);

			if( $revrank == 0 && 
			    $rev_name eq $qseq ) {
			    @pair = ("$hspecies:$hit",
				     sprintf("%.4e",$signif),$len);
			    last;
			}
			last;
		    }
		}
	    }
	    if( @pair ) {
		print $pairout join("\t","$qspecies:$qseq", @pair),"\n";  
	    }
	}
	close($pairout);
	$fwd_dbh = $rev_dbh = undef;
	untie(%pairsr);
	untie(%pairsf);
    }
}
