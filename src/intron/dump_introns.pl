#!/usr/bin/perl -w
use strict;
use Env;
use Bio::DB::GFF;
use Bio::SeqIO;
use Bio::PrimarySeq;

my $db;
my ($dsn,$user,$password);
$user = $USER;
my $DEBUG = 0;
my $outname;
my @features;
use Getopt::Long;
GetOptions(
	   'u|user:s'      => \$user,
	   'p|pass:s'      => \$password,
	   'd|db|dbname:s' => \$db,
	   'dsn:s'         => \$dsn,
	   'o|out:s'       => \$outname,
	   'f|features:s'  => \@features,
	   'v|verbose+'     => \$DEBUG,
	   );


if( ! defined $db ) {
    die("Must provide a valid bio::db::gff db name or directory with -d/-db/--dbname\n");
}

if( ! $user || ! $password || ! $dsn ) {
    my $confile = File::Spec->catfile($HOME,".my.cnf");
    if( -e $confile ) {
        open(IN,$confile ) || die $!;
        my $ready = 0;
        while(<IN>) {
            if(/\[\s*(\S+)\s*\]/ ) {
                $ready = 1 if $1 eq 'mysql';
            } elsif( $ready ) {
                chomp;
                s/^\s+//;
                my ($key,$value) = split(/\s*\=\s*/,$_);
                if( $key eq 'user' ) {
                    $user = $value;
                } elsif( $key =~ /^pass/ ) {
                    $password = $value;
                } elsif( $key =~ /^host/ ) {
                    $dsn = "dbi:mysql:hostname=$value";
                }           
            }
        }
    }
    
}
$dsn .= ";database=$db";

my $dbh = Bio::DB::GFF->new(-dsn  => $dsn,
                            -user => $user,
                            -pass => $password,
                            -aggregators => [qw(coding match transcript)]);
if( ! @features ) {
    die("no features provided with -f\n");
} else {
    warn("requesting @features\n");
}
my $seqstream = $dbh->get_seq_stream(-type => $features[0]);
my $out;
if( $outname ) {
    $out = Bio::SeqIO->new(-file   => ">$outname",
                          -format => 'fasta');
} else {
    $out = Bio::SeqIO->new(-file   => ">$db.introns.fa",
			   -format => 'fasta');

}

my @genes;
F: while( my $feature = $seqstream->next_feature ) {    
    my $group = $feature->group;
    my $i = 1;
    next unless $feature->CDS;
    warn("$group $feature\n") if $DEBUG;

    my ($fstrand,$mixed);
    my @locs =   map { $_->[0] }
		   # sort so that most negative is first basically to order
		   # the features on the opposite strand 5'->3' on their strand
		   # rather than they way most are input which is on the fwd strand
		   
		   sort { $a->[1] <=> $b->[1] } # Yes Tim, Schwartzian transformation
		   map { 
		       $fstrand = $_->strand unless defined $fstrand;
		       $mixed = 1 if defined $_->strand && $fstrand != $_->strand;
		       [ $_, $_->start * ($_->strand || 1)];	    
		   } $feature->CDS ;
    my $lastloc;
    for my $cds (  @locs ) {
	if( $lastloc ) {
	    warn("$lastloc $cds\n") if $DEBUG;
	    my ($istart,$iend);
	    if( $cds->strand > 0 ) {
		($istart,$iend) = ( $lastloc->end+1,
				    $cds->start-1);
	    } else {
		($iend,$istart) = ( $cds->start+1,
				    $lastloc->end-1);
	    }
	    my $intron = $dbh->segment($feature->seq_id, $istart => $iend);
	    my $iseq = Bio::PrimarySeq->new(-seq => $intron->seq,
					    -display_id => "$group\_i$i",
					    -description => 
					    sprintf("%s:%s",
						    $feature->seq_id,
						    $istart < $iend ?
						    "$istart..$iend" :
						    "complement($iend..$istart)"
						    ));
	    $out->write_seq($iseq);
	    $i++;
	}
	$lastloc = $cds;
    }
}
