#!/usr/bin/perl -w

=head1 NAME

make_splicesite - make a genbank file for AUGUSTUS

=head1 SYNOPSIS

   -d,--db            databasename or fa dir [required]
   -s/--species       species name or abbreviation [required]
   -v,--version       GFF version [default 3]
   -h/--help          display this usage help

=cut

use strict;
use Getopt::Long;
use List::Util qw(shuffle);
use Env;
use Bio::DB::GFF;
use Bio::Tools::GFF;
use Bio::SeqIO;
use Bio::Seq;
use Bio::SeqFeature::Generic;

use constant CODONLEN => 3;

my $REGCT = 1;

my ($species,$db,@gff);
my $gversion = 3;
my $window = 1000;
my $DEBUG = 0;
my $MinIntron = 5;
my ($dsn,$user,$password);
$user = $USER;
my @features;
GetOptions('d|db|dbname:s'    => \$db,
	   'u|user:s'         => \$user,
	   'p|pass|password:s'=> \$password,
	   'dsn:s'            => \$dsn,
	   'v|version:s'      => \$gversion,
	   's|species:s'      => \$species,
	   'w|window:i'       => \$window,
	   'mi:i'             => \$MinIntron,
	   'debug'            => \$DEBUG,
	   'f|feature:s'      => \@features,
	   'h|help'           => sub { system('perldoc', $0);
				       exit;
				   },
	   );

if( !defined $species ) {
    die("Must provide a species name (abbrev) with -s/--species\n");
}

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
my $out = Bio::SeqIO->new(-file   => ">$species.train.gb",
			  -format => 'genbank');

my $out_test = Bio::SeqIO->new(-file   => ">$species.test.gb",
			       -format => 'genbank');

my @genes;
F: while( my $feature = $seqstream->next_feature ) {    
    my $group = $feature->group;
    my ($start,$end)= sort { $a <=> $b } ( $feature->start,
					   $feature->end);
    $start -= $window;
    $start = 1 if $start < 1;
    $end += $window;
    my $length = $dbh->segment($feature->seq_id)->length;
    $end = $length if $end > $length;
    ($start,$end) = ($end,$start) if $feature->strand < 0;
    my $segment = $dbh->segment($feature->seq_id, $start => $end);    
    my @features = $segment->features;

    my $seq = Bio::Seq->new(-seq        => $segment->seq,
			    -display_id => "region-$REGCT",
			    -description => sprintf("%s:%d..%d",
						    $segment->seq_id,
						    $segment->start,
						    $segment->end));
    $seq->add_SeqFeature(Bio::SeqFeature::Generic->new
			 (-start       => 1,
			  -end         => $seq->length,
			  -primary_tag => 'source',
			  -strand      => 1));
    my $seqid    = $seq->display_id;
    my $slength  = $seq->length;
    for my $f ( @features ) {
	next unless $f->group eq $group;
	next if( $f->primary_tag eq 'gene');
	$f->ref($segment);	
	if( $f->primary_tag eq 'coding' || $f->primary_tag eq 'transcript') {
	    
	    my ($fstrand,$mixed);
	    my @locs = map { $_->[0] }
	    # sort so that most negative is first basically to order
	    # the features on the opposite strand 5'->3' on their strand
	    # rather than they way most are input which is on the fwd strand
	    
	    sort { $a->[1] <=> $b->[1] } # Yes Tim, Schwartzian transformation
	    map { 
		$fstrand = $_->strand unless defined $fstrand;
		$mixed = 1 if defined $_->strand && $fstrand != $_->strand;
		[ $_, $_->start * ($_->strand || 1)];	    
	    } $f->location->each_Location;
	    my $split = Bio::Location::Split->new(-seq_id => $seqid);
	    for my $loc ( @locs ) {
		my ($start,$end) = sort { $a <=> $b } ( $loc->start,$loc->end);
		next F if( $start == $end || 
			   ! defined $start || 
			   ! defined $end );
		$split->add_sub_Location(  Bio::Location::Simple->new
					   (-start => $start,
					    -end   => $end,
					    -strand=> 1));
	    }
	    my $groupstr = $group;
	    $groupstr =~ s/^Gene://;
	    $f = Bio::SeqFeature::Generic->new(-location    => $split,
					       -primary_tag => 'CDS',
					       -seq_id      => $seqid,
					       -tag => { 'gene' => $groupstr});
	    $seq->add_SeqFeature($f);
	}
    }
						       
    push @genes, $seq;
    $REGCT++;
    last if $DEBUG && $REGCT > 10;
}
@genes = shuffle @genes;
warn("Found ",scalar @genes, " genes\n");
for my $g ( splice(@genes,0,scalar @genes/2) ) {    
    $out->write_seq($g);
}

for my $g ( @genes ) {
    $out_test->write_seq($g);
}
