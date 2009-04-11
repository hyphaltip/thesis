#!/usr/bin/perl -w
# $Id: fgblast,v 1.1 2005/06/02 19:22:42 jes12 Exp $

=head1 NAME

fgblast -- A script that will use a Blast-formatted database and a
gbrowse database and allow you to show graphicaly your blast hits

=head1 SYNOPSIS

This script allows you to display graphicaly the results of a blast search using a gbrowse database and some blast formated databases.


=head1 DESRIPTION

This script should be placed in your cgi-bin directory.

The script can then be called by going to the following address:

http://my.home.com/cgi-bin/gblast/yeast_chr

Assuming my.home.com is the hostname of your server, cgi-bin is your cgi bin directory where gblast resides, and yeast_chr 
is the database that you wish to blast against.

=head1 INSTALLATION

Place the gblast script in your cgi-bin directory.  Change any settings as described below.


Change any settings in the CONSTANTS portion of this script

  web_server		- This is the internet addressable hostname of this webserver (leave out http://)
  cgi_home    		- This is your cgi-bin directory where you placed this script
  gblast_bin		- The name of this script
  gbrowse_cgi_dir 	- This is the directory which holds the gbrowse executables
  gbrowse_database_post	- This is if you have a postfix that you place on your gbrowse databases
  orf_type		- The Class that you use for your open reading frames
  blast_db_dir		- The directory that holds your blast databases


You must then create your blast databases.  This script accounts for 4 different blast databases that are each prefixed with the name
of the gbrowse database.  For example, if my gbrowse database was called yeast_chr, you would create the following blast databaes

  yeast_chr			- A blast database of the nucleotide sequence of the contigs
  yeast_chr_orfs_nt		- A blast database of the nucleotides sequence of all of the open reading frames
  yeast_chr_orfs_aa		- A blast database of the amino acid sequence of all of the open reading frames
  yeast_chr_unused_reads_nt 	- A blast database of all of the reads that are not used (are not represented by the gbrowse database)

These files must be placed in the correct location as defined by the blast_db_dir variable.

The following subroutines must then be changed to reflect your system.

  get_orf_offset		- This will return the start position of the open reading frame related to it's reference sequence
  get_orf_contig		- This will return the reference sequence name that a particular orf resides on
  get_read_link			- This will return the link used for the unused_reads_nt database

=head1 AUTHOR

Michael Cipriano <mcipriano@mbl.edu>.

Updated by Jason Stajich <jason.stajich-at-duke.edu>

This is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut

use lib '/var/www/fungal/lib';
use CGI qw(:all);
use CGI::Pretty;
use DBI;
use Bio::Seq;
use Bio::SearchIO;
use Bio::SearchIO::Writer::HTMLResultWriter;
use Bio::Factory::FTLocationFactory;
use Bio::SeqIO;
use IO::String;
use File::Temp qw(tempfile);
use strict;
use warnings;

# CONSTANTS - CHANGE THESE TO YOUR PERSONAL SETTINGS
use constant { ImageWidth => 800,
	       web_server => 'http://fungal.genome.duke.edu',
	       HspSepSMax => 500,
	       HspSepQMax => 500,
	       cgi_home   => '/cgi-bin/',
	       gblast_bin => 'fgblast',
	       blast_exe_dir   => '/var/www/fungal/bin',
	       blast_db_dir    => '/var/www/fungal/blastdb',
	       gbrowse_cgi_dir => '/cgi-bin/',
	       db_conf_file    => '/var/www/fungal/conf/fgblast.conf',
	       tmpdir          => '/tmp/webblast',
	       DEBUG           => 0,
	   };

my $gblast_cgi = cgi_home . gblast_bin;
my $gbrowse_cgi = gbrowse_cgi_dir . 'gbrowse';
my $gbrowse_img = gbrowse_cgi_dir . 'gbrowse_img';

# Gbrowse DATABASE settings

# DO NOT CHANGE ANYTHING BELOW HERE UNLESS BROKEN


$ENV{'BLASTDB'}     = blast_db_dir;
$ENV{'BLASTFILTER'} = blast_exe_dir.'/filter';
$ENV{'BLASTMAT'}    = blast_exe_dir .'/matrix';

my $confh;
open($confh, db_conf_file) || die("cannot open ",db_conf_file,": $!");
my %dbs;
my %attributes;
my %reverse;
while(<$confh>) {
    chomp;
    my ($type,$name,$db) = split(/\t/,$_);
    $dbs{$db} = "$type $name";
    $reverse{"$type $name"} = $db;
}

close($confh);
print header();

if( !param('DATALIB') ) {
    print start_html(-title=> 'Fungal Genome BLASTs at Duke',
		     -head => [Link({-rel=>"stylesheet",
				     -type=>"text/css",
				     -href=>"/css/blast.css"})]),
    start_multipart_form( -name => 'MainBlastForm',
			  -action=> $gblast_cgi, 
			  -method=>'POST');
    
    print h3('Database and Program Options:'),p(),a({-href=>"/blast/docs/blast_program.html"},'Program'),
    popup_menu(-name    =>"PROGRAM", 
	       -values  =>[qw(blastp blastn blastx tblastn)],
	       -default =>param('program')),
    "Databases\n",    
    popup_menu(-name=>"DATALIB", 
	       -values  => [sort values %dbs],
	       -labels  =>\%dbs,
	       -default => param('database') );
    
print "\n",input({-type    =>'checkbox',
		  -value   =>"yes", 
		  -checked => 1 ,
		  -name    =>"show_img"}), 
    'Overlay Hits over Genome Image',
p(), 'Enter sequence below (most standard formats accepted but FASTA suggested)', 
    br, textarea({-name  =>"SEQUENCE",
		  -rows  =>12,
		  -cols  =>100}, param('seq'));
print '
<br>
Or load it from disk 
<p>
Set subsequence: From
&nbsp;&nbsp<input TYPE="text" NAME="QUERY_FROM" VALUE="" SIZE="10">
&nbsp;&nbsp&nbsp;&nbsp To
<input TYPE="text" NAME="QUERY_TO" VALUE="" SIZE="10">
<p>
<input TYPE="button" VALUE="Clear sequence" onClick="MainBlastForm.SEQUENCE.value="";MainBlastForm.QUERY_FROM.value="";MainBlastForm.QUERY_TO.value="";MainBlastForm.SEQUENCE.focus();">
<input TYPE="submit" VALUE="Search">
<hr>
 
The query sequence is 
<a href="http://www.ncbi.nlm.nih.gov/blast/html/blastcgihelp.shtml#filter">filtered</a> 
for low complexity regions by default.
<br>
<a href="http://www.ncbi.nlm.nih.gov/blast/html/blastcgihelp.shtml#filter">Filter</a>
 <input TYPE="checkbox" VALUE="L" NAME="FILTER" CHECKED> Low complexity</input>
<p>
 <input TYPE="checkbox" VALUE="Y" NAME="POSTSW" CHECKED> Post Process with Smith-Waterman (BLASTP)</input>
<p>

<a href="http://www.ncbi.nlm.nih.gov/blast/html/blastcgihelp.shtml#expect">Expect</a>
<select name ="EXPECT">
    <option> 1e-50</option>
    <option> 1e-25</option>
    <option> 1e-10</option>
    <option> 0.0001</option> 
    <option> 0.01 </option>
    <option selected> 1 </option>
    <option> 10 </option>
    <option>100 </option>
    <option>1000 </option>
</select>

<a href="http://www.ncbi.nlm.nih.gov/blast/html/blastcgihelp.shtml#Matrix">Matrix</a>
<select name ="MAT_PARAM">
    <option value ="PAM30"> PAM30 </option>
    <option value ="PAM70"> PAM70 </option> 
    <option value ="BLOSUM80"> BLOSUM80 </option>
    <option selected value ="BLOSUM62"> BLOSUM62 </option>
    <option value ="BLOSUM45"> BLOSUM45 </option>
</select>
<p>
</select>
<p>
<hr>
<p>
<input TYPE="button" VALUE="Clear sequence" onClick="MainBlastForm.SEQUENCE.value="";MainBlastForm.SEQUENCE.focus();">
<input TYPE="submit" VALUE="Search">
</form>
		</td>
</table>
<p><i>Powered by the <a href="http://blast.wustl.edu/">WU-Blast Programs</a> and <a href="http://www.bioperl.org">BioPerl</a>.</i></p>
		</td>
		<td width="4%"></td>
		<td width="48%" valign="top">
			<table border=0 width="100%">
			<tr>
				<td>

				
				</td>
			</tr>
			</table>		
		</td>
	</tr>

</table>
';

} else {
    my $prog = lc param('PROGRAM');
    my $exe = File::Spec->catfile(blast_exe_dir,$prog);
    my $stringio = IO::String->new(uc(param('SEQUENCE')));
    my $db = param('DATALIB');
    $db = $reverse{$db};
    my @params = ("-d \"$db\"",
		  '-cpus=1',
		  '-notes', '-warnings',
		  '-links', '-hspsepsmax='.HspSepSMax,
		  '-hspsepqmax='.HspSepQMax,
		  '-kap',
		  '-span1', '-topcomboN=5');
    
    my $writerhtml = Bio::SearchIO::Writer::HTMLResultWriter->new;
    $writerhtml->start_report(\&my_start_report);
    $writerhtml->title(\&my_title);
    $writerhtml->hit_link_align(\&my_hit_link_align);
    $writerhtml->hit_link_desc(\&my_hit_link_desc);
    if(param('FILTER'))
    {
	if($prog eq 'blastn' ) {
	    push @params, '-wordmask dust';
	} else {
	    push @params, '-wordmask seg+xnu';
	}
    }
    if( param('POSTSW') && $prog eq 'blastp') {
	push @params, '-postsw';
    }

    if(param('ALIGNMENTS')) {
	push @params, 'B='.param('ALIGNMENTS');
    }

    if(param('DESCRIPTIONS'))
    {
	push @params, 'V='.param('DESCRIPTIONS');
    }

    if(param('EXPECT'))
    {
	push @params, 'E='.param('EXPECT');
    }
#    if(param('MAT_PARAM'))
#    {
#	push @params, 'matrix='.param('MAT_PARAM');
#    }

    my $seqin;
    eval {
	$seqin = Bio::SeqIO->new( '-fh' =>$stringio);
    };
    
    if($@) {
	print h2(b("Sequence Not in recognized format"));
    } else {
	my $fh;
	mkdir(tmpdir) unless -d tmpdir;
	my ($tfh,$tmpfile) = tempfile(DIR => tmpdir, 
				      SUFFIX => '.blast', 
				      UNLINK => 1);
	my ($tfh_s,$tmpfile_s) = tempfile(DIR => tmpdir, 
					  SUFFIX => '.fa', 
					  UNLINK => 1);
	close($tfh);
	
	my $seqio = Bio::SeqIO->new(-format =>'fasta',
				    -fh     => $tfh_s);
	while(my $input = $seqin->next_seq()) {
	    $seqio->write_seq($input);
	}
	undef($seqio);
	close($tfh_s);
	my $cmd = "$exe -i $tmpfile_s -o $tmpfile ".join(" ",@params);
	#warn("run is $cmd\n");
	`$cmd >& /dev/null`;
	my $parser = Bio::SearchIO->new(-format => 'blast', 
					-file   =>$tmpfile);
	while( my $result = $parser->next_result ) {	
	    print $writerhtml->to_string($result);
	}
	# End IF Succedes    
    } # End If a recognized sequence
    close($stringio);
    undef $stringio;
}

	
sub my_start_report
{
    my $report = shift;
    my $str = <<EOF
<html>
<HEAD><title>Fungal Genome Resources at Duke -- BLAST</title>
<link REL="stylesheet" TYPE="text/css" href="/css/blast.css">
<link REL="icon" TYPE="image/ico" href="/favicon.ico">
</HEAD>
<body bgcolor="wheat">
EOF
;
    return $str;
}

sub my_title
{
    my $result = shift;
    return sprintf("<h2>%s Query of %s against %s</h2><br>", $result->algorithm,
		   $result->query_name,
		   param('DATALIB'));
}


sub my_hit_link_align
{
    my $self = shift;
    my $hit = shift;
    my $result = shift;
    my $min_val = 100000000000000;
    my $max_val = 0;

    my $hsp_string;
    my $db = param('DATALIB');
    my $offset = 0;
    my $contig_name;
    my $multiplier = 1;
    my (@h,$organism);
    if($db =~ /CDS/i)
    {
	($organism,$contig_name,
	 $offset) = &get_orf_contig_offset($hit->name(),
					   $hit->description);
	
    } elsif($db =~ /AA/i) {
	($organism,$contig_name,
	 $offset) = &get_orf_contig_offset($hit->name(),
					   $hit->description);
	$multiplier = 3;
    } else {

	$contig_name = $hit->name();
	if( $contig_name =~ s/(\S+):// ) {
	    $organism = $1;
	} else {
	    ($organism) = split(/_/,$contig_name);
	}
    }
    my %groups;
    while(my $hsp = $hit->next_hsp) {
	my $link = $hsp->links;
	$link =~ s/[\(\)]//g;
	push @{$groups{$link}}, $hsp;
    }
    my $linknum = scalar keys %groups;
    my @ret;
    while( my ($grp,$members) = each %groups ) {
	for my $hsp ( @$members ) {
	    my ($start,$end) = ($hsp->hit->start, $hsp->hit->end);
	    if($start > $end ) {
		warn("start is $start, end is $end for ", $hsp->hit->seq_id, "\n") if DEBUG;
	    }
	    $end   *= $multiplier;
	    $start *= $multiplier;
	    $hsp_string .= sprintf("+%d-%d",($start + $offset),( $end + $offset));
	    # Find area to center around
	    if( ($end+$offset) > $max_val) {
	    $max_val = $end + $offset;
	}
	    if( ($start+$offset) < $min_val) {
		$min_val = $start + $offset;
	    }
	}
	$min_val = $min_val - 100;
	$min_val = 1 if $min_val < 1;
	$max_val = $max_val +100;
	my $gbrowse_organism_cgi = $gbrowse_cgi .'/'.$organism;
#    my $return_string = '<a href="' . $gbrowse_organism_cgi . '?name=' . $contig_name . ':' . $min_val . '..' . $max_val . ';add=' . $contig_name . '+%22Blast%20Hit%22+Match' . $hsp_string . '">' ;
#    my $link_string = $return_string;
#    $return_string .=  $hit->name() . '</a>';;
#    if(param('show_img'))
#    {
#	my $gbrowse_organism_img = $gbrowse_img .'/'.$organism;
#	$return_string .= "<center>$link_string<img src=\"" . $gbrowse_organism_img . "?name=" . $contig_name .  ":$min_val..$max_val;width=800;add=$contig_name+%22Blast%20Hit%22+Match" . $hsp_string . '"></a></center>';
#	}
	
	my $url = sprintf("%s/%s?name=%s:%d..%d;add=%s+%%22Blast%%20Hit%%22+Match%s",
			  $gbrowse_cgi, $organism, $contig_name,
			  $min_val, $max_val,
			  $contig_name,$hsp_string);
	warn("$url for $organism $contig_name\n") if DEBUG;
	my $return_string = a({-href => $url},$hit->name);
	$return_string .=" Link_group:$grp" if $linknum > 1;
	if(param('show_img')) {
	    $return_string .= center(a({-href=>$url},
				       img({-class => 'noborder',
					    -src=> sprintf("%s/%s?name=%s:%d..%d;width=%d;add=%s+%%22Blast%%20Hit%%22+Match%s",
							   $gbrowse_img, 
							   $organism,$contig_name, 
							   $min_val,$max_val,
							   ImageWidth,
							   $contig_name,
							   $hsp_string)})
				       ));
	}
	push @ret,$return_string.br();
    }
    return join("\n", @ret);
}

sub my_hit_link_desc {
    my($self, $hit, $result) = @_;
    my (@h,$organism,$contig_name,$gene);
    my $db = param('DATALIB');
    if($db =~ /CDS|AA/i) {
	my ($pref,$gene) = split(/:/,$hit->name());
	($organism,@h) = split(/_/,$pref);
	if( @h > 1 ) {
	    $organism .= "_$h[0]";
	}
	$contig_name = $gene;
    } else {
	$contig_name = $hit->name();
	if( $contig_name =~ s/(\S+):// ) {
	    $organism = $1;
	} else {
	    ($organism) = split(/_/,$contig_name);
	}
    }

    my $url = sprintf("%s/%s?name=%s",
		      $gbrowse_cgi, $organism, $contig_name, 
		      );
    return a({-href=>$url},$hit->name);
}

sub get_orf_contig_offset
{
    # This must be changed to reflect your method of finding out the offset of the orf start position
    # Mine is encoded in the DESCRIPTION part of the FASTA header
    # >ORFNAME CHROM:LOCATIONSTRING cdslen=CDSLEN
    # Like this 
    # >mgri_brd:MG03167.4 mgri_2.631:join(complement(47648..47743),complement(45977..47434)) cdslen=1767
 
    my $orfid = shift;
    my $orfdesc = shift;    
    my ($prefix,$gene) = split(/:/,$orfid);
    my ($organism,$strain,$src) = split(/_/,$prefix);
    if( $src ) {
	$organism .= "_$strain";
    } else { $src = $strain }

    $orfdesc =~ s/^\s+//;
    $orfdesc =~ s/\s+cdslen\S+//;
    $orfdesc =~ s/\s+//g;
    # here I get the orf location from the header, parse it so I can map ORF/AA location
    # to contig location
    my ($contig_loc) = split(/\s+/,$orfdesc);
    my ($contig,$location) = split(/:/,$contig_loc);
    my $offset = Bio::Factory::FTLocationFactory->from_string($location)->start;
    return ($organism,$contig,$offset);
}

sub get_gbrowse_dbname {
    my $hitname = shift;
    my ($prefix) = split(/:/,$hitname,2);
    my ($sp,$strain,$src) = split(/_/,$prefix);
    if( $src ) {
	$sp .= "_$strain";
    }
    $sp;
}
