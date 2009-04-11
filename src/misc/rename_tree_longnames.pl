#!/usr/bin/perl -w
use strict;
use Bio::TreeIO;

my %map = (
	   'atha' => 'Arabidopsis thaliana',
	   'ddis' => 'Dictyostelium discoideum',
	   'hsap' => 'Homo sapiens',
	   'rory' => 'Rhizopus oryzae',
	   'umay' => 'Ustilago maydis',
	   'cneo' => 'Cryptococcus neoformans',
	   'pchr' => 'Phanerochaete chrysosporium',
	   'ccin' => 'Coprinus cinereus',
	   'spom' => 'Schizosaccharomyces pombe',
	   'snod' => 'Stagonospora nodorum',
	   'anid' => 'Aspergillus nidulans',
	   'afum' => 'Aspergillus fumigatus',
	   'ater' => 'Aspergillus terreus',
	   'hcap' => 'Histoplasma capsulatum',
	   'cimm' => 'Coccidioides immitis',
	   'uree' => 'Uncinocarpus reesii',
	   'bcin' => 'Botrytis cinerea',
	   'sscl' => 'Sclerotinia sclerotiorum',
	   'fgra' => 'Fusarium graminearum',
	   'fver' => 'Fusarium verticillioides',
	   'mgri' => 'Magnaporthe grisea',
	   'pans' => 'Podospora anserina',
	   'cglo' => 'Chaetomium globosum',
	   'ncra' => 'Neurospora crassa',
	   'ylip' => 'Yarrowia lipolytica',
	   'dhan' => 'Debaryomyces hansenii',
	   'scer' => 'Saccharomyces cerevisiae',
	   'calb' => 'Candida albicans',
	   'cgla' => 'Candida glabrata',
	   'cgui' => 'Candida guilliermondii',
	   'clus' => 'Candida lusitaniae',
	   'ctro' => 'Candida tropicalis',
	   'klac' => 'Kluyveromyces lactis',
	   'agos' => 'Ashbya gossypii',
	   'scas' => 'Saccharomyces castellii',
	   'sbay' => 'Saccharomyces bayanus',
	   'spar' => 'Saccharomyces paradoxus',
	   'smik' => 'Saccharomyces mikatae',
	   'skud' => 'Saccharomyces kudriavzevii',
	   'sklu' => 'Saccharomyces kluyveri',
	   'kwal' => 'Kluyveromyces waltii',
	   'tree' => 'Trichoderma reesi',
	   'ptro' => 'Pan troglodytes',
	   'hsap' => 'Homo sapiens',
	   'dmel' => 'Drosophila melanogaster',
	   'cele' => 'Caenorhabdtis elegans',
	   'cbri' => 'Caenorhabdtis briggsae',
	   'rnor' => 'Rattus norvegicus',
	   'mmus' => 'Mus musculus',
	   'ggal' => 'Gallus gallus',
	   'drer' => 'Danio rerio',
	   'frub' => 'Fugu rubripes',
	   'tnig' => 'Tetraodon nigroviridis',
	   'cfam' => 'Canis familiaris',
	   'cint' => 'Ciona intestinalis',
	   'xtro' => 'Xenopus tropicalis',
	   'amel' => 'Apis mellifera',
	   'cdub' => 'Candida dubliniensis',
	   'agam' => 'Anopholes gambiae',
	   
	   );
my $in = Bio::TreeIO->new(-format => 'newick',
			  -file   => shift);
my $out = Bio::TreeIO->new(-format => 'newick');

while( my $tree = $in->next_tree ) {
    for my $node ( $tree->get_nodes ) {
	my $id = $node->id;
	next unless $id;
	my ($id_pref,@rest) = split(/_/,$id);
	my $name = $map{$id};
	unless( $name ) {
	    if( $name = $map{$id_pref} ) {
		if( @rest > 1 ) {
		    $name .= " $rest[-1]";
		}
	    } else {
		warn("no name for $id_pref\n") if @rest;
	    }
	}
	if( $name ) {
	    $name =~ s/\s/_/g;
	    $node->id($name);
	}

    }
    
    $out->write_tree($tree);
}
