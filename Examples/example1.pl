#!/usr/bin/perl

use Bio::Tools::CodonOptTable;
use Data::Dumper;

my $seqobj = Bio::Tools::CodonOptTable->new ( -seq => 'ATGGGGTGGGCACCATGCTGCTGTCGTGAATTTGGGCACGATGGTGTACGTGCTCGTAGCTAGGGTGGGTGGTTTG',
				   -id  => 'GeneFragment-12',
				   -accession_number => 'X78121',
				   -alphabet => 'dna',
				   -is_circular => 1,
				   -genetic_code => 1,
				   );

my $myCodons = $seqobj->rscu_rac_table();

if($myCodons)
    {
		for my $each_aa (@$myCodons)
		{
			print "AminoAcid  : ",$each_aa->{'aa_name'},"\t";
			print "Codon      : ",$each_aa->{'codon'},"\t";
			print "Frequency  : ",$each_aa->{'frequency'},"\t";
			print "RSCU Value : ",$each_aa->{'rscu'},"\t"; #Relative Synonymous Codons Uses
			print "RAC Value  : ",$each_aa->{'rac'},"\t"; #Relative Adaptiveness of a Codon
			print "\n";
		}
    }

# To get the prefered codon list based on RSCU & RAC Values
my $prefered_codons = $seqobj->prefered_codon($myCodons);

    while ( my ($amino_acid, $codon) = each(%$prefered_codons) ) {
        print "AminoAcid : $amino_acid \t Codon : $codon\n";
    }

# To generate a Graph between RSCU & RAC
$seqobj->generate_graph($myCodons,"myoutput.gif");


