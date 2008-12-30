#!/usr/bin/perl

use Bio::Tools::CodonOptTable;
use Data::Dumper;

my $seqobj = Bio::Tools::CodonOptTable->new ( -seq => 'ATGGGGTGGGCACCATGCTGCTGTCGTGAATTTGGGCACGATGGTGTACGTGCTCGTAGCTAGGGTGGGTGGTTTG',
				   -id  => 'GeneFragment-12',
				   -accession_number => 'X78121',
				   -alphabet => 'dna',
				   -is_circular => 1
				   );

my $myCodons = $seqobj->rscu_rac_table();
if($myCodons)
    {
        for my $each_aa (@$myCodons)
	{
	    print "Codon      : ",$each_aa->[1]->{'codon'},"\t";
	    print "Frequency  : ",$each_aa->[1]->{'frequency'},"\t";
	    print "AminoAcid  : ",$each_aa->[1]->{'aa_name'},"\t";
	    print "RSCU Value : ",$each_aa->[1]->{'rscu'},"\t"; #Relative Synonymous Codons Uses
	    print "RAC Value  : ",$each_aa->[1]->{'rac'},"\t"; #Relative Adaptiveness of a Codon
	    print "\n";
	}
    }
# to generate a Graph between RSCU & RAC
$seqobj->generate_graph($myCodons,"myoutput.gif");


