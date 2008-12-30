#!/usr/bin/perl

use Bio::Tools::CodonOptTable;
use Data::Dumper;

#get from NCBI
my $seqobj = Bio::Tools::CodonOptTable->new(-ncbi_id => "J00522");

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


