#!/usr/bin/perl

use Bio::Tools::CodonOptTable;
use Data::Dumper;

#get from NCBI
my $seqobj = Bio::Tools::CodonOptTable->new(-ncbi_id => "J00522",
					    -genetic_code => 1,
					   );

my $myCodons = $seqobj->rscu_rac_table();
if($myCodons)
    {
        for my $each_aa (@$myCodons)
        {
            print "Codon      : ",$each_aa->{'codon'},"\t";
            print "Frequency  : ",$each_aa->{'frequency'},"\t";
            print "AminoAcid  : ",$each_aa->{'aa_name'},"\t";
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
    
# to generate a Graph between RSCU & RAC
$seqobj->generate_graph($myCodons,"myoutput.gif");


