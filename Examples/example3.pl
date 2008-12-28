#!/usr/bin/perl

use Bio::Tools::CodonOptTable;
use Data::Dumper;

#get from NCBI
my $seqobj = Bio::Tools::CodonOptTable->new(-ncbi_id => "J00522");

my $myCodons = $seqobj->rscu_rac_table();
if($myCodons)
    {
        for my $a (@$myCodons)
                {
                   print "Codon : ",$a->{'codon'},"\t";
                   print "Frequency : ",$a->{'frequency'},"\t";
                   print "Amino_Acid : ",$a->{'aa_name'},"\t";
                   print "RSCU Value: ",$a->{'rscu'},"\t";
                   print "RAC Value: ",$a->{'rac'},"\t";
                   print "\n";
                }
    }


