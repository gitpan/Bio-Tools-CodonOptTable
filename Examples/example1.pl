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
        for my $a (@$myCodons)
                {
                   print "Codon : ",$a->{'codon'},"\t";
                   print "Frequency : ",$a->{'frequency'},"\t";
                   print "Amino acid : ",$a->{'aa_name'},"\t";
                   print "RSCU Value: ",$a->{'rscu'},"\t";
                   print "RAC Value: ",$a->{'rac'},"\t";
                   print "\n";
                }
    }


