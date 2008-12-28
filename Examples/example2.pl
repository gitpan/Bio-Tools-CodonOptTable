#!/usr/bin/perl

use Bio::Tools::CodonOptTable;
use Data::Dumper;

#read from file
my $seqobj = Bio::Tools::CodonOptTable->new(-file => "contig.fasta",
                                         -format => 'Fasta');

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


