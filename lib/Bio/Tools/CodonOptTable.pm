package Bio::Tools::CodonOptTable;

use warnings;
use strict;
use Data::Dumper;

use Bio::Root::RootI;
use Bio::Root::Root;
use Bio::PrimarySeq;
use Bio::SeqIO;
use Bio::Tools::SeqStats;
use Bio::Tools::CodonTable;
use Bio::DB::GenBank;

=head1 NAME

Bio::Tools::CodonOptTable - A more elaborative way to check the codons quality

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

use vars qw(@ISA %Amnioacid);

@ISA = ('Bio::Root::Root', 'Bio::SeqIO', 'Bio::PrimarySeq');
use base qw(Exporter);

our @EXPORT = qw(
		    new
                    rscu_rac_table
		);

my %Amnioacid = (
        'A'=>'Ala',
        'R'=>'Arg',
        'N'=>'Asn',
        'D'=>'Asp',
        'C'=>'Cys',
        'Q'=>'Gln',
        'E'=>'Glu',
        'G'=>'Gly',
        'H'=>'His',
        'I'=>'Ile',
        'L'=>'Leu',
        'K'=>'Lys',
        'M'=>'Met',
        'F'=>'Phe',
        'P'=>'Pro',
        'S'=>'Ser',
        'T'=>'Thr',
        'W'=>'Trp',
        'Y'=>'Tyr',
        'V'=>'Val',
        'B'=>'Asx',
        'Z'=>'glx',
        'X'=>'Xaa'
    );

=head1 SYNOPSIS

We produces each codon frequency,
	    Relative Synonymous Codons Uses and
	    Relative Adaptiveness of a Codon

that will help you to calculate the Codon Adaptation Index (CAI) of a gene, to see the gene expression level.

Perhaps a little code snippet.

    use Bio::Tools::CodonOptTable;

    my $seqobj = Bio::Tools::CodonOptTable->new ( -seq => 'ATGGGGTGGGCACCATGCTGCTGTCGTGAATTTGGGCACGATGGTGTACGTGCTCGTAGCTAGGGTGGGTGGTTTG',
				   -id  => 'GeneFragment-12',
				   -accession_number => 'Myseq1',
				   -alphabet => 'dna',
				   -is_circular => 1
				   );

    #If you wanna read from file
    my $seqobj = Bio::Tools::CodonOptTable->new(-file => "contig.fasta",
                                             -format => 'Fasta');
    
    #If you have Accession number and want to get file from NCBI
    my $seqobj = Bio::Tools::CodonOptTable->new(-ncbi_id => "J00522");
    
    my $myCodons = $seqobj->rscu_rac_table();
    
    if($myCodons)
	{
	    for my $a (@$myCodons)
		    {
		       print "Codon 	: ",$a->{'codon'},"\t";
		       print "Frequency : ",$a->{'frequency'},"\t";
		       print "AminoAcid : ",$a->{'aa_name'},"\t";
		       print "RSCU	: ",$a->{'rscu'},"\t"; #Relative Synonymous Codons Uses
		       print "RAC	: ",$a->{'rac'},"\t"; #Relative Adaptiveness of a Codon
		       print "\n";
		    }
	}
    ...
    
=head1 METHODS
    Title   : new
    Usage1   : $seq    = Bio::Tools::CodonOptTable->new( -seq => 'ATGGGGGTGGTGGTACCCT',
					      -id  => 'human_id',
					      -accession_number => 'AL000012',
					      );
					      
    Usage2   : $seq    = Bio::Tools::CodonOptTable->new( -file => 'myseq.fasta',
					      -format => 'fasta',
					      );
					      
    Usage3   : $seq    = Bio::Tools::CodonOptTable->new( -ncbi_id => 'J00522');
					      
    Function: Returns a new primary seq object from
	      basic constructors, being a string for the sequence
	      and strings for id and accession_number.
	    
    Returns : a new Bio::PrimarySeq object
    
    Args    : -seq         	=> sequence string
	      -display_id  	=> display id of the sequence (locus name) 
	      -accession_number => accession number
	      -primary_id  	=> primary id (Genbank id)
	      -desc        	=> description text
	      -alphabet    	=> molecule type (dna,rna,protein)
	      -id          	=> alias for display id
	      -file		=> file location
	      -format		=> file format
	      -ncbi_id		=> NCBI accession number
	      
    Note    : IF you are reading sequence from file it will call _read_localfile method
	      IF you are fetching file form NCBI it will call _read_remotefile method

=head2 METHODS
    Title   : calculate_rscu
					      
    Function: Calculate the RSCU(Relative Synonymous Codons Uses).
	    	      
    Note    : The formula is used in the following references.
	    http://www.pubmedcentral.nih.gov/articlerender.fcgi?tool=pubmed&pubmedid=3547335

=head2 METHODS
    Title   : calculate_rac
					      
    Function: Calculate the RAC(Relative Adaptiveness of a Codon).
	    	      
    Note    : The formula is used in the following references.
	    http://www.pubmedcentral.nih.gov/articlerender.fcgi?tool=pubmed&pubmedid=3547335

=cut

sub new {
    my($class,@args) = @_;
    my $self = $class->SUPER::new(@args);
    my $seqobj;
    my($seq,$id,$acc,$pid,$desc,$alphabet,$given_id,$is_circular,$file,$format,$ncbi_id) =
	$self->_rearrange([qw(SEQ
			      DISPLAY_ID
			      ACCESSION_NUMBER
			      PRIMARY_ID
			      DESC
			      ALPHABET
			      ID
			      IS_CIRCULAR
			      FILE
			      FORMAT
                              NCBI_ID
			      )],
			  @args);
    if($file && $format)
    {
	($seq,$id,$alphabet) = _read_localfile($file,$format);
    }
    if($ncbi_id)
    {
        ($seq,$id,$desc,$alphabet) = _read_remotefile($ncbi_id);   
    }
    $seqobj = Bio::PrimarySeq->new ( -seq => $seq,
                                     -id  => $id,
                                     -accession_number => $acc,
                                     -display_id => $given_id,
                                     -desc => $desc,
                                     -primary_id => $pid,
                                     -alphabet => $alphabet,
                                     -is_circular => $is_circular
                                   );
    return bless $seqobj, $class;
}

sub _read_remotefile
{
    my $ncbi_id = $_[0];
    
    my($seq,$id,$desc,$alphabet);
    
    my $retrivefile = new Bio::DB::GenBank(-retrievaltype => 'tempfile' , 
                                           -format => 'Fasta');;
    
    my $fetchedfile = $retrivefile->get_Stream_by_acc($ncbi_id);
    my $seq_data    =  $fetchedfile->next_seq;
    
    $seq        = $seq_data->seq;
    $id         = $seq_data->primary_id;
    $desc       = $seq_data->desc;
    $alphabet   = $seq_data->alphabet;
    
    return($seq,$id,$desc,$alphabet);
}

sub _read_localfile
{
    my ($file,$format) = @_;
    
    my($seq,$id,$alphabet);
    
    my $inputstream = Bio::SeqIO->new(-file => $file,-format => $format);
    my $input       = $inputstream->next_seq();
    
    $seq        = $input->seq();
    $id         = $input->id;
    $alphabet   = $input->alphabet;
    
    return($seq,$id,$alphabet);
}

sub rscu_rac_table
{
    my $seqobj = $_[0];
    
    my $seq_stats  =  Bio::Tools::SeqStats->new(-seq=>$seqobj);
    my $codons = $seq_stats-> count_codons();
    
    my $rscu_rac = map_codon_iupac($codons);
    
    return $rscu_rac;
}

sub map_codon_iupac
{
    my $codons = $_[0];
    
    my $myCodonTable   = Bio::Tools::CodonTable->new();
    
    my (@myCodons, @maxfreq_in_aa);
    my %frequency_table_aa;

    foreach my $single_codon (keys %$codons)
    {
        my $aa_name_abri = $myCodonTable->translate($single_codon); 
        my $aa_name = $Amnioacid{$aa_name_abri};
	push @myCodons, {
		'codon' => $single_codon,
		'frequency'   => $codons->{$single_codon},
		'aa_name'    => $aa_name,
	};
	if(!defined($frequency_table_aa{$aa_name}) || ($frequency_table_aa{$aa_name} < $codons->{$single_codon}))
	{
	    $frequency_table_aa{$aa_name} = $codons->{$single_codon};
	}
    }
    &calculate_rscu(\@myCodons,\%frequency_table_aa);
}

sub calculate_rscu
{
    my($codons,$max_codons) = @_;
    
    my($freq,$rscu,$rac,@myCodons);
    my %rscu_max_table;
    
    foreach my $each_codon (@$codons)
    {
	my $amino = $each_codon->{'aa_name'};
	my $freq  = $each_codon->{'frequency'};
	my $count = 0;
	my $all_freq_aa = 0;
	if($amino)
	{
	    foreach my $goforall (@$codons)
	    {
		if(defined $amino && $goforall->{'aa_name'} eq $amino)
		{
		    $all_freq_aa+= $goforall->{'frequency'};
		    $count++;
		}
	    }	    
	    $rscu = $count*$freq/$all_freq_aa;
	}
	if(!defined($rscu_max_table{$amino}) || ($rscu_max_table{$amino} < $rscu))
	{
	    $rscu_max_table{$amino} = $rscu;
	}
	push @myCodons, {
		'codon' => $each_codon->{'codon'},
		'frequency'   => $freq,
		'aa_name'    => $amino,
		'rscu'	=> $rscu,
		'total_aa_comb' => $count,
		'all_fre_aa' => $all_freq_aa,
	};
    }
    &calculate_rac(\@myCodons,\%rscu_max_table)
}

sub calculate_rac
{
    my($codons,$max_rscu) = @_;
    my($freq,$rac,@myCodons);
    
    foreach my $each_codon (@$codons)
    {
	my $amino = $each_codon->{'aa_name'};
	my $rscu  = $each_codon->{'rscu'};
	if($amino)
	{
	    my $max = $max_rscu->{$amino};
	    $rac = $rscu/$max;
	    push @myCodons, {
		'codon' 	=> $each_codon->{'codon'},
		'frequency'   	=> $each_codon->{'frequency'},
		'aa_name'    	=> $amino,
		'rscu'		=> $rscu,
		'rac' 		=> $rac,
	    };
	}
    }
    # A CAI of 1.0 is considered to be perfect in the desired expression organism, and a
    # CAI of >0.9 is regarded as very good, in terms of high gene expression level.
    return \@myCodons;
}

=head1 AUTHOR

Rakesh Kumar Shardiwal, C<< <rakesh.shardiwal at gmail.com> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-bio-tools-codonopttable at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Bio-Tools-CodonOptTable>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.


=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Bio::Tools::CodonOptTable


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=Bio-Tools-CodonOptTable>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/Bio-Tools-CodonOptTable>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/Bio-Tools-CodonOptTable>

=item * Search CPAN

L<http://search.cpan.org/dist/Bio-Tools-CodonOptTable/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 COPYRIGHT & LICENSE

Copyright 2008 Rakesh Kumar Shardiwal, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.


=cut

1; # End of Bio::Tools::CodonOptTable
