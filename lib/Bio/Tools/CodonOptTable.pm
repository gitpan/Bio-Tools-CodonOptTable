package Bio::Tools::CodonOptTable;

use strict;
use Data::Dumper;

use Bio::Root::RootI;
use Bio::Root::Root;
use Bio::PrimarySeq;
use Bio::SeqIO;
use Bio::Tools::SeqStats;
use Bio::Tools::CodonTable;
use Bio::DB::GenBank;
use GD::Graph::bars;

=head1 NAME

Bio::Tools::CodonOptTable - A more elaborative way to check the codons usage.

=head1 VERSION

Version 0.06

=cut

our $VERSION = '0.06';

use vars qw(@ISA %Amnioacid);

@ISA = ('Bio::Root::Root', 'Bio::SeqIO', 'Bio::PrimarySeq');
use base qw(Exporter);

our @EXPORT = qw(
        		    new
                    rscu_rac_table
                    prefered_codon
                    generate_graph
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


=head1 DESCRIPTION

The purpose of this module is to show codon usage.

We produces each codon frequency,
	    Relative Synonymous Codons Uses and
	    Relative Adaptiveness of a Codon table and bar graph
that will help you to calculate the Codon Adaptation Index (CAI) of a gene, to see the gene expression level.

Relative Synonymous Codons Uses(RSCU) values are the number of times a particular codon is observed, relative to the number of times
that the codon would be observed in the absence of any codon usage bias.

In the absence of any codon usage bias, the RSCU value would be 1.00.
A codon that is used less frequently than expected will have a value of less than 1.00 and vice versa for a codon that is used more frequently than expected.

Genetics Code: NCBI takes great care to ensure that the translation for each coding sequence (CDS) present in GenBank records is correct. Central to this effort is careful checking on the taxonomy of each record and assignment of the correct genetic code (shown as a /transl_table qualifier on the CDS in the flat files) for each organism and record. This page summarizes and references this work.
http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi

=cut

=head1 SYNOPSIS

    use Bio::Tools::CodonOptTable;

    my $seqobj = Bio::Tools::CodonOptTable->new ( -seq => 'ATGGGGTGGGCACCATGCTGCTGTCGTGAATTTGGGCACGATGGTGTACGTGCTCGTAGCTAGGGTGGGTGGTTTG',
                                                -id  => 'GeneFragment-12',
                                                -accession_number => 'Myseq1',
                                                -alphabet => 'dna',
                                                -is_circular => 1,
                                                -genetic_code => 1,
				   );

    B<#If you wanna read from file>
    my $seqobj = Bio::Tools::CodonOptTable->new(-file => "contig.fasta",
                                             -format => 'Fasta',
                                             -genetic_code => 1,
                                             );

    B<#If you have Accession number and want to get file from NCBI>
    my $seqobj = Bio::Tools::CodonOptTable->new(-ncbi_id => "J00522",
                                                -genetic_code => 1,);

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
    
    B<# To get the prefered codon list based on RSCU & RAC Values >
    my $prefered_codons = $seqobj->prefered_codon($myCodons);

    while ( my ($amino_acid, $codon) = each(%$prefered_codons) ) {
        print "AminoAcid : $amino_acid \t Codon : $codon\n";
    }
    
    B<# To produce a graph between RSCU & RAC>
    # Graph output file extension should be GIF, we support GIF only
    
    $seqobj->generate_graph($myCodons,"myoutput.gif");
    
    ...
=cut

=head1 METHODS

=head2 Constructor

    Title   : new
    
    Usage1   : $seq    = Bio::Tools::CodonOptTable->new( -seq => 'ATGGGGGTGGTGGTACCCT',
                                                        -id  => 'human_id',
                                                        -accession_number => 'AL000012',
                                                        -genetic_code => 1,
					      );
					      
    Usage2   : $seq    = Bio::Tools::CodonOptTable->new( -file => 'myseq.fasta',
                                                        -format => 'fasta',
                                                        -genetic_code => 1,
					      );
					      
    Usage3   : $seq    = Bio::Tools::CodonOptTable->new( -ncbi_id => 'J00522',
                                                        -genetic_code => 1,
                                                    );
					      
    Function: Returns a new primary seq object from
	      basic constructors, being a string for the sequence
	      and strings for id and accession_number.
	    
    Returns : a new Bio::PrimarySeq object
    
    Args    :   -seq         	    => sequence string
                -display_id  	    => display id of the sequence (locus name) 
                -accession_number   => accession number
                -primary_id  	    => primary id (Genbank id)
                -desc        	    => description text
                -alphabet    	    => molecule type (dna,rna,protein)
                -id          	    => alias for display id
                -file               => file location
                -format             => file format
                -ncbi_id            => NCBI accession number
                -genetic_code       => 1 (Default)
                                       # See this list to know more about genetic_code
                                       # http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
	      
    Note    : IF you are reading sequence from file it will call _read_localfile method
              IF you are fetching file form NCBI it will call _read_remotefile method

=head2 Calculate RSCU
    
    Title   : calculate_rscu
					      
    Function: Calculate the RSCU(Relative Synonymous Codons Uses).
	    	      
    Note    : The formula is used in the following references.
	    http://www.pubmedcentral.nih.gov/articlerender.fcgi?tool=pubmed&pubmedid=3547335

=head2 Calculate RAC
    
    Title   : calculate_rac
					      
    Function: Calculate the RAC(Relative Adaptiveness of a Codon).
	    	      
    Note    : The formula is used in the following references.
	    http://www.pubmedcentral.nih.gov/articlerender.fcgi?tool=pubmed&pubmedid=3547335

=head2 Get Prefered Codons based on RAC & RSCU
    
    Title   : prefered_codon
					      
    Function: Give you prefered codons list.

=head2 Produce RSCU & RAC Graph
    
    Title   : generate_graph
					      
    Function: Produce a bar graph between RAC(Relative Adaptiveness of a Codon) & RSCU(Relative Synonymous Codons Uses).

=cut


=pod
    Constructor
=cut
sub new {
    my($class,@args) = @_;
    
    my $self = $class->SUPER::new(@args);
    my $seqobj;
    
    my($seq,$id,$acc,$pid,$desc,$alphabet,$given_id,$is_circular,$file,$format,$ncbi_id,$genetic_code) =
	$self->_rearrange([qw(
                            SEQ
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
                            GENETIC_CODE
                    )],@args);
    
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

=pod
    Function to read remote file
=cut
sub _read_remotefile
{
    my $ncbi_id = $_[0];
    
    my($seq,$id,$desc,$alphabet);
    
    my $retrivefile = new Bio::DB::GenBank(-retrievaltype => 'tempfile' , 
                                           -format => 'Fasta');
    
    my $fetchedfile = $retrivefile->get_Stream_by_acc($ncbi_id);
    my $seq_data    =  $fetchedfile->next_seq;
    
    $seq        = $seq_data->seq;
    $id         = $seq_data->primary_id;
    $desc       = $seq_data->desc;
    $alphabet   = $seq_data->alphabet;
    
    return($seq,$id,$desc,$alphabet);
}
=pod
    Function to read local file
=cut
sub _read_localfile
{
    my ($file,$format) = @_;
    
    my($seq,$id,$alphabet);
    
    my $inputstream = Bio::SeqIO->new(-file => $file,
                                      -format => $format);
    my $input       = $inputstream->next_seq();
    
    $seq        = $input->seq();
    $id         = $input->id;
    $alphabet   = $input->alphabet;
    
    return($seq,$id,$alphabet);
}
=pod
    Function to produce rscu and rac table
=cut
sub rscu_rac_table
{
    my ($seqobj) = @_;
    
    my $seq_stats  =  Bio::Tools::SeqStats->new(-seq=>$seqobj);
    my $codons = $seq_stats-> count_codons();
    
    my $rscu_rac = map_codon_iupac($codons);
    
    my @sorted_codons_by_aa = sort { $a->{aa_name} cmp $b->{aa_name} } @$rscu_rac;
    
    return \@sorted_codons_by_aa;
}
=pod
    Function to map codon and iupac names
=cut
sub map_codon_iupac
{
    my ($codons,$genetics_code) = @_;
    
    my $myCodonTable   = Bio::Tools::CodonTable->new();
    
    # change codon table if require
    if(!$genetics_code) { $genetics_code = 1; }
    $myCodonTable->id($genetics_code);
    
    my (@myCodons, @maxfreq_in_aa);
    my %frequency_table_aa;

    foreach my $single_codon (keys %$codons)
    {
        my $aa_name_abri    = $myCodonTable->translate($single_codon); 
        my $aa_name         = $Amnioacid{$aa_name_abri};
        push @myCodons, {
            'codon'         => $single_codon,
            'frequency'     => $codons->{$single_codon},
            'aa_name'       => $aa_name,
        };
        if(!defined($frequency_table_aa{$aa_name}) ||
          ($frequency_table_aa{$aa_name} < $codons->{$single_codon})
          )
        {
            $frequency_table_aa{$aa_name} = $codons->{$single_codon};
        }
    }
    calculate_rscu(\@myCodons,\%frequency_table_aa);
}
=pod
    Function to calculate rscu
=cut
sub calculate_rscu
{
    my($codons,$max_codons) = @_;
    
    my($freq,$rscu,$rac,@myCodons);
    my %rscu_max_table;
    
    foreach my $each_codon (@$codons)
    {
        my $amino       = $each_codon->{'aa_name'};
        my $freq        = $each_codon->{'frequency'};
        my $count       = 0;
        my $all_freq_aa = 0;
        if($amino)
        {
            foreach my $goforall (@$codons)
            {
                if( $amino && ($goforall->{'aa_name'} eq $amino))
                {
                    $all_freq_aa += $goforall->{'frequency'};
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
            'codon'         => $each_codon->{'codon'},
            'frequency'     => $freq,
            'aa_name'       => $amino,
            'rscu'	        => $rscu,
            'total_aa_comb' => $count,
            'all_fre_aa'    => $all_freq_aa,
        };
    }
    calculate_rac(\@myCodons,\%rscu_max_table)
}
=pod
    Function to calculate rac value
=cut
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
            'frequency' => $each_codon->{'frequency'},
            'aa_name'   => $amino,
            'rscu'		=> sprintf("%.5f", $rscu),
            'rac' 		=> sprintf("%.5f", $rac),
            };
        }
    }
    # A CAI of 1.0 is considered to be perfect in the desired expression organism, and a
    # CAI of >0.9 is regarded as very good, in terms of high gene expression level.
    return (\@myCodons);
}
=pod
    Function to get prefered codons
=cut
sub prefered_codon
{
    my($self,$codons) = @_;
    my $prefered_codon;
    for my $each_aa (@$codons)
	{
		my $aa_name 	= $each_aa->{'aa_name'};
		my $rscu 		= $each_aa->{'rscu'} ;
		my $codon 		= $each_aa->{'codon'};
		my $frequency 	= $each_aa->{'frequency'};
		
		if(!defined($prefered_codon->{$aa_name}) || ($prefered_codon->{$aa_name} < $rscu))
		{
			$prefered_codon->{$aa_name} = $codon;
		}
	}
    return $prefered_codon;
}
=pod
    Function to generate graph
=cut
sub generate_graph
{
    my($self,$codons,$output_file) =@_;

    my (@x_axis_labels,@rscu,@rac,@x_axis_values,@codons_table,@codon_freq);
    my $y_axis_max 		    = 5;
    my @category_colours 	= qw(red dgreen);
    my $bar_graph 		    = new GD::Graph::bars(1000,500);
    
    foreach my $each_aa (@$codons) 
	{
	    if($each_aa->{'aa_name'})
	    {
		push(@codons_table,$each_aa->{'aa_name'}."(".$each_aa->{'codon'}.")");
		push(@codon_freq,$each_aa->{'frequency'});
		push(@x_axis_labels,$each_aa->{'codon'}."(".$each_aa->{'frequency'}.")"."-".$each_aa->{'aa_name'});
		push(@rscu,$each_aa->{'rscu'});
		push(@rac,$each_aa->{'rac'});
	    }
	}
    
    my @bar_graph_table;
    push(@bar_graph_table, \@x_axis_labels);
    push(@bar_graph_table, \@rscu);
    push(@bar_graph_table, \@rac);
    
    $bar_graph->set(
        title               => 'Graph Representing : Relative Synonymous Codons Uses and Relative Adaptiveness of a Codon for '.$self->display_id,
        y_label             => 'RSCU and RAC values',   #y-axis label
        y_max_value         => $y_axis_max,             #the max value of the y-axis
        y_min_value         => 0,                       #the min value of y-axis, note set below 0 if negative values are required
        y_tick_number       => 20,                      #y-axis scale increment
        y_label_skip        => 1,                       #label every other y-axis marker
        box_axis            => 0,                       #do not draw border around graph
        line_width          => 2,                       #width of lines
        legend_spacing      => 5,                       #spacing between legend elements
        legend_placement    =>'RC',                     #put legend to the centre right of chart
        dclrs               => \@category_colours,      #reference to array of category colours for each line
        bgclr               => 'red',
        long_ticks          => 0,
        tick_length         => 3,
        x_labels_vertical   => 1,
    ) || die "\nFailed to create line graph: $bar_graph->error()";
        
    $bar_graph->set_legend('RSCU Value', 'RAC Value');
    
    my $plot 	= $bar_graph->plot(\@bar_graph_table);
    
    my $line_file = $output_file;
    open(GPH, ">$line_file") || die ("\nFailed to save graph to file: $line_file. $!");
    binmode(GPH);
    print GPH $plot->gif();
    close(GPH);
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

Lalchand Kumawat <lalchand82@gmail.com> 
Rajneesh Kumar Sharma <biorajneesh@gmail.com>

=head1 COPYRIGHT & LICENSE

Copyright 2008 Rakesh Kumar Shardiwal, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.


=cut

1; # End of Bio::Tools::CodonOptTable

__END__