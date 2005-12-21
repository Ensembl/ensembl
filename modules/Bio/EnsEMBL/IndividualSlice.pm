#
# Ensembl module for Bio::EnsEMBL::IndividualSlice
#
#
# Copyright Team Ensembl
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::IndividualSlice - SubClass of the Slice. Represents the slice of the genome 
for a certain individual (applying the alleles for this individual)

=head1 SYNOPSIS

   $sa = $db->get_SliceAdaptor;

   $slice = $sa->fetch_by_region('chromosome', 'X', 1_000_000, 2_000_000);

   $individualSlice = $slice->get_by_Individual($individual_name);

   #get the sequence from the Individual Slice: will contain IUPAC codes for SNPs and Ensembl ambiguity codes for indels
   my $seq = $individualSlice->seq();
   print $seq;

   #get a subSlice of the Strain
   my $subSlice_individual = $individualSlice->sub_Slice(5_000,8_000,1)

   #compare two different individuals in the same Slice
   my $sliceIndividual2 = $slice->get_by_Individual($individual_name2);
   my $differences = $individualSlice->get_all_differences_IndividualSlice($sliceIndividual2);
   foreach my $af (@{$differences}){
      print "There is a difference between $individual_name and $individual_name2 at ", $af->start,"-",$af->end,
            " with allele ", $af->allele_string(),"\n";
   }

=head1 DESCRIPTION

A IndividualSlice object represents a region of a genome for a certain individual.  It can be used to retrieve
sequence or features from a individual.

=head1 CONTACT

This modules is part of the Ensembl project http://www.ensembl.org

Questions can be posted to the ensembl-dev mailing list:
ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

package Bio::EnsEMBL::IndividualSlice;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Mapper;
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate warning);

use Data::Dumper;
@ISA = qw(Bio::EnsEMBL::Slice);


=head2 new

    Arg [1..N]  : List of named arguments
                  Bio::EnsEMBL::CoordSystem COORD_SYSTEM
                  string SEQ_REGION_NAME,
                  int    START,
                  int    END,
                  string VERSION (optional, defaults to '')
                  int    STRAND, (optional, defaults to 1)
                  Bio::EnsEMBL::DBSQL::SliceAdaptor ADAPTOR (optional)
    Arg[N+1]      : string $individual_name
    Example     : $individualSlice = Bio::EnsEMBL::IndividualSlice->new(-coord_system => $cs,
									-start => 1,
									-end => 10000,
									-strand => 1,
									-seq_region_name => 'X',
									-seq_region_length => 12e6,
									-individual_name => $individual_name);
    Description : Creates a new Bio::EnsEMBL::IndividualSlice object that will contain a shallow copy of the
                  Slice object, plus additional information such as the individual this Slice refers to
                  and listref of Bio::EnsEMBL::Variation::AlleleFeatures of differences with the
                  reference sequence
    ReturnType  : Bio::EnsEMBL::IndividualSlice
    Exceptions  : none
    Caller      : general

=cut

sub new{
    my $caller = shift;
    my $class = ref($caller) || $caller;

    #create the IndividualSlice object as the Slice, plus the individual attribute
    my ($individual_name) = rearrange(['INDIVIDUAL'],@_);

    my $self = $class->SUPER::new(@_);

    $self->{'individual_name'} = $individual_name;

    if(!$self->adaptor()) {
	warning('Cannot get new IndividualSlice features without attached adaptor');
	return '';
    }
    my $variation_db = $self->adaptor->db->get_db_adaptor('variation');

    unless($variation_db) {
	warning("Variation database must be attached to core database to " .
		"retrieve variation information" );
	return '';
    }
    #get the AlleleFeatures in the Individual
    my $af_adaptor = $variation_db->get_AlleleFeatureAdaptor;
    
    if( $af_adaptor ) {
	#set the adaptor to retrieve data from genotype table
	$af_adaptor->from_IndividualSlice(1);
	#get the Individual for the given strain
	my $ind_adaptor = $variation_db->get_IndividualAdaptor;

	if ($ind_adaptor){
	    my $individual = $ind_adaptor->fetch_all_by_name($self->{'individual_name'}); #ignore individuals with same name
	    #check that there is such individual in the database
	    if (defined $individual->[0]){
		#get all the VariationFeatures in the $individual and the Slice given
		my $allele_features = $af_adaptor->fetch_all_by_Slice_Individual($self,$individual->[0]);
		
		$self->{'alleleFeatures'} = $allele_features;
		return $self;
	    }
	    else{ 
		warning("Strain not in the database");
		return '';
	    }
	}
	else{
	    warning("Not possible to retrieve IndividualAdaptor from the variation database");
	    return '';
	}
    } else {
	warning("Not possible to retrieve AlleleFeatureAdaptor from variation database");
	return '';
    }
}

=head2 individual_name

    Arg [1]     : (optional) string $individual_name
    Example     : my $individual_name = $individualSlice->individual_name();
    Description : Getter/Setter for the name of the individual in the slice
    ReturnType  : string
    Exceptions  : none
    Caller      : general

=cut

sub Individual{
   my $self = shift;
   if (@_){
       $self->{'individual_name'} = shift @_;
   }
   return $self->{'individual_name'};
}
