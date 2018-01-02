=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::DB::ExternalFeatureFactoryI -
Legacy Abstract interface for External Feature
Factories. Bio::EnsEMBL::External::ExternalFeatureAdaptor should be used
instead if possible.

=head1 SYNOPSIS

  $external_ff = new ImplementingExternalFeatureFactoryClass;

  $database_adaptor = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host   => 'blah',
    -dbname => 'other',
    -pass   => 'pass'
  );

  # alternatively, you can add external databases to an obj once made
  $database_adaptor->add_ExternalFeatureFactory($external_ff);

  # now the ExternalFeatureFactory has been added, Ensembl RawContigs
  # and Slices will now have ExternalFeatures on them
  $contig =
    $db_adaptor->get_RawContigAdaptor->fetch_by_name('AC00056.00001');
  @external = @{ $contig->get_all_ExternalFeatures() };

  # this works on Slices as well
  $slice =
    $db_adaptor->get_SliceAdaptor->fetch_by_chr_start_end( '12', 10000,
    30000 );
  @external = @{ $slice->get_all_ExternalFeatures() };

=head1 DESCRIPTION

This is a legacy class.  It is included only for backwards
compatibility with ExternalFeatureFactories which are presumably
still used to place data into ensembl.  It is recommended that if
you wish to create EnsEMBL features externally that you use the
Bio::EnsEMBL::External::ExternalFeatureAdaptor instead.

This object defines the abstract interface for External Database access
inside Ensembl. The aim is that one can attach an External Database
which will generate Sequence Features and these Sequence Features will
be accessible along side all the internal Ensembl sequence features, for
drawing, EMBL dumping etc. In particular, the external database does not
have to worry about the transformation of the Sequence Feature objects
into VirtualContigs.

Sequence Features have to be defined in one of two coordinate systems:
Original EMBL/GenBank coordinates of a particular sequnence version or
the Ensembl contig coordinates. This means you have to calculate your
sequence features in one these two coordinate systems

The methods that have to be implemented are:

  get_External_SeqFeatures_contig( $ensembl_contig_identifier,
    $sequence_version, $start, $end );

  get_External_SeqFeatures_clone( $embl_accession_number,
    $sequence_version, $start, $end );

The semantics of this method is as follows:

  $ensembl_contig_identifier - the ensembl contig id (external id).
  $sequence_version - embl/genbank sequence version
  $embl_accession_number - the embl/genbank accession number

The $start/$end can be ignored, but methods can take advantage of it.
This is so that ensembl can ask for features only on a region of DNA,
and if desired, the external database can respond with features only in
this region, rather than the entire sequence.

The hope is that the second method could potentially have a very complex
set of mappings of other embl_accession numbers to one embl_accession
number and provide the complex mapping.

The methods should return Sequence Features with the following spec:

  a) must implement the Bio::SeqFeatureI interface.

  b) must accept "set" calls on 

  start,end,strand

  to provide coordinate transformation of the feature.

  c) must be unique in-memory objects, ie, the implementation is not
  allowed to cache the sequence feature in its entirity. Two separate
  calls to get_External_SeqFeatures_contig must be able to separately
  set start,end,strand information without clobbering each other. The
  other information, if so wished, can be cached by each SeqFeature
  holding onto another object, but this is left to the implementor to
  decide on the correct strategy.

  d) must return an unique identifier when called with method id.

You must implement both functions. In most cases, one function will
always return an empty list, whereas the other function will actually
query the external database.

The second way of accessing the External Database from Ensembl is using
unique internal identifiers in that database. The method is:

  get_SeqFeature_by_id($id);

It should return exactly one Sequence Feature object of the same type as
above.

=head1 METHODS

=cut

package Bio::EnsEMBL::DB::ExternalFeatureFactoryI;

use strict;
use warnings;

use Bio::EnsEMBL::External::ExternalFeatureAdaptor;
use vars qw(@ISA);

@ISA = ( 'Bio::EnsEMBL::External::ExternalFeatureAdaptor' );


=head2 coordinate_systems

  Arg [1]    : none
  Example    : none
  Description: This method is present to make the ExternalFeatureFactory 
               interface behave as an ExternalFeatureAdaptor. It is for
               backwards compatibility.
  Returntype : none
  Exceptions : none
  Caller     : internal

=cut

sub coordinate_systems {
  my $self = shift;
  return qw(CONTIG);
}


=head2 fetch_all_by_contig_name

  Arg [1]    : none
  Example    : none
  Description: This method is present to make the ExternalFeatureFactory 
               interface behave as an ExternalFeatureAdaptor. It is for
               backwards compatibility.
  Returntype : none
  Exceptions : none
  Caller     : internal

=cut

sub fetch_all_by_contig_name {
   my ($self, $contig_name) = @_;

   unless($self->db) {
     $self->throw('DB attribute not set.  This value must be set for the ' .
		  'ExternalFeatureFactory to function correctly');
   }

   my @features = ();

   my $ctg = $self->db->get_RawContigAdaptor->fetch_by_name($contig_name);
   my $clone = $ctg->clone;
   my $version = $clone->version;
   my $ctg_length = $ctg->length;
   
   #get contig features
   push @features, $self->get_Ensembl_SeqFeatures_contig($ctg->name, 
							 $version, 
							 1,
							 $ctg_length);

   #get clone features
   my $clone_start = $ctg->embl_offset;
   my $clone_end   = $clone_start + $ctg_length - 1;
   my @clone_features = $self->get_Ensembl_SeqFeatures_clone($clone->id,
							     $version,
							     $clone_start,
							     $clone_end);

   #change clone coordinates to contig coordinates
   my ($start, $end); 
   foreach my $f (@clone_features) {
     $start = $f->start - $clone_start + 1;
     $end   = $f->end   - $clone_start + 1;

     #skip features outside the contig
     next if($end < 1 || $start > $ctg_length);

     $f->start($start);
     $f->end($end);

     push @features, $f;
   }

   return \@features;
}

=head2 get_Ensembl_SeqFeatures_contig

 Title   : get_Ensembl_SeqFeatures_contig (Abstract)
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Ensembl_SeqFeatures_contig{
   my ($self) = @_;

   $self->warn("Abstract method get_External_SeqFeatures_contig " .
	    "encountered in base class. Implementation failed to complete it");

}

=head2 get_Ensembl_SeqFeatures_clone

 Title   : get_Ensembl_SeqFeatures_clone (Abstract)
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Ensembl_SeqFeatures_clone{
   my ($self) = @_;
   
   $self->warn("Abstract method get_Ensembl_SeqFeatures_clone " .
	    "encountered in base class. Implementation failed to complete it");

}

=head2 get_Ensembl_Genes_clone

 Title   : get_Ensembl_Genes_clone
 Function: returns Gene objects in clone coordinates from a gene id
 Returns : An array of Gene objects
 Args    : clone id

=cut

sub get_Ensembl_Genes_clone {
    my $self = @_;

    return;
}

=head2 get_SeqFeature_by_id

 Title   : get_SeqFeature_by_id (Abstract)
 Usage   : 
 Function: Return SeqFeature object for any valid unique id  
 Example :
 Returns : 
 Args    : id as determined by the External Database


=cut

       
sub get_SeqFeature_by_id {
   my ($self) = @_;
   $self->warn("Abstract method get_SeqFeature_by_id  encountered " .
	       "in base class. Implementation failed to complete it");
}


1;







