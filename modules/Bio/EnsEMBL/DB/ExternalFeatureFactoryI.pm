
#
# BioPerl module for Bio::EnsEMBL::DB::ExternalFeatureFactoryI
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DB::ExternalFeatureFactoryI - Abstract interface for
databases providing external features

=head1 SYNOPSIS

   # you have to make an Object which conforms to this interface
   $mydb = SomeDatabase::SomeInterface->new();

   # you can make an obj with external databases
   $dbobj = Bio::EnsEMBL::DBSQL::Obj->new( -host => 'blah',
					   -dbname => 'other',
					   # other arguments
					   -external => [ $mydb ] );

   # alternatively, you can add external databases to an obj once made
   $dbobj->add_ExternalFeatureFactory($mydb);

   # now the ExternalFeatureFactory has been added, Ensembl RawContigs
   # and virtualcontigs will now have ExternalFeatures() on them

   # could get a contig in any manner
   $contig = $dbobj->get_Contig('AC00056.00001');
   @external = $contig->get_all_ExternalFeatures();

   # this works on VirtualContigs as well

   $vc = Bio::EnsEMBL::DB::VirtualContig ( -focuscontig => $contig,
					   -focusposition => 10,
					   -ori => 1,
					   -left => 1000,
					   -right => 1000
					   );


   @external = $vc->get_all_ExternalFeatures();

   

=head1 DESCRIPTION

This object defines the abstract interface for External Database
access inside Ensembl. The aim is that one can attach an External
Database which will generate Sequence Features and these Sequence
Features will be accessible along side all the internal Ensembl
sequence features, for drawing, EMBL dumping etc. In particular, the 
external database does not have to worry about the transformation of
the Sequence Feature objects into VirtualContigs.

Sequence Features have to be defined in one of two coordinate systems:
Original EMBL/GenBank coordinates of a particular sequnence version or the
Ensembl contig coordinates. This means you have to calculate your sequence
features in one these two coordinate systems

The methods that have to be implemented are:

    get_External_SeqFeatures_contig(
        $ensembl_contig_identifier,$sequence_version,$start,$end);

    get_External_SeqFeatures_clone(
        $embl_accession_number,$sequence_version,$start,$end);

The semantics of this method is as follows:

    $ensembl_contig_identifier - the ensembl contig id (external id).
    $sequence_version - embl/genbank sequence version
    $embl_accession_number - the embl/genbank accession number

The $start/$end can be ignored, but methods can take advantage of it.
This is so that ensembl can ask for features only on a region of DNA,
and if desired, the external database can respond with features only
in this region, rather than the entire sequence.

The hope is that the second method could potentially have a very
complex set of mappings of other embl_accession numbers to one
embl_accession number and provide the complex mapping.

The methods should return Sequence Features with
the following spec:

    a) must implement the Bio::SeqFeatureI interface.

    b) must accept "set" calls on 

    start,end,strand

    to provide coordinate transformation of the feature.

    c) must be unique in-memory objects, ie, the implementation is not
    allowed to cache the sequence feature in its entirity. Two
    separate calls to get_External_SeqFeatures_contig must be able to
    separately set start,end,strand information without clobbering
    each other. The other information, if so wished, can be cached by
    each SeqFeature holding onto another object, but this is left to
    the implementor to decide on the correct strategy.

    d) must return an unique identifier when called with method id. 

You must implement both functions. In most cases, one function will
always return an empty list, whereas the other function will actually
query the external database.

The second way of accessing the External Database from Ensembl is
using unique internal identifiers in that database. The method is:

    get_SeqFeature_by_id($id);

It should return exactly one Sequence Feature object of the same type
as above.


=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  vsns-bcd-perl@lists.uni-bielefeld.de          - General discussion
  vsns-bcd-perl-guts@lists.uni-bielefeld.de     - Technically-oriented discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Ewan Birney

Email birney@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DB::ExternalFeatureFactoryI;
use Bio::Root::RootI;
use vars qw(@ISA);

@ISA = ( 'Bio::Root::RootI');

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

   $self->warn("Abstract method get_External_SeqFeatures_contig encountered in base class. Implementation failed to complete it");

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
   
   $self->warn("Abstract method get_Ensembl_SeqFeatures_clone encountered in base class. Implementation failed to complete it");

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
   $self->warn("Abstract method get_SeqFeature_by_id  encountered in base class. Implementation failed to complete it");


}

1;







