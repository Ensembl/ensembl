
#
# BioPerl module for ProteinAdaptor
#
# Cared for by Emmanuel Mongin <mongin@ebi.ac.uk>
#
# Copyright Emmanuel Mongin
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

ProteinAdaptor - DESCRIPTION of Object

=head1 SYNOPSIS

use Bio::EnsEMBL::DBSQL::DBAdaptor;


my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor( -user   => 'root', 
                                             -db     => 'pog' , 
                                             -host   => 'caldy' , 
                                             -driver => 'mysql' );

my $protein_adaptor = $db->get_ProteinAdaptor();

my $protein = $protein_adaptor->fetch_by_transcript_id(1234);



=head1 DESCRIPTION

An adaptor which allows for the creation of Protein Objects from transcript
or translation ids. There is no Protein table in the EnsEMBL core database
and possibly this functionality should be moved into the translation object.

=head1 CONTACT

mongin@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal 
methods are usually preceded with a _

=cut


# Let the code begin...
package Bio::EnsEMBL::DBSQL::ProteinAdaptor;
use vars qw(@ISA);
use strict;


use Bio::EnsEMBL::DBSQL::BaseAdaptor;

use Bio::EnsEMBL::Protein;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);


=head2 fetch_by_transcript_id

  Arg [1]    : int $transid
               The unique internal identifier of this proteins transcript 
  Example    : $protein = $protein_adaptor->fetch_by_transcript_id(1234);
  Description: (formerly fetch_Protein_by_transcriptId) Retrieves a protein
               object via the internal database identifier of its transcript
  Returntype : Bio::EnsEMBL::Protein
  Exceptions : none
  Caller     : protview?

=cut

sub fetch_by_transcript_id{
   my ($self,$transid) = @_;
   my $query = "SELECT	t.translation_id 
		FROM    translation as t, 
			transcript_stable_id as s 
		WHERE	s.stable_id = '$transid' 
		AND	t.transcript_id = s.transcript_id";
   my $sth = $self->prepare($query);
   $sth->execute();
   my @row = $sth->fetchrow;
   
   return $self->fetch_by_translation_id($row[0]);
 }


=head2 fetch_by_translation_stable_id

  Arg [1]    : int $transid
               the stable identifier of the translation of the desired protein
  Example    : $prot = $pa->fetch_by_translation_stable_id('ENSP00000278194');
  Description: Creates a protein object using the translation table of the DB
  Returntype : Bio::EnsEMBL::Protein
  Exceptions : none
  Caller     : protview?

=cut

sub fetch_by_translation_stable_id {
   my ($self,$transid) = @_;
   my $query = "SELECT	translation_id 
		FROM	translation_stable_id
		WHERE	stable_id = '$transid'";
   my $sth = $self->prepare($query);
   $sth->execute();
   my @row = $sth->fetchrow;
   return $self->fetch_by_translation_id($row[0]);
 }



=head2 fetch_by_translation_id

  Arg [1]    : int $translation_id
               the unique DB identifier for the translation corresponding to 
               the desired protein
  Example    : $prot = $prot_adaptor->fetch_by_translation_id($id);
  Description: Retrieves a protein object from a database 
               (formerly fetch_Protein_by_dbid)
  Returntype : Bio::EnsEMBL::Protein 
  Exceptions : thrown if transcript or gene cannot be retrieved from database,
               or if the amino acid sequence cannot be obtained 
  Caller     : protview?

=cut

sub fetch_by_translation_id {
   my ($self, $translation_id) = @_;

   #Get the transcript id from the translation id 
   my $query = "SELECT	ts.transcript_id, 
			ts.gene_id 
		FROM	translation tr, transcript ts  
	        WHERE	translation_id = '$translation_id'
                AND     tr.transcript_id = ts.transcript_id";

   my $sth = $self->prepare($query);
   $sth ->execute();

   my ($transcript_id, $gene_id) = $sth->fetchrow();

   if (!defined $transcript_id) {
       $self->throw("$translation_id does not have a transcript_id");
   }

   if (!defined $gene_id) {
       $self->throw("$translation_id does not have a gene_id");
   }

   #Get Transcript object ( will allow us to get the aa sequence of the protein
   my $ta = $self->db()->get_TranscriptAdaptor();
   my $transcript = $ta->fetch_by_dbID($transcript_id);
 
   my $ga = $self->db()->get_GeneAdaptor();
   my $gene = $ga->fetch_by_dbID($gene_id);
   
   #Get all of the family (at the Transcript level), not implemented yet
   #my $family = $self->fetch_Family_by_dbid($id);

   #Get the polypeptide sequence out of the transcript object   
   my $sequence = $transcript->translate->seq;

   #Calculate the length of the Peptide   
   my $length = length($sequence);
   if ($length == 0) {
     $self->throw("Transcript " . $transcript->stable_id() .
		  " does not have an amino acid sequence"); 
   }
  
   #Define the moltype
   my $moltype = "protein";
   my $meta_obj = $self->db->get_MetaContainer();
   my $species = $meta_obj->get_Species();
   my $desc = "Protein predicted by Ensembl";

   #Create the Protein object
   my $protein = Bio::EnsEMBL::Protein->new ( -seq =>$sequence,
					  -accession_number => $translation_id,
					  -display_id => $translation_id,
					  -primary_id => $translation_id,
					  -id => $translation_id,
					  -desc => $desc,
					  -moltype => $moltype
					      );

   #Set up the adaptor handler for the protein object
   $protein->adaptor($self);

   #Add the species object to protein object
   $protein->species($species);

   $protein->transcript($transcript);
   $protein->gene($gene);

   return $protein;
}

1;


