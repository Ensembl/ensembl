#
# Ensembl module for Bio::EnsEMBL::DBSQL::AttributeAdaptor
#
# Copyright (c) 2003 EnsEMBL
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::AttributeAdaptor - Provides database interaction for
Bio::EnsEMBL::Attribute objects.


=head1 SYNOPSIS

  #$db is a Bio::EnsEMBL::DBSQL::DBAdaptor object:
  $attribute_adaptor = $db->get_AttributeAdaptor();

  $attributes = $attribute_adaptor->fetch_all_by_MiscFeature( $feature );

  $attributes = $attribute_adaptor->fetch_all_by_Slice( $slice );

  $attribute_adaptor->store_on_Slice( $slice, \@attributes );

  $attribute_adaptor->store_on_MiscFeature( $misc_feature, \@attributes )


=head1 DESCRIPTION


=head1 CONTACT

This modules is part of the Ensembl project http://www.ensembl.org

Questions can be posted to the ensembl-dev mailing list:
ensembl-dev@ebi.ac.uk

=cut

use strict;
use warnings;

package Bio::EnsEMBL::DBSQL::AttributeAdaptor;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Attribute;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);


=head2 new

  Arg [...]  : Superclass args.  See Bio::EnsEMBL::DBSQL::BaseAdaptor
  Description: Instantiates a Bio::EnsEMBL::DBSQL::AttributeAdaptor
  Returntype : Bio::EnsEMBL::AttributeAdaptor
  Exceptions : none
  Caller     : DBAdaptor
  Status     : Stable

=cut


sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);


  # cache creation could go here
  return $self;
}



=head2 fetch_all_by_MiscFeature

  Arg [1]    : Bio::EnsEMBL::MiscFeature $mf
  Example    : @attributes = @{$attrib_adaptor->fetch_all_by_MiscFeature($mf)};
  Description: Fetches all attributes for a given MiscFeature
  Returntype : Bio::EnsEMBL::Attribute
  Exceptions : throw if incorrect arguments
               throw if provided MiscFeature does not have a dbID
  Caller     : MiscFeature
  Status     : Stable

=cut

sub fetch_all_by_MiscFeature {
  my $self = shift;
  my $mf   = shift;

  if(!ref($mf) || !$mf->isa('Bio::EnsEMBL::MiscFeature')) {
    throw('MiscFeature argument is required.');
  }

  my $mfid = $mf->dbID();

  if(!defined($mfid)) {
    throw("MiscFeature must have dbID.");
  }

  my $sth = $self->prepare("SELECT at.code, at.name, at.description, " .
                           "       ma.value " .
                           "FROM misc_attrib ma, attrib_type at " .
                           "WHERE ma.misc_feature_id = ? " .
                           "AND   at.attrib_type_id = ma.attrib_type_id");

  $sth->bind_param(1,$mfid,SQL_INTEGER);
  $sth->execute();

  my $results = $self->_obj_from_sth($sth);

  $sth->finish();

  return $results;
}


=head2 fetch_all_by_Gene

  Arg [1]    : Bio::EnsEMBL::Gene $gene
  Example    : @attributes = @{ $attrib_adaptor->fetch_all_by_Gene($gene) };
  Description: Fetches all attributes for a given Gene
  Returntype : Bio::EnsEMBL::Attribute
  Exceptions : throw if incorrect arguments
               throw if provided MiscFeature does not have a dbID
  Caller     : Gene
  Status     : Stable

=cut

sub fetch_all_by_Gene {
  my $self = shift;
  my $gene = shift;

  if(!ref($gene) || !$gene->isa('Bio::EnsEMBL::Gene')) {
    throw('Gene argument is required.');
  }

  my $gid = $gene->dbID();

  if(!defined($gid)) {
    throw("Gene must have dbID.");
  }

  my $sth = $self->prepare("SELECT at.code, at.name, at.description, " .
                           "       ga.value " .
                           "FROM gene_attrib ga, attrib_type at " .
                           "WHERE ga.gene_id = ? " .
                           "AND   at.attrib_type_id = ga.attrib_type_id");

  $sth->execute($gid);

  my $results = $self->_obj_from_sth($sth);

  $sth->finish();

  return $results;
}


=head2 fetch_all_by_Transcript

  Arg [1]    : Bio::EnsEMBL::Transcript $tr
  Example    : @attributes = @{$attrib_adaptor->fetch_all_by_Transcript( $tr )};
  Description: Fetches all attributes for a given Transcript
  Returntype : Bio::EnsEMBL::Attribute
  Exceptions : throw if incorrect arguments
               throw if provided MiscFeature does not have a dbID
  Caller     : Transcript
  Status     : Stable

=cut

sub fetch_all_by_Transcript {
  my $self = shift;
  my $tr   = shift;

  if(!ref($tr) || !$tr->isa('Bio::EnsEMBL::Transcript')) {
    throw('Transcript argument is required.');
  }

  my $trid = $tr->dbID();

  if(!defined($trid)) {
    throw("Transcript must have dbID.");
  }

  my $sth = $self->prepare("SELECT at.code, at.name, at.description, " .
                           "       ta.value " .
                           "FROM transcript_attrib ta, attrib_type at " .
                           "WHERE ta.transcript_id = ? " .
                           "AND   at.attrib_type_id = ta.attrib_type_id");

  $sth->bind_param(1,$trid,SQL_INTEGER);
  $sth->execute();

  my $results = $self->_obj_from_sth($sth);

  $sth->finish();

  return $results;
}


=head2 fetch_all_by_Translation

  Arg [1]    : Bio::EnsEMBL::Translation $tl
  Example    : @attributes = @{$attrib_adaptor->fetch_all_by_Translation( $tl )};
  Description: Fetches all attributes for a given Translation
  Returntype : Bio::EnsEMBL::Attribute
  Exceptions : throw if incorrect arguments
               throw if provided Translation does not have a dbID
  Caller     : Transcript
  Status     : Stable

=cut

sub fetch_all_by_Translation {
  my $self = shift;
  my $tl   = shift;

  if(!ref($tl) || !$tl->isa('Bio::EnsEMBL::Translation')) {
    throw('Translation argument is required.');
  }

  my $tlid = $tl->dbID();

  if(!defined($tlid)) {
    throw("Translation must have dbID.");
  }

  my $sth = $self->prepare("SELECT at.code, at.name, at.description, " .
                           "       ta.value " .
                           "FROM translation_attrib ta, attrib_type at " .
                           "WHERE ta.translation_id = ? " .
                           "AND   at.attrib_type_id = ta.attrib_type_id");

  $sth->bind_param(1,$tlid,SQL_INTEGER);
  $sth->execute();

  my $results = $self->_obj_from_sth($sth);

  $sth->finish();

  return $results;
}




=head2 fetch_all_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
  Example    : @attributes = @{$attrib_adaptor->fetch_all_by_Slice($slice)};
  Description: Fetches all attributes for a given sequence region (which the
               passed in slice is on)
  Returntype : Bio::EnsEMBL::Attribute
  Exceptions : throw if incorrect arguments
               throw if cannot get seq_region_id from provided Slice
  Caller     : Slice
  Status     : Stable

=cut

sub fetch_all_by_Slice {
  my $self = shift;
  my $slice = shift;

  if(!ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
    throw('Slice argument is required.');
  }

  my $seq_region_id = $slice->get_seq_region_id();
  if(!defined($seq_region_id)) {
    throw("Could not get seq_region_id for provided slice: ".$slice->name());
  }

  my $sth = $self->prepare("SELECT at.code, at.name, at.description, " .
                           "       sra.value " .
                           "FROM seq_region_attrib sra, attrib_type at " .
                           "WHERE sra.seq_region_id = ? " .
                           "AND   at.attrib_type_id = sra.attrib_type_id");

  $sth->bind_param(1,$seq_region_id,SQL_INTEGER);
  $sth->execute();

  my $results = $self->_obj_from_sth($sth);
  $sth->finish();

  return $results;
}



=head2 store_on_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
  Arg [2]    : listref of Bio::EnsEMBL::Attribute objects $attribs
  Example    : $attribute_adaptor->store_on_Slice($slice, \@attribs);
  Description: Stores a set of attributes on a sequence region given a
               Slice object which is on the seq_region for which attributes are
               being stored.
  Returntype : none
  Exceptions : throw if $slice argument not provided
  Caller     : general
  Status     : Stable

=cut

sub store_on_Slice {
  my $self     = shift;
  my $slice    = shift;
  my $attribs = shift;

  if(!ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
    throw("Slice argument expected.");
  }

  if(ref($attribs) ne 'ARRAY') {
    throw("Reference to list of Bio::EnsEMBL::Attribute objects " .
          "argument expected.");
  }

  my $seq_region_id = $slice->get_seq_region_id();

  if(!$seq_region_id) {
    throw("Could not get seq_region_id for provided slice: ".$slice->name());
  }

  my $sth = $self->prepare( "INSERT into seq_region_attrib ".
                            "SET seq_region_id = ?, attrib_type_id = ?, ".
                            "value = ? " );

  foreach my $at ( @$attribs ) {
    if(!ref($at) && $at->isa('Bio::EnsEMBL::Attribute')) {
      throw("Reference to list of Bio::EnsEMBL::Attribute objects " .
            "argument expected.");
    }
    my $atid = $self->_store_type( $at );
    $sth->bind_param(1,$seq_region_id,SQL_INTEGER);
    $sth->bind_param(2,$atid,SQL_INTEGER);
    $sth->bind_param(3,$at->value,SQL_VARCHAR);
    $sth->execute();
  }

  return;
}




=head2 store_on_MiscFeature

  Arg [1]    : Bio::EnsEMBL::MiscFeature $feature
  Example    : $attribute_adaptor->store_on_MiscFeature($my_misc_feature,
                                                        $attributes)
  Description: Stores all the attributes on the misc feature. 
               Will duplicate things if called twice.
  Returntype : none
  Exceptions : throw on incorrect arguments
               throw if provided feature is not stored in this database
  Caller     : general, MiscFeatureAdaptor
  Status     : Stable

=cut

sub store_on_MiscFeature {
  my $self       = shift;
  my $feature    = shift;
  my $attributes = shift;

  if(!ref($feature) || !$feature->isa('Bio::EnsEMBL::MiscFeature')) {
    throw("MiscFeature argument expected");
  }

  if(ref($attributes) ne 'ARRAY') {
    throw("Reference to list of Bio::EnsEMBL::Attribute objects argument " .
          "expected");
  }

  my $db = $self->db();
  if(!$feature->is_stored($db)) {
    throw("MiscFeature is not stored in this DB - cannot store attributes.");
  }

  my $feature_id = $feature->dbID();

  my $sth = $self->prepare( "INSERT into misc_attrib ".
			    "SET misc_feature_id = ?, attrib_type_id = ?, ".
			    "value = ? " );

  for my $attrib ( @$attributes ) {
    if(!ref($attrib) && $attrib->isa('Bio::EnsEMBL::Attribute')) {
      throw("Reference to list of Bio::EnsEMBL::Attribute objects " .
            "argument expected.");
    }
    my $atid = $self->_store_type( $attrib );
    $sth->bind_param(1,$feature_id,SQL_INTEGER);
    $sth->bind_param(2,$atid,SQL_INTEGER);
    $sth->bind_param(3,$attrib->value,SQL_VARCHAR);
    $sth->execute();
  }

  return;
}


=head2 store_on_Gene

  Arg [1]    : Bio::EnsEMBL::Gene $gene
  Arg [2]    : listref Bio::EnsEMBL::Attribute $attribs
  Example    : $attribute_adaptor->store_on_Gene($my_gene, $attributes)
  Description: Stores all the attributes on the Gene. 
               Will duplicate things if called twice.
  Returntype : none
  Exceptions : throw on incorrect arguments
               throw if provided Gene is not stored in this database
  Caller     : general, GeneAdaptor
  Status     : Stable

=cut

sub store_on_Gene {
  my $self = shift;
  my $gene = shift;
  my $attributes = shift;

  if(!ref($gene) || !$gene->isa('Bio::EnsEMBL::Gene')) {
    throw("Gene argument expected");
  }

  if(ref($attributes) ne 'ARRAY') {
    throw("Reference to list of Bio::EnsEMBL::Attribute objects argument " .
          "expected");
  }

  my $db = $self->db();
  if(!$gene->is_stored($db)) {
    throw("Gene is not stored in this DB - cannot store attributes.");
  }

  my $gene_id = $gene->dbID();

  my $sth = $self->prepare( "INSERT into gene_attrib ".
			    "SET gene_id = ?, attrib_type_id = ?, ".
			    "value = ? " );

  for my $attrib ( @$attributes ) {
    if(!ref($attrib) && $attrib->isa('Bio::EnsEMBL::Attribute')) {
      throw("Reference to list of Bio::EnsEMBL::Attribute objects " .
            "argument expected.");
    }
    my $atid = $self->_store_type( $attrib );
    $sth->execute( $gene_id, $atid, $attrib->value() );
  }

  return;
}


=head2 store_on_Transcript

  Arg [1]    : Bio::EnsEMBL::Transcript $transcript
  Arg [2]    : listref Bio::EnsEMBL::Attribute $attribs
  Example    : $attribute_adaptor->store_on_Transcript($my_transcript,
                                                        $attributes)
  Description: Stores all the attributes on the Transcript. 
               Will duplicate things if called twice.
  Returntype : none
  Exceptions : throw on incorrect arguments
               throw if provided Transcript is not stored in this database
  Caller     : general, TranscriptAdaptor
  Status     : Stable

=cut

sub store_on_Transcript {
  my $self       = shift;
  my $transcript    = shift;
  my $attributes = shift;

  if(!ref($transcript) || !$transcript->isa('Bio::EnsEMBL::Transcript')) {
    throw("Transcript argument expected");
  }

  if(ref($attributes) ne 'ARRAY') {
    throw("Reference to list of Bio::EnsEMBL::Attribute objects argument " .
          "expected");
  }

  my $db = $self->db();
  if(!$transcript->is_stored($db)) {
    throw("Transcript is not stored in this DB - cannot store attributes.");
  }

  my $transcript_id = $transcript->dbID();

  my $sth = $self->prepare( "INSERT into transcript_attrib ".
			    "SET transcript_id = ?, attrib_type_id = ?, ".
			    "value = ? " );

  for my $attrib ( @$attributes ) {
    if(!ref($attrib) && $attrib->isa('Bio::EnsEMBL::Attribute')) {
      throw("Reference to list of Bio::EnsEMBL::Attribute objects " .
            "argument expected.");
    }
    my $atid = $self->_store_type( $attrib );
    $sth->bind_param(1,$transcript_id,SQL_INTEGER);
    $sth->bind_param(2,$atid,SQL_INTEGER);
    $sth->bind_param(3,$attrib->value,SQL_VARCHAR);
    $sth->execute();
  }

  return;
}



=head2 store_on_Translation

  Arg [1]    : Bio::EnsEMBL::Translation $translation
  Arg [2]    : listref Bio::EnsEMBL::Attribute $attribs
  Example    : $attribute_adaptor->store_on_Translation($my_translation,
                                                        $attributes)
  Description: Stores all the attributes on the Translation. 
               Will duplicate things if called twice.
  Returntype : none
  Exceptions : throw on incorrect arguments
               throw if provided Translation is not stored in this database
  Caller     : general, TranslationAdaptor
  Status     : Stable

=cut

sub store_on_Translation {
  my $self        = shift;
  my $translation = shift;
  my $attributes  = shift;

  if(!ref($translation) || !$translation->isa('Bio::EnsEMBL::Translation')) {
    throw("Translation argument expected");
  }

  if(ref($attributes) ne 'ARRAY') {
    throw("Reference to list of Bio::EnsEMBL::Attribute objects argument " .
          "expected");
  }

  my $db = $self->db();
  if(!$translation->is_stored($db)) {
    throw("Translation is not stored in this DB - cannot store attributes.");
  }

  my $translation_id = $translation->dbID();

  my $sth = $self->prepare( "INSERT into translation_attrib ".
			    "SET translation_id = ?, attrib_type_id = ?, ".
			    "value = ? " );

  for my $attrib ( @$attributes ) {
    if(!ref($attrib) && $attrib->isa('Bio::EnsEMBL::Attribute')) {
      throw("Reference to list of Bio::EnsEMBL::Attribute objects " .
            "argument expected.");
    }
    my $atid = $self->_store_type( $attrib );
    $sth->bind_param(1,$translation_id,SQL_INTEGER);
    $sth->bind_param(2,$atid,SQL_INTEGER);
    $sth->bind_param(3,$attrib->value,SQL_VARCHAR);
    $sth->execute();
  }

  return;
}


=head2 remove_from_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               The Slice to remove attributes from
  Arg [2]    : (optional) dbID $attrib_type_id
               Database id of attributes to remove
  Example    : $attribute_adaptor->remove_from_Slice($slice);
  Description: Removes all attributes which are stored in the database on
               the provided Slice.
  Returntype : none
  Exceptions : throw on incorrect arguments
               throw if cannot obtain seq_region_id from provided Slice
  Caller     : general
  Status     : Stable

=cut

sub remove_from_Slice {
  my $self = shift;
  my $slice = shift;
  my $attrib_type_id = shift;

  if(!ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
    throw("Bio::EnsEMBL::Slice argument expected.");
  }

  my $srid = $slice->get_seq_region_id();

  if(!$srid) {
    throw("Could not get seq_region_id for provided slice: ".$slice->name());
  }
  
  if (defined $attrib_type_id) {
    my $sth = $self->prepare("DELETE FROM seq_region_attrib " .
        "WHERE seq_region_id = ? AND attrib_type_id = ?");
    $sth->bind_param(1,$srid,SQL_INTEGER);
    $sth->bind_param(2,$attrib_type_id,SQL_INTEGER);
    $sth->execute();
    $sth->finish;
  } else {
    my $sth = $self->prepare("DELETE FROM seq_region_attrib " .
        "WHERE seq_region_id = ?");
    $sth->bind_param(1,$srid,SQL_INTEGER);
    $sth->execute();
    $sth->finish();
  }

  return;
}


=head2 remove_from_MiscFeature

  Arg [1]    : Bio::EnsEMBL::MiscFeature $mf
               The MiscFeature to remove attributes from
  Example    : $attribute_adaptor->remove_from_MiscFeature($mf);
  Description: Removes all attributes which are stored in the database on
               the provided MiscFeature.
  Returntype : none
  Exceptions : throw on incorrect arguments
               throw if MiscFeature is not stored in this database
  Caller     : general
  Status     : Stable

=cut

sub remove_from_MiscFeature {
  my $self = shift;
  my $mf   = shift;

  if(!ref($mf) || !$mf->isa('Bio::EnsEMBL::MiscFeature')) {
    throw("Bio::EnsEMBL::MiscFeature argument is required.");
  }

  my $db = $self->db();

  if(!$mf->is_stored($db)) {
    throw("MiscFeature is not stored in this database.");
  }

  my $sth = $db->dbc->prepare("DELETE FROM misc_attrib " .
                         "WHERE misc_feature_id = ?");

  $sth->bind_param(1,$mf->dbID,SQL_INTEGER);
  $sth->execute();

  $sth->finish();

  return;
}


=head2 remove_from_Gene

  Arg [1]    : Bio::EnsEMBL::Gene $gene
               The Gene to remove attributes from
  Arg [2]    : (optional) dbID $attrib_type_id
               Database id of attributes to remove
  Example    : $attribute_adaptor->remove_from_Gene($gene);
  Description: Removes all attributes which are stored in the database on
               the provided Gene.
  Returntype : none
  Exceptions : throw on incorrect arguments
               throw if gene not stored
  Caller     : general
  Status     : Stable

=cut

sub remove_from_Gene {
  my $self = shift;
  my $gene = shift;
  my $attrib_type_id = shift;

  if(!ref($gene) || !$gene->isa('Bio::EnsEMBL::Gene')) {
    throw("Bio::EnsEMBL::Gene argument expected.");
  }

  my $db = $self->db();

  if (! $gene->is_stored($db)) {
    throw("Gene is not stored in the database.");
  }
  
  if (defined $attrib_type_id) {
    my $sth = $self->prepare("DELETE FROM gene_attrib " .
        "WHERE gene_id = ? AND attrib_type_id = ?");
    $sth->execute($gene->dbID, $attrib_type_id);
    $sth->finish;
  } else {
    my $sth = $self->prepare("DELETE FROM gene_attrib " .
        "WHERE gene_id = ?");
    $sth->execute($gene->dbID);
    $sth->finish();
  }

  return;
}


=head2 remove_from_Transcript

  Arg [1]    : Bio::EnsEMBL::Transcript $transcript
               The Transcript to remove attributes from
  Arg [2]    : (optional) dbID $attrib_type_id
               Database id of attributes to remove
  Example    : $attribute_adaptor->remove_from_Transcript($transcript);
  Description: Removes all attributes which are stored in the database on
               the provided Transcript.
  Returntype : none
  Exceptions : throw on incorrect arguments
               throw if transcript not stored
  Caller     : general
  Status     : Stable

=cut

sub remove_from_Transcript {
  my $self = shift;
  my $transcript = shift;
  my $attrib_type_id = shift;

  if(!ref($transcript) || !$transcript->isa('Bio::EnsEMBL::Transcript')) {
    throw("Bio::EnsEMBL::Transcript argument expected.");
  }

  my $db = $self->db();

  if (! $transcript->is_stored($db)) {
    throw("Transcript is not stored in the database.");
  }
  
  if (defined $attrib_type_id) {
    my $sth = $self->prepare("DELETE FROM transcript_attrib " .
        "WHERE transcript_id = ? AND attrib_type_id = ?");
    $sth->bind_param(1,$transcript->dbID,SQL_INTEGER);
    $sth->bind_param(2,$attrib_type_id,SQL_INTEGER);
    $sth->execute();
    $sth->finish;
  } else {
    my $sth = $self->prepare("DELETE FROM transcript_attrib " .
        "WHERE transcript_id = ?");
    $sth->bind_param(1,$transcript->dbID,SQL_INTEGER);
    $sth->execute();
    $sth->finish();
  }

  return;
}



=head2 remove_from_Translation

  Arg [1]    : Bio::EnsEMBL::Translation $translation
               The Translation to remove attributes from
  Arg [2]    : (optional) dbID $attrib_type_id
               Database id of attributes to remove
  Example    : $attribute_adaptor->remove_from_Transcript($transcript);
  Description: Removes all attributes which are stored in the database on
               the provided Transcript.
  Returntype : none
  Exceptions : throw on incorrect arguments
               throw if transcript not stored
  Caller     : general
  Status     : Stable

=cut

sub remove_from_Translation {
  my $self = shift;
  my $translation = shift;
  my $attrib_type_id = shift;

  if(!ref($translation) || !$translation->isa('Bio::EnsEMBL::Translation')) {
    throw("Bio::EnsEMBL::Translation argument expected.");
  }

  my $db = $self->db();

  if (! $translation->is_stored($db)) {
    throw("Transcript is not stored in the database.");
  }
  
  if (defined $attrib_type_id) {
    my $sth = $self->prepare("DELETE FROM translation_attrib " .
        "WHERE translation_id = ? AND attrib_type_id = ?");
    $sth->bind_param(1,$translation->dbID,SQL_INTEGER);
    $sth->bind_param(2,$attrib_type_id,SQL_INTEGER);
    $sth->execute();
    $sth->finish;
  } else {
    my $sth = $self->prepare("DELETE FROM translation_attrib " .
        "WHERE translation_id = ?");
    $sth->bind_param(1,$translation->dbID,SQL_INTEGER);
    $sth->execute();
    $sth->finish();
  }

  return;
}






# _store_type

sub _store_type {
  my $self = shift;
  my $attrib = shift;

  my $sth1 = $self->prepare
    ("INSERT IGNORE INTO attrib_type set code = ?, name = ?, ".
     "description = ?" );


  $sth1->bind_param(1,$attrib->code,SQL_VARCHAR);
  $sth1->bind_param(2,$attrib->name,SQL_VARCHAR);
  $sth1->bind_param(3,$attrib->description,SQL_LONGVARCHAR);

  my $rows_inserted =  $sth1->execute();

  my $atid = $sth1->{'mysql_insertid'};

  if($rows_inserted == 0) {
    # the insert failed because the code is already stored
    my $sth2 = $self->prepare
      ("SELECT attrib_type_id FROM attrib_type " .
       "WHERE code = ?");
    $sth2->bind_param(1,$attrib->code,SQL_VARCHAR);
    $sth2->execute();
    ($atid) = $sth2->fetchrow_array();

    $sth2->finish();

    if(!$atid) {
      throw("Could not store or fetch attrib_type code [".$attrib->code."]\n" .
	    "Wrong database user/permissions?");
    }
  }

  $sth1->finish();


  return $atid;
}




sub _obj_from_sth {
  my $self = shift;
  my $sth = shift;

  my ($code, $name, $desc, $value);
  $sth->bind_columns(\$code, \$name, \$desc, \$value);

  my @results;
  while($sth->fetch()) {
    push @results, Bio::EnsEMBL::Attribute->new(-CODE => $code,
                                                -NAME => $name,
                                                -DESCRIPTION => $desc,
                                                -VALUE => $value);
  }

  return \@results;
}


1;
