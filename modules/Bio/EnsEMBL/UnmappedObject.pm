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

Bio:EnsEMBL::UnmappedObject - A object representing why a particular entity
was NOT mapped to the ensembl.

=head1 SYNOPSIS

  use Bio::EnsEMBL::UnmappedObject;

  my $uo = Bio::EnsEMBL::UnmappedObject->new(
    -type           => 'xref',
    -analysis       => $analysis,
    -external_db_id => 4100,
    -identifier     => "Q12345",
    -query_score    => 45.5,
    -target_score   => 29.2,
    -ensembl_id     => 122346,
    -ensembl_type   => "Translation",
    -summary        => "match failed for exonerate",
    -full_desc      => "match failed for the xref exonerate run "
      . "as match was below threshold of 90"
  );

=head1 DESCRIPTION

UnmappedObjects represent entities NOT mapped to ensembl. Therefore this
should help users to find out why certain accessions etc can not be
found.

=head1 METHODS

=cut



package Bio::EnsEMBL::UnmappedObject;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Storable;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Storable);


=head2 new

  Arg [TYPE]           : the type of mapping i.e. 'xref','cDNA' 
  Arg [ANALYSIS]       : Analysis object. 
  Arg [EXTERNAL_DB_ID] : id for the external db id this identifier is from
  Arg [IDENTIFIER]     : name of the identifier i.e. accession
  Arg [QUERY_SCORE]    : (optional) The query score
  Arg [TARGET_SCORE]   : (optional) The target score
  Arg [SUMMARY]        : The summary reason for not mapping.
  Arg [FULL_DESC]      : The Full description of why it did not map.
  Arg [ENSEMBL_ID]     : (optional) internal ensembl id for the best match
  Arg [ENSEMBL_OBJECT_TYPE] : (optional) the type of object for the best match
  Example              : see SYNOPSIS
  Returntype           : Bio::EnsEMBL::UnmappedObject
  Exceptions           : If any of the none optional args are missing
  Caller               : general
  Status               : At Risk

=cut

sub new {
  my $caller = shift;

  #allow constructor to be called as class or object method
  my $class = ref($caller) || $caller;



  my ($dbID, $unmapped_reason_id, $type, $analysis, $ex_db_id, $identifier, 
      $query_score, $target_score, $summary, $full_desc,
      $ensembl_id, $ensembl_object_type, $adaptor ) = 
    rearrange([qw(UNMAPPED_OBJECT_ID UNMAPPED_REASON_ID TYPE ANALYSIS 
		  EXTERNAL_DB_ID IDENTIFIER QUERY_SCORE TARGET_SCORE 
		  SUMMARY FULL_DESC ENSEMBL_ID ENSEMBL_OBJECT_TYPE ADAPTOR)], @_);

  my $self = $class->SUPER::new(@_);
  if(defined($analysis)) {
    if(!ref($analysis) || !$analysis->isa('Bio::EnsEMBL::Analysis')) {
      throw('-ANALYSIS argument must be a Bio::EnsEMBL::Analysis not '.
            $analysis);
    }
  }
  else{
    throw('-ANALYSIS argument must be given');
  }
  $self->{'analysis'} = $analysis;
  $self->{'dbID'} = $dbID if (defined($dbID));
  $self->{'description'} = $full_desc   || throw('FULL_DESC must be given');
  $self->{'summary'} = $summary         || throw('SUMMARY must be given');
  $self->{'type'} = $type               || throw('TYPE must be given');
  $self->{'external_db_id'} = $ex_db_id;

  if (lc($type) eq "xref") {
    throw('EXTERNAL_DB_ID must be given') if ! defined $ex_db_id;
  }

  $self->{'identifier'} = $identifier   || throw('IDENTIFIER must be given');
  $self->{'query_score'} = $query_score  if(defined($query_score));
  $self->{'target_score'} = $target_score  if(defined($target_score));
  $self->{'ensembl_id'} = $ensembl_id  if(defined($ensembl_id));
  $self->{'ensembl_object_type'} = $ensembl_object_type  
    if(defined($ensembl_object_type));
  $self->{'unmapped_reason_id'} = $unmapped_reason_id 
    if(defined($unmapped_reason_id));
  $self->adaptor($adaptor)  if(defined($adaptor));
  return $self;
}


=head2 description

  Arg [1]     : (optional) * to be set to
  Example     : print $unmappedObject->description."\n";
  Description : Basic getter/setter for description
  ReturnType  : String
  Exceptions  : none
  Caller      : general
  Status      : At Risk

=cut

sub description{
  my $self = shift;

  if(@_) {
    my $des = shift;
    $self->{'description'} = $des;
  }

  return $self->{'description'};
}

=head2 summary

  Arg [1]     : (optional) summary to be set to
  Example     : print $unmappedObject->summary."\n";
  Description : Basic getter/setter for summary
  ReturnType  : String
  Exceptions  : none
  Caller      : general
  Status      : At Risk

=cut

sub summary{
  my $self = shift;

  if(@_) {
    my $des = shift;
    $self->{'summary'} = $des;
  }

  return $self->{'summary'};
}

=head2 type

  Arg [1]     : (optional) type to be set to
  Example     : print $unmappedObject->type."\n";
  Description : Basic getter/setter for type
  ReturnType  : String
  Exceptions  : none
  Caller      : general
  Status      : At Risk

=cut

sub type{
  my $self = shift;

  if(@_) {
    my $arg = shift;
    $self->{'type'} = $arg;
  }

  return $self->{'type'};
}

=head2 ensembl_object_type

  Arg [1]     : (optional) ensembl object type to be set to
  Example     : print $unmappedObject->ensembl_object_type."\n";
  Description : Basic getter/setter for ensembl object type
  ReturnType  : String
  Exceptions  : none
  Caller      : general
  Status      : At Risk

=cut

sub ensembl_object_type{
  my $self = shift;

  if(@_) {
    my $arg = shift;
    $self->{'ensembl_object_type'} = $arg;
  }

  return $self->{'ensembl_object_type'};
}

=head2 ensembl_id

  Arg [1]     : (optional) ensembl id to be set to
  Example     : print $unmappedObject->ensembl_id."\n";
  Description : Basic getter/setter for ensembl id
  ReturnType  : String
  Exceptions  : none
  Caller      : general
  Status      : At Risk

=cut

sub ensembl_id{
  my $self = shift;

  if(@_) {
    my $arg = shift;
    $self->{'ensembl_id'} = $arg;
  }

  return $self->{'ensembl_id'};
}

=head2 external_db_id

  Arg [1]     : (optional) external_db_id to be set to
  Example     : print $unmappedObject->external_db_id."\n";
  Description : Basic getter/setter for external_db_id
  ReturnType  : int
  Exceptions  : none
  Caller      : general
  Status      : At Risk

=cut

sub external_db_id{
  my $self = shift;

  if(@_) {
    my $arg = shift;
    $self->{'external_db_id'} = $arg;
  }
  
  return $self->{'external_db_id'};
}

=head2 external_db_name

  Example     : print $unmappedObject->external_db_name()."\n";
  Description : Basic getter for external_db_name
  ReturnType  : int
  Exceptions  : none
  Caller      : general
  Status      : At Risk

=cut

sub external_db_name{
  my $self = shift;

  my $handle = $self->adaptor;
  if(defined($handle) and defined($self->{'external_db_id'})){
    my $sth = $handle->prepare("select db_name from external_db where external_db_id = ".$self->{'external_db_id'});
    $sth->execute();
    my $name;
    $sth->bind_columns(\$name);
    $sth->fetch();
    return $name;
  }
  return "";
}


sub stable_id{
  my ($self) = shift;
  
  my $handle = $self->adaptor;
  if(defined($handle)){
    my $sql = "select stable_id from ".lc($self->{'ensembl_object_type'})." where ".
      lc($self->{'ensembl_object_type'})."_id = ".
	$self->{'ensembl_id'};
    my $sth = $handle->prepare($sql);
    $sth->execute();
    my $name;
    $sth->bind_columns(\$name);
    $sth->fetch();
    return $name;
  }
  return "";
}
					     
 #  my $adaptor;
#  if($self->{'ensembl_object_type'} eq  "Transcript"){
#    $adaptor= $self->adaptor->db->get_TranscriptAdaptor();
#  }
#  elsif($self->{'ensembl_object_type'} eq  "Translation"){
#    $adaptor= $self->adaptor->db->get_TranslationAdaptor();
#  }
#  elsif($self->{'ensembl_object_type'} eq  "Gene"){
#    $adaptor= $self->adaptor->db->get_GeneAdaptor();
#  } 
#  else{
#    return undef;
#  }
#  my $object = $adaptor->fetch_by_dbID($self->{'ensembl_id'});
#  if(defined($object)){
#    return $object->stable_id;
#  }
#  else{
#    return undef;
#  }
#}

 
=head2 identifier

  Arg [1]     : (optional) identifier to be set to
  Example     : print $unmappedObject->identifier."\n";
  Description : Basic getter/setter for identifier
  ReturnType  : String
  Exceptions  : none
  Caller      : general
  Status      : At Risk

=cut

sub identifier{
  my $self = shift;

  if(@_) {
    my $arg = shift;
    $self->{'identifier'} = $arg;
  }

  return $self->{'identifier'};
}

=head2 query_score

  Arg [1]     : (optional) query_score to be set to
  Example     : print $unmappedObject->query_score."\n";
  Description : Basic getter/setter for query_score
  ReturnType  : float
  Exceptions  : none
  Caller      : general
  Status      : At Risk

=cut

sub query_score{
  my $self = shift;

  if(@_) {
    my $arg = shift;
    $self->{'query_score'} = $arg;
  }

  return $self->{'query_score'};
}

=head2 target_score

  Arg [1]     : (optional) target_score to be set to
  Example     : print $unmappedObject->target_score."\n";
  Description : Basic getter/setter for target_score
  ReturnType  : float
  Exceptions  : none
  Caller      : general
  Status      : At Risk

=cut

sub target_score{
  my $self = shift;

  if(@_) {
    my $arg = shift;
    $self->{'target_score'} = $arg;
  }

  return $self->{'target_score'};
}

=head2 unmapped_reason_id

  Arg [1]     : (optional) unmapped_reason_id to be set to
  Example     : print $unmappedObject->unmapped_reason_id."\n";
  Description : Basic getter/setter for unmapped_reason_id
  ReturnType  : int
  Exceptions  : none
  Caller      : general
  Status      : At Risk

=cut

sub unmapped_reason_id{
  my $self = shift;

  if(@_) {
    my $arg = shift;
    $self->{'unmapped_reason_id'} = $arg;
  }

  return $self->{'unmapped_reason_id'};
}

=head2 analysis

  Arg [1]     : (optional) analysis to be set to
  Example     : print $unmappedObject->analysis->logic_name."\n";
  Description : Basic getter/setter for analysis
  ReturnType  : Bio::EnsEMBL::Analysis
  Exceptions  : none
  Caller      : general
  Status      : At Risk

=cut

sub analysis {
  my $self = shift;

  if(@_) {
    my $an = shift;
    if(defined($an) && (!ref($an) || !$an->isa('Bio::EnsEMBL::Analysis'))) {
      throw('analysis argument must be a Bio::EnsEMBL::Analysis');
    }
    $self->{'analysis'} = $an;
  }

  return $self->{'analysis'};
}

1;
