#
# Ensembl module for Bio::EnsEMBL::UnmappedObject
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio:EnsEMBL::UnmappedObject - A object representing why a particular entity
was NOT mapped to the ensembl.

=head1 SYNOPSIS

use Bio::EnsEMBL::UnmappedObject;

my $uo = Bio::EnsEMBL::UnmappedObject->new (
        -type           => 'xref',
        -analysis       => $analysis,
        -external_db_id => 4100,
        -identifier     => "Q12345",
        -query_score    => 45.5,
        -target_score   => 29.2,
        -ensembl_id     => 122346,
        -ensembl_type   => 'Translation",
        -summary        => "match failed for exonerate",
        -full_desc      => "match failed for the xref exonerate run as mtch was below threshold of 90");


=head1 DESCRIPTION

UnmappedObjects represent entities NOT mapped to ensembl. Therefore this should help
users to find out why certain accessions etc can not be found.

This module is part of the Ensembl project http://www.ensembl.org

=head1 CONTACT

Post comments/questions to the ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut



use strict;
use warnings;

package Bio::EnsEMBL::UnmappedObject;

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
      $ensembl_id, $ensembl_object_type) = 
    rearrange([qw(UNMAPPED_OBJECT_ID UNMAPPED_REASON_ID TYPE ANALYSIS 
		  EXTERNAL_DB_ID IDENTIFIER QUERY_SCORE TARGET_SCORE 
		  SUMMARY FULL_DESC ENSEMBL_ID ENSEMBL_OBJECT_TYPE)], @_);

  my $self = $class->SUPER::new(@_);
  if($analysis) {
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
  $self->{'external_db_id'} = $ex_db_id || throw('EXTERNAL_DB_ID must be given');
  $self->{'identifier'} = $identifier   || throw('IDENTIFIER must be given');
  $self->{'query_score'} = $query_score  if(defined($query_score));
  $self->{'target_score'} = $target_score  if(defined($target_score));
  $self->{'ensembl_id'} = $ensembl_id  if(defined($ensembl_id));
  $self->{'ensembl_object_type'} = $ensembl_object_type  
    if(defined($ensembl_object_type));
  $self->{'unmapped_reason_id'} = $unmapped_reason_id 
    if(defined($unmapped_reason_id));
  return $self;
}

=head2 new_fast

  Arg [...]  : none
  Example    : $feature = Bio::EnsEMBL::UnmappedObject->new_fast();
  Description: Creates a new Unmapped Object.
  Returntype : Bio::EnsEMBL::UnmappedObject
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub new_fast{
  my $caller = shift;

  #allow constructor to be called as class or object method
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

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
