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

Bio::EnsEMBL::DBSQL::UnmappedObjectAdaptor

=head1 SYNOPSIS

  my $uoa = $database_adaptor->get_UnmappedObjectAdaptor();

  my $missed = @{ $uoa->fetch_all_by_type('xref') };

=head1 DESCRIPTION

Unmapped ObjectAdaptor - An adaptor responsible for the creation,
editing, retrieval of Unmapped Objects. These being the Objects that
where not mapped in a specific process i.e. xref, cDNA, Markers.

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::UnmappedObjectAdaptor;
use vars qw(@ISA);
use strict;


use POSIX;
use Bio::EnsEMBL::Utils::Cache;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::UnmappedObject;
use Bio::EnsEMBL::Analysis;
@ISA = qw(Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor);


=head2 new

  Arg [1]    : list of args @args
               Superclass constructor arguments
  Example    : none
  Description: Constructor which just initializes internal cache structures
  Returntype : Bio::EnsEMBL::DBSQL::UnmappedObjectAdaptor
  Exceptions : none
  Caller     : implementing subclass constructors
  Status     : At Risk

=cut

our %desc_to_id;
sub new {
  my $proto = shift;

  my $class = ref($proto) || $proto;
  my $self = $class->SUPER::new(@_);

  my $sth =
      $self->prepare("
        SELECT 
          unmapped_reason_id, 
          full_description 
        FROM 
          unmapped_reason
        ");

    $sth->execute();

    my ( $id, $desc );
    $sth->bind_columns( \( $id, $desc ) );

    while ( $sth->fetch() ) {
      $desc_to_id{$desc} = $id;
    }

    $sth->finish();


  return $self;
}


# _tables
#  Arg [1]    : none
#  Description: PROTECTED implementation of superclass abstract method
#               returns the names, aliases of the tables to use for queries
#  Returntype : list of listrefs of strings
#  Exceptions : none
#  Caller     : internal
#  Status     : At Risk
sub _tables {
  my $self = shift;

  return (['unmapped_object', 'uo'],
	  ['unmapped_reason', 'ur']);
}


# _columns
#  Arg [1]    : none
#  Example    : none
#  Description: PROTECTED implementation of superclass abstract method
#               returns a list of columns to use for queries
#  Returntype : list of strings
#  Exceptions : none
#  Caller     : internal
#  Status     : At Risk

sub _columns {
  my $self = shift;

  return qw(uo.unmapped_object_id uo.type uo.analysis_id uo.external_db_id
            uo.identifier uo.unmapped_reason_id uo.query_score uo.target_score 
	    uo.ensembl_id uo.ensembl_object_type 
	    ur.summary_description ur.full_description);
}

sub _left_join {
  return ( [
      'unmapped_object', "uo.unmapped_reason_id = ur.unmapped_reason_id"
    ] );
}

=head2 list_dbIDs

  Arg [1]    : none
  Example    : @unmapped_object_ids = @{$unmapped_object_adaptor->list_dbIDs()};
  Description: Gets an array of internal ids for all unmapped_objects in the current db
  Returntype : list of ints
  Exceptions : none
  Caller     : ?
  Status     : Stable

=cut

sub list_dbIDs {
   my ($self) = @_;

   return $self->_list_dbIDs("unmapped_object");
}

=head2 list_unmapped_reasons

  Arg [1]    : none
  Example    : @unmapped_object_reason+ids = 
                   @{$unmapped_object_adaptor->list_unmapped_reasons()};
  Description: Gets an array of internal ids for all unmapped_objects in the current db
  Returntype : list of ints
  Exceptions : none
  Caller     : ?
  Status     : Stable

=cut

sub list_unmapped_reasons {
   my ($self) = @_;

   return $self->_list_dbIDs("unmapped_reason");
}


# _objs_from_sth

#  Arg [1]    : StatementHandle $sth
#  Example    : none
#  Description: PROTECTED implementation of abstract superclass method.
#               responsible for the creation of UnmappedObjects
#  Returntype : listref of Bio::EnsEMBL::UnmappedObjects
#  Exceptions : none
#  Caller     : internal
#  Status     : At Risk

sub _objs_from_sth {
  my ($self, $sth) = @_;

  my($unmapped_object_id, $type, $analysis_id, $external_db_id, $identifier,
     $unmapped_reason_id, $query_score, $target_score, $ensembl_id, 
     $ensembl_object_type, $summary, $full_desc);

  $sth->bind_columns(\$unmapped_object_id,\$type, \$analysis_id, 
		     \$external_db_id, \$identifier, \$unmapped_reason_id, 
		     \$query_score, \$target_score, \$ensembl_id, 
		     \$ensembl_object_type, \$summary, \$full_desc);

  my $analysis_adaptor = $self->db->get_AnalysisAdaptor();

  my @features;
  while($sth->fetch()) {
    my $analysis = $analysis_adaptor->fetch_by_dbID($analysis_id);
    
    #print "$identifier\n";

    push( @features,
          $self->_create_feature(
                         'Bio::EnsEMBL::UnmappedObject', {
                           -unmapped_object_id  => $unmapped_object_id,
                           -unmapped_reason_id  => $unmapped_reason_id,
                           -type                => $type,
                           -analysis            => $analysis,
                           -external_db_id      => $external_db_id,
                           -identifier          => $identifier,
                           -query_score         => $query_score,
                           -target_score        => $target_score,
                           -ensembl_id          => $ensembl_id,
                           -ensembl_object_type => $ensembl_object_type,
                           -summary             => $summary,
                           -full_desc           => $full_desc,
                           -adaptor             => $self
                         } ) );

  }
  return \@features;
}



=head2 store

  Arg [1]    : list of Bio::EnsEMBL::UnmappedObjects @uo
               the unmapped objects to store in the database
  Example    : $ou_adaptor->store(@uo);
  Description: Stores a list of unmapped objects in the database
  Returntype : none
  Exceptions : thrown if no Analysis, or no list of objects to store. 
  Caller     : general
  Status     : Stable

=cut

sub store{
  my ($self,@uos) = @_;

  if( scalar(@uos) == 0 ) {
    throw("Must call store with list of UnmappedObjects");
  }


  my $db = $self->db();
  my $analysis_adaptor = $db->get_AnalysisAdaptor();

  my $sth_reason = $self->prepare
    ("INSERT INTO unmapped_reason (summary_description, full_description)".
     " VALUES (?,?)");

  my $sth_unmapped_object = $self->prepare
	("INSERT INTO unmapped_object (type, analysis_id, external_db_id,
              identifier, unmapped_reason_id, query_score, target_score,
              ensembl_id, ensembl_object_type)".
	    " VALUES (?,?,?,?,?,?,?,?,?)");

  my $sth_fetch_reason =
    $self->prepare(
      "SELECT 
        unmapped_reason_id
      FROM 
        unmapped_reason
      WHERE
        full_description = ? 
        " );

  FEATURE: foreach my $uo ( @uos ) {

    if( !ref $uo || !$uo->isa("Bio::EnsEMBL::UnmappedObject") ) {
      throw("UnmappedObject must be an Ensembl UnmappedObject, " .
            "not a [".ref($uo)."]");
    }
    if($uo->is_stored($db)){
      next;
    } 
    
    my $analysis = $uo->analysis();
    throw("UnmappedObject must have an analysis object.".$uo->analysis."\n") if(!defined($analysis));

    my $analysis_id;
    if($analysis->is_stored($db)) {
      $analysis_id = $analysis->dbID();
    } 
    else {
      $analysis_id = $db->get_AnalysisAdaptor->store($analysis);
    }
    


    # Check if unmapped reason is cached
    # Check if it has been added since the cache was created
    # Try to store it

    if(!defined($desc_to_id{$uo->{'description'}})){
      $sth_reason->bind_param(1,$uo->{'summary'},SQL_VARCHAR);
      $sth_reason->bind_param(2,$uo->{'description'},SQL_VARCHAR);
      
      if(! eval{ $sth_reason->execute(); 1 }){    
        # DBI Trace possible here?
        warning($@); #
        my $msg;
        $msg .= "INSERT INTO unmapped_reason (summary_description, full_description) VALUES (";
        $msg .=  $uo->{'summary'} .','. $uo->{'description'}. ')';
        # Temporary fix for naughty cross-dependency regulation code.
        use Data::Dumper;
        warning("Query: \n$msg");
        print STDERR "UnmappedObject: \n";
        print STDERR Dumper $uo;

        $sth_fetch_reason->execute($uo->{'description'});
        
        my $unmapped_reasons = $sth_fetch_reason->fetchrow_arrayref();
        if(! defined($unmapped_reasons)){
          my $msg = $uo->{'description'}. " unable to store. Check MySQL schema, maybe PK not big enough?";
          throw($msg);
        }
        if(scalar @$unmapped_reasons != 1){
          throw("Multiple results for this description");
        }
        $uo->{'unmapped_reason_id'} = $unmapped_reasons->[0];

      }    
      else{
        # print 'Last mohican: '.$self->last_insert_id ."\n";
        $uo->{'unmapped_reason_id'} = $self->last_insert_id;
      }
    }
    else{
      $uo->{'unmapped_reason_id'} = $desc_to_id{$uo->{'description'}};
    }
    $sth_unmapped_object->bind_param(1,$uo->{'type'},SQL_VARCHAR);
    $sth_unmapped_object->bind_param(2,$uo->analysis->dbID,SQL_INTEGER);
    $sth_unmapped_object->bind_param(3,$uo->{'external_db_id'},SQL_INTEGER);
    $sth_unmapped_object->bind_param(4,$uo->{'identifier'},SQL_VARCHAR);
    $sth_unmapped_object->bind_param(5,$uo->{'unmapped_reason_id'},SQL_VARCHAR);
    $sth_unmapped_object->bind_param(6,$uo->{'query_score'},SQL_DOUBLE);
    $sth_unmapped_object->bind_param(7,$uo->{'target_score'},SQL_DOUBLE);
    $sth_unmapped_object->bind_param(8,$uo->{'ensembl_id'},SQL_INTEGER);
    $sth_unmapped_object->bind_param(9,$uo->{'ensembl_object_type'},SQL_VARCHAR);
    $sth_unmapped_object->execute();
    $uo->dbID($self->last_insert_id('unmapped_object_id', undef, 'unmapped_object'));
  }
  $sth_reason->finish();      
  return;
}


=head2 fetch_all_by_type

  Arg [1]  : string type. The type of unmapped objects
  Example  : @unmapped_object = @{$uoa->fetch_all_by_type('xref')};
  Description : Retrieves all the unmapped object for a particular
                type. e.g. 'xref','cDNA', 'marker'
  Returntype  : Array ref of Bio::EnsEMBL::UnmappedObject
  Exceptions  : none
  Caller      : general
  Status      : At Risk

=cut

sub fetch_all_by_type {
  my ($self, $type) = @_;
  
  unless($type) {
    throw("type argument is required");
  }
  $self->bind_param_generic_fetch($type,SQL_VARCHAR);
  $self->generic_fetch("uo.type = ?");
  
}

=head2 fetch_all_by_analysis

  Arg [1]  : Bio:EnsEMBL::Analysis object
  Arg [2]  : (optional) string database name
  Example  : @unmapped_object = @{$uoa->fetch_all_by_analysis($analysis)};
  Description : Retrieves all the unmapped object for a particular
                analysis type with the the option of a particular
                database type.
  Returntype  : array ref of Bio::EnsEMBL::UnmappedObject
  Exceptions  : thorws if first argument is not an anaylisi object
  Caller      : general
  Status      : At Risk

=cut

sub fetch_all_by_analysis {
  my ($self, $analysis,$dbname) = @_;
  
  unless($analysis) {
    throw("analysis argument is required");
  }
  $self->bind_param_generic_fetch($analysis->dbID,SQL_INTEGER);
  my $constraint = "uo.analysis_id = ?";
  if(defined($dbname)){
    my $db_id =0;
    my $sth = $self->prepare('select external_db_id from external_db where db_name like "'.
			      $dbname.'"');
    $sth->execute;
    $sth->bind_columns(\$db_id);
    $sth->fetch();
    if(!defined($db_id) or $db_id == 0){
      throw("$dbname could not be found in the external database table\n");
    }
    $self->bind_param_generic_fetch($db_id,SQL_INTEGER);
    $constraint .= " AND uo.external_db_id = ?";
  }
  $self->generic_fetch($constraint);
  
}

=head2 fetch_by_identifier

  Arg [1]  : string type. The type of unmapped objects
  Arg [2]  : (optional) string database name
  Example  : @unmapped_object = @{$uoa->fetch_by_identifier('Q123345')};
  Description : Retrieves the unmapped object for a particular
                identifier/accession
  Returntype  : array ref of Bio::EnsEMBL::UnmappedObject
  Exceptions  : none
  Caller      : general
  Status      : At Risk

=cut

sub fetch_by_identifier {
  my ($self, $identifier, $dbname) = @_;
  
  unless($identifier) {
    throw("identifier argument is required");
  }
  $self->bind_param_generic_fetch($identifier,SQL_VARCHAR);
  my $constraint = 'uo.identifier like ?';

  if(defined($dbname)){
    my $db_id =0;
    my $sth = $self->prepare('select external_db_id from external_db where db_name like "'.
			      $dbname.'"');
    $sth->execute;
    $sth->bind_columns(\$db_id);
    $sth->fetch();
    if(!defined($db_id) or $db_id == 0){
      throw("$dbname could not be found in the external database table\n");
    }
    $self->bind_param_generic_fetch($db_id,SQL_INTEGER);
    $constraint .= " AND uo.external_db_id = ?";
  }
  return $self->generic_fetch($constraint);
}

=head2 fetch_all_by_object_type_id

  Arg [1]  : string - The object type of the ensembl object e.g. Gene
  Arg [2]  : int    - The internal dbID of the ensembl object
  Example  : my @unmapped_objects = @{$uoa->fetch_all_by_object_type_id('Gene', 12341)};
  Description : Retrieves the unmapped objects for a particular ensembl object
                This is a base method which should be called by wrapper methods
                defining the correct object type e.g. $uoa->fetch_all_by_Gene($gene)
  Returntype  : array ref of Bio::EnsEMBL::UnmappedObject objects
  Exceptions  : Throws if arguments are not defined
  Caller      : general
  Status      : At Risk

=cut

sub fetch_all_by_object_type_id {
  my ($self, $object_type, $dbid) = @_;
  
  if(! ($object_type && $dbid)){
    throw("object_type and dbid arguments required");
  }
  
  $self->bind_param_generic_fetch($object_type, SQL_VARCHAR);
  $self->bind_param_generic_fetch($dbid,        SQL_INTEGER);

  my $constraint = 'uo.ensembl_object_type=? and uo.ensembl_id=?';

  return $self->generic_fetch($constraint);
}



1;
