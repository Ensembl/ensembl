# Perl module for Bio::EnsEMBL::DBSQL::AnalysisAdaptor
#
# Creator: Arne Stabenau <stabenau@ebi.ac.uk>
# Date of creation: 25.01.2001
# Last modified : 25.01.2001 by Arne Stabenau
#
# Copyright EMBL-EBI 2000
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::AnalysisAdaptor 

=head1 SYNOPSIS

  $analysisAdaptor = $db_adaptor->getAnalysisAdaptor;
  $analysisAdaptor = $analysisobj->getAnalysisAdaptor;


=head1 DESCRIPTION
  
  Module to encapsulate all db access for persistent class Analysis.
  There should be just one per application and database connection.
     

=head1 CONTACT

    Contact Arne Stabenau on implemetation/design detail: stabenau@ebi.ac.uk
    Contact Ewan Birney on EnsEMBL in general: birney@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DBSQL::AnalysisAdaptor;

use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;

use vars qw(@ISA);
use strict;

@ISA = qw( Bio::EnsEMBL::DBSQL::BaseAdaptor);



=head2 new

  Args       : Bio::EnsEMBL::DBSQL::DBAdaptor
  Example    : my $aa = new Bio::EnsEMBL::DBSQL::AnalysisAdaptor();
  Description: Creates a new Bio::EnsEMBL::DBSQL::AnalysisAdaptor object and
               internally loads and caches all the Analysis objects from the 
               database.
  Returntype : Bio::EnsEMBL::DBSQL::AnalysisAdaptor
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBSQL::DBAdaptor

=cut

sub new {
  my ($class, $db) = @_;

  my $self = $class->SUPER::new($db);
  
  #load and cache all of the Analysis objects
  $self->fetch_all;

  return $self;
}


=head2 fetch_all

  Args       : none
  Example    : my @analysis = $analysis_adaptor->fetch_all()
  Description: fetches all of the Analysis objects from the database and caches
               them internally.
  Returntype : listref of Bio::EnsEMBL::Analysis retrieved from the database
  Exceptions : none
  Caller     : AnalysisAdaptor::new

=cut

sub fetch_all {
  my $self = shift;
  my ( $analysis, $dbID );
  my $rowHashRef;

  $self->{_cache} = {};
  $self->{_logic_name_cache} = {};
  
  my $sth = $self->prepare( q {
    SELECT analysis_id, logic_name,
           program, program_version, program_file,
           db, db_version, db_file,
           module, module_version,
           gff_source, gff_feature,
           created, parameters
    FROM   analysis } );
  $sth->execute;

  while( $rowHashRef = $sth->fetchrow_hashref ) {
    my $analysis = $self->_objFromHashref( $rowHashRef  );

    $self->{_cache}->{$analysis->dbID}                    = $analysis;
    $self->{_logic_name_cache}->{lc($analysis->logic_name())} = $analysis;
  }

  my @ana = values %{$self->{_cache}};

  return \@ana;
}

sub deleteObj {
    my( $self ) = @_;
    
    $self->{_cache} = undef;
    $self->{_logic_name_cache} = undef;
    $self->SUPER::deleteObj;
}

=head2 fetch_by_dbID

  Arg [1]    : int $internal_analysis_id - the database id of the analysis 
               record to retrieve
  Example    : my $analysis = $analysis_adaptor->fetch_by_dbID(1);
  Description: Retrieves an Analysis object from the database via its internal
               id.
  Returntype : Bio::EnsEMBL::Analysis
  Exceptions : none
  Caller     : general

=cut

sub fetch_by_dbID {
  my $self = shift;
  my $id = shift;

  if( defined $self->{_cache}->{$id} ) {
    return $self->{_cache}->{$id};
  }


  my $query = q{
    SELECT analysis_id, logic_name,
           program, program_version, program_file,
           db, db_version, db_file,
           module, module_version,
           gff_source, gff_feature,
           created, parameters
    FROM   analysis
    WHERE  analysis_id = ? };

  my $sth = $self->prepare($query);  
  $sth->execute( $id );
  my $rowHashRef = $sth->fetchrow_hashref;
  if( ! defined $rowHashRef ) {
    return undef;
  }

  my $anal = $self->_objFromHashref( $rowHashRef );
  $self->{_cache}->{$anal->dbID} = $anal;
  $self->{_logic_name_cache}->{lc($anal->logic_name())} = $anal;
  return $anal;
}


=head2 fetch_by_logic_name

  Arg [1]    : string $logic_name the logic name of the analysis to retrieve
  Example    : my $analysis = $a_adaptor->fetch_by_logic_name('Eponine');
  Description: Retrieves an analysis object from the database using its unique
               logic name.
  Returntype : Bio::EnsEMBL::Analysis
  Exceptions : none
  Caller     : general

=cut

sub fetch_by_logic_name {
  my $self = shift;
  my $logic_name = shift;
  my $analysis;
  my $rowHash;

  #check the cache for the logic name
  if(defined $self->{_logic_name_cache}{lc($logic_name)}) {
    return $self->{_logic_name_cache}{lc($logic_name)};
  }

  my $sth = $self->prepare( "
    SELECT analysis_id, logic_name,
           program, program_version, program_file,
           db, db_version, db_file,
           module, module_version,
           gff_source, gff_feature,
           created, parameters
    FROM   analysis
    WHERE  logic_name = ?" );
  
  $sth->execute($logic_name);
  my $rowHashRef;
  $rowHashRef = $sth->fetchrow_hashref; 

  unless(defined $rowHashRef) {
    return undef;
  }

  $analysis = $self->_objFromHashref( $rowHashRef );
  
  #place the analysis in the caches, cross referenced by dbID and logic_name
  $self->{_cache}{$analysis->dbID()} = $analysis;
  $self->{_logic_name_cache}{lc($logic_name)} = $analysis;

  return $analysis;
}



=head2 store

  Arg [1]    : Bio:EnsEMBL::Analysis $analysis 
  Example    : $analysis_adaptor->store($analysis);
  Description: stores $analysis in db. Does not if already equiped with dbID.
               Sets created date if not already set. Sets dbID and adaptor
               inside $analysis. Returns dbID.
  Returntype : int dbID of stored analysis
  Exceptions : thrown if analysis argument does not have a logic name
  Caller     : ?

=cut

sub store {

  my $self = shift;
  my $analysis = shift;
  
  if( !defined $analysis || !ref $analysis) {
    $self->throw("called store on AnalysisAdaptor with a [$analysis]");
  }

  $analysis->dbID && $analysis->adaptor && ( $analysis->adaptor() == $self ) && 
    return $analysis->dbID;


  my $dbID;

  if( $dbID = $self->exists( $analysis )) {
    $analysis->adaptor( $self );
    $analysis->dbID( $dbID );
    return $dbID;
  }
 
  if( !defined $analysis->logic_name ) {
    $self->throw("Must have a logic name on the analysis object");
  }

 
  if($analysis->created ) {
    my $sth = $self->prepare( q{
      INSERT INTO analysis
      SET created = ?,
          logic_name = ?,
	  db = ?,
	  db_version = ?,
          db_file = ?,
          program = ?,
          program_version = ?,
          program_file = ?,
	  parameters = ?,
          module = ?,
          module_version = ?,
          gff_source = ?,
          gff_feature = ? } );
    $sth->execute
      ( $analysis->created,
	$analysis->logic_name,
	$analysis->db,
	$analysis->db_version,
	$analysis->db_file,
	$analysis->program,
	$analysis->program_version,
	$analysis->program_file,
	$analysis->parameters,
	$analysis->module,
	$analysis->module_version,
	$analysis->gff_source,
	$analysis->gff_feature
      );
    $dbID = $sth->{'mysql_insertid'};
  } else {
    my $sth = $self->prepare( q{

      INSERT INTO analysis
      SET created = now(),
          logic_name = ?,
	  db = ?,
	  db_version = ?,
          db_file = ?,
          program = ?,
          program_version = ?,
          program_file = ?,
	  parameters = ?,
          module = ?,
          module_version = ?,
          gff_source = ?,
          gff_feature = ? } );

    $sth->execute
      ( $analysis->logic_name,
	$analysis->db,
	$analysis->db_version,
	$analysis->db_file,
	$analysis->program,
	$analysis->program_version,
	$analysis->program_file,
	$analysis->parameters,
	$analysis->module,
	$analysis->module_version,
	$analysis->gff_source,
	$analysis->gff_feature
      );

    $dbID = $sth->{'mysql_insertid'};

    if( $dbID ) {
      $sth = $self->prepare( q{
	SELECT created 
	FROM   analysis
	WHERE  analysis_id = ? } );
      $sth->execute( $dbID );
      $analysis->created( ($sth->fetchrow_array)[0] );
    }
  }
  $self->{_cache}->{$dbID} = $analysis;
  $self->{_logic_name_cache}{lc($analysis->logic_name)} = $analysis;

  $analysis->adaptor( $self );
  $analysis->dbID( $dbID );
  
  return $dbID;
}



=head2 update

  Arg [1]    : Bio::EnsEMBL::Analysis $anal
  Example    : $adaptor->update($anal)
  Description: Updates this analysis in the database
  Returntype :
  Exceptions : thrown if $anal arg is not an analysis object
  Caller     : ?

=cut

sub update {
  my ($self, $analysis) = @_;
  
  if (!defined $analysis || !ref $analysis) {
    $self->throw("called update on AnalysisAdaptor with a [$analysis]");
  }

  $analysis->dbID && ($analysis->adaptor() == $self) or
    return undef;

  my $dbID;

  unless ($dbID = $self->exists($analysis)) {
    return undef;
  }

  my $query = "UPDATE analysis SET ";

  foreach my $m (qw/
    created         logic_name
    db              db_version      db_file
    program         program_version program_file
    parameters      module          module_version
    gff_source      gff_feature/) {
    $query .= " $m = '" . $analysis->$m . "'," if defined $analysis->$m;
  };
  chop $query;
  $query .= " WHERE analysis_id = $dbID";

  my $sth = $self->db->db_handle->do($query);

  return 1;
}



=head2 remove

  Arg [1]    : Bio::EnsEMBL::Analysis $anal
  Example    : $adaptor->remove($anal)
  Description: Removes this analysis from the database
  Returntype :
  Exceptions : thrown if $anal arg is not an analysis object
  Caller     : ?

=cut

sub remove {
  my ($self, $analysis) = @_;
  my $dbID;
  
  if (!defined $analysis || !ref $analysis) {
    $self->throw("called remove on AnalysisAdaptor with a [$analysis]");
  }

  unless ($dbID = $self->exists($analysis)) {
    return undef;
  }

  my $res = $self->db->db_handle->do(qq{
    DELETE from analysis
    WHERE  analysis_id = $dbID
  });
}



=head2 exists

  Arg [1]    : Bio::EnsEMBL::Analysis $anal
  Example    : if($adaptor->exists($anal)) #do something
  Description: Tests whether this Analysis already exists in the database.
               Returns the dbID, if it does, undef if it doesnt.
  Returntype : int or undef
  Exceptions : thrown if $anal arg is not an analysis object
  Caller     : ?

=cut

sub exists {
  my ($self,$anal) = @_;

  unless($anal->isa("Bio::EnsEMBL::Analysis")) {
    $self->throw("Object is not a Bio::EnsEMBL::Analysis");
  } 
  
  
  # objects with already have this adaptor are store here.
  if( $anal->can("adaptor") && defined $anal->adaptor &&
      $anal->adaptor == $self ) {
    if (my $id = $anal->dbID) {
      return $id;
    }
    else {
      $self->throw ("analysis does not have an analysisId");
    }
  }
  
  foreach my $cacheId (keys %{$self->{_cache}}) {
    if ($self->{_cache}->{$cacheId}->compare($anal) >= 0) {
      # $anal->dbID( $cacheId );
      # $anal->adaptor( $self );
      return $cacheId;
    }
  }
  return undef;
}


=head2 _objFromHashref

  Arg [1]    : hashref $rowHash
  Description: Private helper function generates an Analysis object from a 
               mysql row hash reference.
  Returntype : Bio::EnsEMBL::Analysis
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBSQL::AnalsisAdaptor::fetch_* methods

=cut
  
sub _objFromHashref {
  my $self = shift;
  my $rowHash = shift;

  my $analysis = Bio::EnsEMBL::Analysis->new(
      -id              => $rowHash->{analysis_id},
      -adaptor         => $self,
      -db              => $rowHash->{db},
      -db_file         => $rowHash->{db_file},
      -db_version      => $rowHash->{db_version},
      -program         => $rowHash->{program},
      -program_version => $rowHash->{program_version},
      -program_file    => $rowHash->{program_file},
      -gff_source      => $rowHash->{gff_source},
      -gff_feature     => $rowHash->{gff_feature},
      -module          => $rowHash->{module},
      -module_version  => $rowHash->{module_version},
      -parameters      => $rowHash->{parameters},
      -created         => $rowHash->{created},
      -logic_name      => $rowHash->{logic_name}
    );
  
  return $analysis;
}


=head2 db

  Arg [1]    : (optional) Bio::EnsEMBL::DBSQL::DBAdaptor $db
               the database used by this adaptor.
  Example    : my $db = $analysis_adaptor->db()
  Description: Getter/Setter for the database this adaptor uses internally
               to fetch and store database objects.
  Returntype : Bio::EnsEMBL::DBSQL::DBAdaptor
  Exceptions : none
  Caller     : BaseAdaptor::new, general

=cut

sub db {
  my ( $self, $arg )  = @_;
  ( defined $arg ) &&
    ($self->{_db} = $arg);
  $self->{_db};;
}


1;
