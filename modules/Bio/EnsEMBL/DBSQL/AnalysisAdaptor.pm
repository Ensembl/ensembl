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

  $analysisAdaptor = $dbobj->getAnalysisAdaptor;
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
use Time::Local;

use vars qw(@ISA);
use strict;

@ISA = qw( Bio::EnsEMBL::DBSQL::BaseAdaptor Bio::Root::RootI );


=head2 new

  Arg  1    : Bio::EnsEMBL::DBSQL::DBAdaptor $dbadaptor
  Function  : create an AnalysisAdaptor. Caches all Analysis objects from the database.
  Returntype: Bio::EnsEMBL::DBSQL::AnalysisAdaptor
  Exceptions: none
  Caller    : DBAdaptor::get_AnalysisAdaptor

=cut


sub new {
  my $class = shift;
  my $self = bless {},$class;
  
  my $dbobj = shift;

  $self->db( $dbobj );
  $self->fetch_all;
  return $self;
}


=head2 fetch_all

  Args      : none
  Function  : Retrieves all Analysis objects from the database,
              caches them.
  Returntype: list Bio::EnsEMBL::Analysis
  Exceptions: none
  Caller    : $self->new

=cut


sub fetch_all {
  my $self = shift;
  my ( $analysis, $dbID );
  my $rowHashRef;

  $self->{_cache} = {};
  
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
    $self->{_cache}->{$analysis->dbID} = $analysis;
  }

  return values %{$self->{_cache}};
}


=head2 fetch_by_dbID

  Arg 1     : int $internal_analysis_id
  Function  : Retrieves an analysis from database by internal id.
              Returns undef if id not present in db.
  Returntype: Bio::EnsEMBL::Analysis
  Exceptions: none
  Caller    : generally used

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
  return $anal;
}


=head2 fetch_by_newest_logic_name

  Arg  1    : txt $logic_name
  Function  : Retrieve latest Analysis object from db with given logic_name
  Returntype: Bio::EnsEMBL::Analysis
  Exceptions: none
  Caller    : Probably Pipeline control scripts

=cut


sub fetch_by_newest_logic_name {
  my $self = shift;
  my $logic_name = shift;

  my $sth = $self->prepare( q{
    SELECT analysis_id, logic_name,
           program, program_version, program_file,
           db, db_version, db_file,
           module, module_version,
           gff_source, gff_feature,
           created, parameters
    FROM   analysis
    WHERE  logic_name = ?
    ORDER BY created DESC } );
  
  $sth->execute( $logic_name );
  my $rowHashRef = $sth->fetchrow_hashref;
  if( ! defined $rowHashRef ) {
    return undef;
  }

  return $self->_objFromHashref( $rowHashRef );
}


sub fetch_by_logic_name {
  my $self = shift;
  my $logic_name = shift;
  my @result;
  my $analysis;
  my $rowHash;

  my $sth = $self->prepare( q{
    SELECT analysis_id, logic_name,
           program, program_version, program_file,
           db, db_version, db_file,
           module, module_version,
           gff_source, gff_feature,
           created, parameters
    FROM   analysis
    WHERE  logic_name = ?
    ORDER BY created DESC } );
  
  $sth->execute( $logic_name );
  my $rowHashRef;
  while( $rowHashRef = $sth->fetchrow_hashref ) {
       $analysis = $self->_objFromHashref( $rowHashRef );
    if( defined $analysis ) {
      push( @result, $analysis );
    }
  }
  return @result;
}


=head2 store

  Arg  1    : Bio:EnsEMBL::Analysis $analysis
  Function  : stores $analysis in db. Doesn if already equppied with dbID
              Sets created date if not already set. Sets dbID and adaptor
              inside $analysis. Returns dbID.
  Returntype: int
  Exceptions: none
  Caller    : Every store that links to analysis object

=cut


sub store {

  my $self = shift;
  my $analysis = shift;
  
  if( !defined $analysis || !ref $analysis) {
    $self->throw("called store on AnalysisAdaptor with a [$analysis]");
  }

  $analysis->dbID && return $analysis->dbID;
  my $dbID;
 
  if( !defined $analysis->logic_name ) {
    $self->throw("Must have a logic name on the analysis object");
  }

 
  if( defined $analysis->created ) {
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

    if( defined $dbID ) {
      $sth = $self->prepare( q{
	SELECT created 
	FROM   analysis
	WHERE  analysis_id = ? } );
      $sth->execute( $dbID );
      $analysis->created( ($sth->fetchrow_array)[0] );
    }
  }
  $self->{_cache}->{$dbID} = $analysis;

  if( $analysis->can( "adaptor" )) {
    $analysis->adaptor( $self );
    $analysis->dbID( $dbID );
  }
  
  return $dbID;
}

=head2 exists

 Title   : exists
 Usage   : $adaptor->exists($anal)
 Function: Tests whether this Analysis already exists in the database
 Example :
 Returns : analysisId or undef
 Args    : Bio::EnsEMBL::Analysis

=cut

sub exists {
    my ($self,$anal) = @_;

    $self->throw("Object is not a Bio::EnsEMBL::Analysis") unless $anal->isa("Bio::EnsEMBL::Analysis");
    
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

=head2 mysql2Unixtime

  Title    : mysql2Unixtime
  Usage    : no object function, yet.
  Function : Calculates the unix time from mysql time string
  Example  :
  Returns  : a unix time
  Args     : "2001-01-09 18:02:13" somthing like this

=cut


sub mysql2Unixtime {
  my $sqltime = shift;
  
  my ($year,$month,$mday,$hour,$min,$sec ) = ( $sqltime =~ /(\d+)-(\d+)-(\d+)\s+(\d+):(\d+):(\d+)/ );
  my $time = timelocal( $sec, $min, $hour, $mday, $month, $year );
}

  
sub _objFromHashref {
  my $self = shift;
  my $rowHash = shift;

  my $analysis = Bio::EnsEMBL::Analysis->new(
      -id              => $rowHash->{analysis_id},
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

sub db {
  my ( $self, $arg )  = @_;
  ( defined $arg ) &&
    ($self->{_db} = $arg);
  $self->{_db};;
}

sub prepare {
  my ( $self, $query ) = @_;
  $self->db->prepare( $query );
}


sub deleteObj {
  my $self = shift;
  my @dummy = values %{$self};
  foreach my $key ( keys %$self ) {
    delete $self->{$key};
  }
  foreach my $obj ( @dummy ) {
    eval {
      $obj->deleteObj;
    }
  }
}



sub create_tables {
  my $self = shift;

  my $sth = $self->prepare( "drop table if exists analysis" );
  $sth->execute();

  $sth = $self->prepare( qq{
    CREATE TABLE analysis (
      analysis_id int(10) unsigned DEFAULT '0' NOT NULL auto_increment,
      created datetime DEFAULT '0000-00-00 00:00:00' NOT NULL,
      logic_name varchar(40) not null,
      db varchar(120),
      db_version varchar(40),
      db_file varchar(120),
      program varchar(80),
      program_version varchar(40),
      program_file varchar(40),
      parameters varchar(80),
      module varchar(80),
      module_version varchar(40),
      gff_source varchar(40),
      gff_feature varchar(40),
      PRIMARY KEY (analysis_id)
    )
  } );
  $sth->execute();
}


1;
