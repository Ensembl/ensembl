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
use Time::Local;

use vars qw(@ISA);
use strict;

@ISA = qw( Bio::EnsEMBL::DBSQL::BaseAdaptor Bio::Root::RootI );

sub new {
  my $class = shift;
  my $self = bless {},$class;
  
  my $dbobj = shift;

  $self->db( $dbobj );
  $self->fetch_all;
  return $self;
}

=head2 fetch_all

  Title   : fetch_all
  Usage   : @analyses = $self->fetch_all;
  Function: retrieves all analyses from db;
  Returns : List of Bio::EnsEMBL::Analysis
  Args    : -

=cut

sub fetch_all {
  my $self = shift;
  my ( $analysis, $dbID );
  my $rowHashRef;

  $self->{_cache} = {};
  
  my $sth = $self->prepare( q {
    SELECT analysisId, logic_name,
           program,program_version,program_file,
           db,db_version,db_file,
           module,module_version,
           gff_source,gff_feature,
           created, parameters
    FROM analysisprocess } );
  $sth->execute;

  while( $rowHashRef = $sth->fetchrow_hashref ) {
    my $analysis = $self->_objFromHashref( $rowHashRef  );
    $self->{_cache}->{$analysis->dbID} = $analysis;
  }

  return values %{$self->{_cache}};
}

=head2 fetch_by_dbID

  Title   : fetch_by_dbID
  Usage   : my $analysis = $adaptor->fetch_by_dbID
  Function: Retrieves an analysis from database by internal id
  Returns : throws exception when something goes wrong.
            undef if the id is not in the db.
  Args    : 

=cut

sub fetch_by_dbID {
  my $self = shift;
  my $id = shift;

  if( defined $self->{_cache}->{$id} ) {
    return $self->{_cache}->{$id};
  }

  my $sth = $self->prepare( q{
    SELECT analysisId, logic_name,
           program,program_version,program_file,
           db,db_version,db_file,
           module,module_version,
           gff_source,gff_feature,
           created, parameters
    FROM analysisprocess
    WHERE analysisId = ? } );
  
  $sth->execute( $id );
  my $rowHashRef = $sth->fetchrow_hashref;
  if( ! defined $rowHashRef ) {
    return undef;
  }

  my $anal = $self->_objFromHashref( $rowHashRef );
  $self->{_cache}->{$anal->dbID} = $anal;
  return $anal;
}

sub fetch_by_newest_logic_name {
  my $self = shift;
  my $logic_name = shift;

  my $sth = $self->prepare( q{
    SELECT analysisId, logic_name,
           program,program_version,program_file,
           db,db_version,db_file,
           module,module_version,
           gff_source,gff_feature,
           created, parameters
    FROM analysisprocess
    WHERE logic_name = ?
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
    SELECT analysisId, logic_name,
           program,program_version,program_file,
           db,db_version,db_file,
           module,module_version,
           gff_source,gff_feature,
           created, parameters
    FROM analysisprocess
    WHERE logic_name = ?
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

# store makes dbID for analysis object
# sets the creation time in created if it wasnt set before
sub store {

  my $self = shift;
  my $analysis = shift;
  $analysis->dbID && return $analysis->dbID;
  my $dbID;
  
  if( defined $analysis->created ) {
    my $sth = $self->prepare( q{
      INSERT INTO analysisprocess
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

      INSERT INTO analysisprocess
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
	FROM analysisprocess
	WHERE analysisId = ? } );
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
 Returns : true if this or more detailed exists.
 Args    : Bio::EnsEMBL::Analysis

=cut

sub exists {
    my ($self,$anal) = @_;

    
    $self->throw("Object is not a Bio::EnsEMBL::Analysis") unless $anal->isa("Bio::EnsEMBL::Analysis");
    
    # objects with already have this adaptor are store here.
    if( $anal->can("adaptor") && defined $anal->adaptor &&
      $anal->adaptor == $self ) {
      return 1;
    }

    while(  my ( $cacheId, $cacheAna) = each %{$self->{_cache}} ) {
      if( $cacheAna->compare( $anal ) >= 0 ) {
        return 1;
      }
    }
    return 0;
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

  my $analysis = Bio::EnsEMBL::Analysis->new
    ( -id => $rowHash->{analysisId},
      -db => $rowHash->{db},
      -db_file => $rowHash->{db_file},
      -program => $rowHash->{program},
      -program_version => $rowHash->{program_version},
      -program_file => $rowHash->{program_file},
      -gff_source => $rowHash->{gff_source},
      -gff_feature => $rowHash->{gff_feature},
      -module => $rowHash->{module},
      -module_version => $rowHash->{module_version},
      -parameters => $rowHash->{parameters},
      -created => $rowHash->{created},
      -logic_name => $rowHash->{logic_name}
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

  my $sth = $self->prepare( "drop table if exists analysisprocess" );
  $sth->execute();

  $sth = $self->prepare( qq{
    CREATE TABLE analysisprocess (
      analysisId int(10) unsigned DEFAULT '0' NOT NULL auto_increment,
      created datetime DEFAULT '0000-00-00 00:00:00' NOT NULL,
      logic_name varchar(40) not null,
      db varchar(40),
      db_version varchar(40),
      db_file varchar(80),
      program varchar(40),
      program_version varchar(40),
      program_file varchar(40),
      parameters varchar(80),
      module varchar(80),
      module_version varchar(40),
      gff_source varchar(40),
      gff_feature varchar(40),
      PRIMARY KEY (analysisId)
    )
  } );
  $sth->execute();
}

1;


