#
# EnsEMBL module for DBEntry
#
# Cared for by Arne Stabenau <stabenau@ebi.ac.uk>
#
# Copyright EMBL/EBI 2001
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::DBEntry - Module to collect information about an
 external reference

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 CONTACT

  Arne Stabenau <stabenau@ebi.ac.uk>
  Ewan Birney <birney@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DBEntry;

use Bio::Annotation::DBLink;
use vars qw(@ISA $AUTOLOAD);
use strict;


@ISA = qw( Bio::EnsEMBL::Root Bio::Annotation::DBLink );


sub new_fast {
  my $class = shift;
  my $hashref = shift;

  return bless $hashref, $class;
}

sub new {
  my ($class, @args) = @_;

  my $self = bless {},$class;

  my ( $adaptor, $dbID, $primary_id, $version,
       $dbname, $release, $display_id  ) = $self->_rearrange
	 ( [ qw { ADAPTOR
		DBID
		PRIMARY_ID
		VERSION
		DBNAME
		RELEASE
		DISPLAY_ID
	      }], @args );

  if( defined $adaptor ) { $self->adaptor( $adaptor )}
  if( defined $dbID ) { $self->dbID( $dbID ) }
  if( defined $primary_id ) { $self->primary_id( $primary_id ) }
  if( defined $version ) { $self->version( $version ) }
  if( defined $dbname ) { $self->dbname( $dbname ) }
  if( defined $release) { $self->release( $release ) }
  if( defined $display_id) { $self->display_id( $display_id ) }
  $self->{_synonyms} = [];;

  return $self;
}


=head2 primary_id

  Arg [1]    : string $primary_id
  Example    : none
  Description: get/set for attribute primary_id
               its this objects primary id in the external database
  Returntype : string
  Exceptions : none
  Caller     : general

=cut


sub primary_id {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    $self->{_primary_id} = $arg;
  } 
  return $self->{_primary_id};
}


=head2 display_id

  Arg [1]    : string $display_id
  Example    : none
  Description: get/set for attribute display_id
               This objects preferred display string. That can be the same
               as primary id or ensembl specific.
  Returntype : string
  Exceptions : none
  Caller     : general

=cut



sub display_id{
   my ( $self, $arg ) = @_;
   if( defined $arg ) {
       $self->{_display_id} = $arg;
   } 
   return $self->{_display_id};

}


=head2 dbname

  Arg [1]    : string $dbname
  Example    : none
  Description: get/set for attribute dbname
  Returntype : string
  Exceptions : none
  Caller     : general

=cut


sub dbname {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    $self->{_dbname} = $arg;
  } 
  return $self->{_dbname};
}


=head2 database

  Args       : none
  Example    : none
  Description: additional get for the dbname to make it compliant to
               some interface
  Returntype : string
  Exceptions : none
  Caller     : general

=cut


sub database {
  my $self = shift;
  return $self->dbname();
}


=head2 optional_id

  Args       : none
  Example    : none
  Description: additional get for the display_id to make it compliant to
               some interface
  Returntype : string
  Exceptions : none
  Caller     : general

=cut



sub optional_id {
  my $self = shift;
  return $self->display_id;
}


=head2 release

  Arg [1]    : string $release
  Example    : none
  Description: get/set for attribute release
  Returntype : string
  Exceptions : none
  Caller     : general

=cut



sub release {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    $self->{_release} = $arg;
  } 
  return $self->{_release};
}


sub adaptor {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    $self->{_adaptor} = $arg;
  } 
  return $self->{_adaptor};
}


sub dbID {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    $self->{_dbID} = $arg;
  } 
  return $self->{_dbID};
}



=head2 version

  Arg [1]    : string $version
  Example    : none
  Description: get/set for attribute version
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub version {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    $self->{_version} = $arg;
  } 
  return $self->{_version};
}




=head2 description

  Arg [1]    : string $description
  Example    : none
  Description: get/set for attribute description
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub description {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    $self->{_description} = $arg;
  } 
  return $self->{_description};
}


=head2 add_synonym

  Arg  1     : string $synonym
  Example    : none
  Description: adding a synonynm for the external object under which it is 
               also known
  Returntype : none
  Exceptions : none
  Caller     : general

=cut


sub add_synonym {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    push( @{$self->{_synonyms}}, $arg );
  }
}


=head2 get_synonyms

  Args       : none
  Example    : @synonyms = @{$db_entry->get_all_synonyms()};
  Description: get a list of synonym added to this object
  Returntype : list reference of strings 
  Exceptions : none
  Caller     : general

=cut

sub get_all_synonyms {
  my $self = shift;
  return $self->{_synonyms};
}


=head2 flush_synonyms

  Args       : none
  Example    : none
  Description: remove all synonyms from this object
  Returntype : none
  Exceptions : none
  Caller     : general

=cut

sub flush_synonyms {
  my $self = shift;
  $self->{_synonyms} = [];
}


=head2 status

  Arg [1]    : string $status
  Example    : none
  Description: get/set for attribute status
  Returntype : string
  Exceptions : none
  Caller     : general

=cut


sub status{
 my ( $self, $arg ) = @_;
   if( defined $arg ) {
       $self->{_status} = $arg;
   } 
   return $self->{_status};
}



=head2 comment

  Args       : none
  Example    : none
  Description: additional get for description to comply with bioperl
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

#Cheat to comply with bioperl
sub comment {
    my ($self) = @_;
    if ($self) {
	return $self->description();
    }
}



=head2 get_synonyms

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use get_all_synonyms instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_synonyms {
  my $self = shift;

  $self->warn("get_synonyms has been rename get_all_synonyms\n" . caller);
  return $self->get_all_synonyms;
}

1;
