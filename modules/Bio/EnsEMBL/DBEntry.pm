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

Bio::EnsEMBL::DBEntry - Module to collect information about an external reference

=head1 SYNOPSIS

=head1 DESCRIPTION

This module stores information about external references to EnsEMBL objects

=head1 CONTACT

Post questions to the EnsEMBL developer mailing list: <ensembl-dev@ebi.ac.uk> 

=head1 METHODS

=cut



package Bio::EnsEMBL::DBEntry;

use Bio::EnsEMBL::Storable;
use Bio::Annotation::DBLink;

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(deprecate);

use vars qw(@ISA $AUTOLOAD);
use strict;


@ISA = qw( Bio::EnsEMBL::Storable Bio::Annotation::DBLink );



=head2 new_fast

  Arg [1]    : hash reference $hashref
  Example    : none
  Description: A very quick constructor that requires internal knowledge of
               the class. This is used in speed critical sections of the code
               where many objects need to be created quickly.
  Returntype : Bio::EnsEMBL::DBEntry
  Exceptions : none
  Caller     : ?
  Status     : Stable

=cut

sub new_fast {
  my $class = shift;
  my $hashref = shift;
  bless $hashref, $class;
  $hashref->{synonyms} = [];
  
  return $hashref;
}



=head2 new

  Args [...] : list of named parameters 
  Example    : my $dbentry = new Bio::EnsEMBL::DBEntry(-adaptor => $adaptor,
						       -primary_id => $pid,
						       -version => $version,
						       -dbname  => $dbname,
						       -release => $release,
						       -display_id => $did,
                                                       -description => $description,
						       -primary_id_linkable =>$primary_id_linkable,
						       -display_id_linkable =>$display_id_linkable,
						       -priority => $priority,
						       -db_display_name => $db_display_name,
						       -info_type => $info_type,
						       -info_text => $info_text);
  Description: Creates a new DBEntry object
  Returntype : Bio::EnsEMBL::DBEntry
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBEntryAdaptor
  Status     : At Risk
             : due to 'PRIMARY_ID_LINKABLE','DISPLAY_ID_LINKABLE','PRIORITY'
             : being under development
             : If you don't use any of these the sub can be considered Stable

=cut

sub new {
  my ($class, @args) = @_;

  my $self = bless {},$class;

  my ( $adaptor, $dbID, $primary_id, $version,
       $dbname, $release, $display_id, $description,
       $primary_id_linkable, $display_id_linkable, $priority, 
       $db_display_name, $info_type, $info_text) =
    rearrange ( ['ADAPTOR','DBID','PRIMARY_ID','VERSION',
                 'DBNAME','RELEASE','DISPLAY_ID','DESCRIPTION',
		 'PRIMARY_ID_LINKABLE','DISPLAY_ID_LINKABLE','PRIORITY',
		 'DB_DISPLAY_NAME', 'INFO_TYPE', 'INFO_TEXT'], @args );

  $self->{'adaptor'} = $adaptor;
  $self->{'dbID'}    = $dbID;

  if( defined $primary_id ) { $self->primary_id( $primary_id ) }
  if( defined $version ) { $self->version( $version ) } else
    { $self->version( "" ); }
  if( defined $dbname ) { $self->dbname( $dbname ) }
  if( defined $release) { $self->release( $release ) }
  if( defined $display_id) { $self->display_id( $display_id ) }
  if( defined $description) { $self->description($description) }
  if( defined $primary_id_linkable) { $self->primary_id_linkable($primary_id_linkable) }
  if( defined $display_id_linkable) { $self->display_id_linkable($display_id_linkable) }
  if( defined $priority) { $self->priority($priority) }
  if( defined $db_display_name) { $self->db_display_name($db_display_name) }
  if( defined $info_type) { $self->db_display_name($info_type) }
  if( defined $info_text) { $self->db_display_name($info_text) }
  $self->{synonyms} = [];;

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
  Status     : Stable

=cut

sub primary_id {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    $self->{primary_id} = $arg;
  } 

  return $self->{primary_id};
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
  Status     : Stable

=cut

sub display_id{
   my ( $self, $arg ) = @_;
   if( defined $arg ) {
       $self->{display_id} = $arg;
   } 
   return $self->{display_id};

}



=head2 dbname

  Arg [1]    : string $dbname
  Example    : none
  Description: get/set for attribute dbname
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub dbname {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    $self->{dbname} = $arg;
  } 
  return $self->{dbname};
}



=head2 database

  Args       : none
  Example    : none
  Description: additional get for the dbname to make it compliant to
               some interface
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

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
  Status     : Stable

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
  Status     : Stable

=cut

sub release {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    $self->{release} = $arg;
  } 
  return $self->{release};
}


=head2 version

  Arg [1]    : string $version
  Example    : none
  Description: get/set for attribute version
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub version {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    $self->{version} = $arg;
  } 
  return $self->{version};
}




=head2 description

  Arg [1]    : string $description
  Example    : none
  Description: get/set for attribute description
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub description {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    $self->{description} = $arg;
  } 
  return $self->{description};
}


=head2 primary_id_linkable

  Arg [1]    : boolean $primary_id_linkable
  Example    : none
  Description: get/set for attribute primary_id_linkable
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk
             : due to it being under development

=cut

sub primary_id_linkable {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    $self->{primary_id_linkable} = $arg;
  } 
  return $self->{primary_id_linkable};
}

=head2 display_id_linkable

  Arg [1]    : boolean $display_id_linkable
  Example    : none
  Display_id_linkable: get/set for attribute display_id_linkable
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk
             : due to it being under development

=cut

sub display_id_linkable {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    $self->{display_id_linkable} = $arg;
  } 
  return $self->{display_id_linkable};
}


=head2 priority

  Arg [1]    : int $priority
  Example    : none
  Priority   : get/set for attribute priority
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk
             : due to it being under development

=cut

sub priority {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    $self->{priority} = $arg;
  } 
  return $self->{priority};
}


=head2 db_display_name

  Arg [1]    : String $db_display_name
  Example    : none
  Db_display_name: get/set for attribute db_display_name
  Returntype : string; has "Projected " prepended if info_type='PROJECTION'
  Exceptions : none
  Caller     : general

=cut

sub db_display_name {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    $self->{db_display_name} = $arg;
  }

  my $name;
  if ($self->{info_type} eq "PROJECTION") {
    $name = "Projected " . $self->{db_display_name};
  } else {
    $name =  $self->{db_display_name};
  }

  return $name;
}


=head2 info_type

  Arg [1]    : String $info_type
  Example    : none
  Info_type: get/set for attribute info_type
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub info_type {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    $self->{info_type} = $arg;
  } 
  return $self->{info_type};
}


=head2 info_text

  Arg [1]    : String $info_text
  Example    : none
  Info_text: get/set for attribute info_text
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub info_text {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    $self->{info_text} = $arg;
  } 
  return $self->{info_text};
}


=head2 add_synonym

  Arg  1     : string $synonym
  Example    : none
  Description: adding a synonynm for the external object under which it is 
               also known
  Returntype : none
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut


sub add_synonym {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    push( @{$self->{synonyms}}, $arg );
  }
}


=head2 get_all_synonyms

  Args       : none
  Example    : @synonyms = @{$db_entry->get_all_synonyms()};
  Description: get a list of synonym added to this object
  Returntype : list reference of strings 
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_synonyms {
  my $self = shift;
  return $self->{synonyms};
}


=head2 flush_synonyms

  Args       : none
  Example    : none
  Description: remove all synonyms from this object
  Returntype : none
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub flush_synonyms {
  my $self = shift;
  $self->{synonyms} = [];
}


=head2 status

  Arg [1]    : string $status
  Example    : none
  Description: get/set for attribute status
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut


sub status{
 my ( $self, $arg ) = @_;
   if( defined $arg ) {
       $self->{status} = $arg;
   } 
   return $self->{status};
}



=head2 comment

  Args       : none
  Example    : none
  Description: additional get for description to comply with bioperl
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

#Cheat to comply with bioperl
sub comment {
    my ($self) = @_;
    if ($self) {
	return $self->description();
    }
}


=head1 DEPRECATED METHODS

=cut

=head2 get_synonyms

  Description: DEPRECATED use get_all_synonyms instead

=cut

sub get_synonyms {
  my $self = shift;

  deprecate("get_synonyms has been renamed get_all_synonyms.");
  return $self->get_all_synonyms;
}

1;
