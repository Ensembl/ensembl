package Bio::EnsEMBL::DBEntry;

=head1 NAME

Bio::EnsEMBL::DBEntry -
Object representing an external reference (xref)

=head1 SYNOPSIS


=head1 DESCRIPTION

This object holds information about external references (xrefs) to Ensembl
objects.

=head1 METHODS


=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http://www.ensembl.org/info/about/code_licence.html for details.

=head1 AUTHOR

Arne Stabenau <stabenau@ebi.ac.uk>, Ensembl core API team

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<ensembl-dev@ebi.ac.uk>

=cut

use strict;
use warnings;
no warnings qw(uninitialized);

use Bio::EnsEMBL::Storable;
use Bio::Annotation::DBLink;

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(deprecate);

our @ISA = qw(Bio::EnsEMBL::Storable Bio::Annotation::DBLink);


=head2 new_fast

  Arg [1]    : Hashref $hashref - hash reference to bless as new DBEntry object
  Description: A very quick constructor that requires internal knowledge of
               the class. This is used in speed critical sections of the code
               where many objects need to be created quickly.
  Returntype : Bio::EnsEMBL::DBEntry
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub new_fast {
  my $class = shift;
  my $hashref = shift;

  bless $hashref, $class;

  return $hashref;
}


=head2 new

  Args [...] : list of named parameters 
  Example    : my $dbentry = new Bio::EnsEMBL::DBEntry(
                    -adaptor => $adaptor,
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
                    -info_text => $info_text,
                    -type => $type,
                    -secondary_db_name => $secondary_db_name,
                    -secondary_db_table => $secondary_db_table
                    -linkage_annotation => $object_xref_text);
  Description: Creates a new DBEntry object
  Returntype : Bio::EnsEMBL::DBEntry
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBEntryAdaptor
  Status     : At Risk
               Due to 'PRIMARY_ID_LINKABLE','DISPLAY_ID_LINKABLE','PRIORITY',
              'INFO_TYPE', 'INFO_TEXT', ''DB_DISPLAY_NAME', 'TYPE',
              'SECONDARY_DB_NAME', 'SECONDARY_DB_TABLE'
               being under development - if you don't use any of these the
               method can be considered Stable

=cut

sub new {
  my ($class, @args) = @_;
  
  my $self = bless {},$class;

  my ( $adaptor, $dbID, $primary_id, $version,
       $dbname, $release, $display_id, $description,
       $primary_id_linkable, $display_id_linkable, $priority,
       $db_display_name, $info_type, $info_text, $type,
       $secondary_db_name, $secondary_db_table, $link_annotation) =
    rearrange ( ['ADAPTOR','DBID','PRIMARY_ID','VERSION',
                 'DBNAME','RELEASE','DISPLAY_ID','DESCRIPTION',
		 'PRIMARY_ID_LINKABLE','DISPLAY_ID_LINKABLE','PRIORITY',
		 'DB_DISPLAY_NAME', 'INFO_TYPE', 'INFO_TEXT', 'TYPE',
                 'SECONDARY_DB_NAME', 'SECONDARY_DB_TABLE', 'LINKAGE_ANNOTATION'], @args );

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
  if( defined $info_type) { $self->info_type($info_type) }
  if( defined $info_text) { $self->info_text($info_text) }
  if( defined $type) { $self->type($type) }
  if( defined $secondary_db_name) { $self->secondary_db_name($secondary_db_name) }
  if( defined $secondary_db_table) { $self->secondary_db_table($secondary_db_table) }
  $self->linkage_annotation($link_annotation) if defined $link_annotation;


  return $self;
}


=head2 primary_id

  Arg [1]    : (optional) String $arg - value to set
  Example    : none
  Description: Getter/setter for attribute 'primary_id'.
               This is the object's primary id in the external database.
  Returntype : String
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

  Arg [1]    : (optional) String $arg - value to set
  Example    : none
  Description: Getter/setter for attribute 'display_id'.
               The object's preferred display name. This can be the same
               as primary_id or ensembl-specific.
  Returntype : String
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


=head2 optional_id

  Args       : none
  Example    : none
  Description: Additional getter for attribute 'display_id'.
               The object's preferred display name.
               Only include for BioPerl interface compliance, please use
               $self->display_id().
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub optional_id {
  my $self = shift;
  return $self->display_id;
}


=head2 dbname

  Arg [1]    : (optional) String $arg - value to set
  Example    : none
  Description: Getter/setter for attribute 'dbname'.
               The name of the external database.
  Returntype : String
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
  Description: Additional getter for attribute 'dbname'.
               The name of the external database.
               Only include for BioPerl interface compliance, please use
               $self->dbname().
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub database {
  my $self = shift;
  return $self->dbname();
}


=head2 release

  Arg [1]    : (optional) String $arg - value to set
  Example    : none
  Description: Getter/setter for attribute 'release'.
               The external database release name.
  Returntype : String
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

  Arg [1]    : (optional) String $arg - value to set
  Example    : none
  Description: Getter/setter for attribute 'version'.
               The object's version in the external database.
  Returntype : String
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

  Arg [1]    : (optional) String $arg - value to set
  Example    : none
  Description: Getter/setter for attribute 'description'.
               The object's description.
  Returntype : String
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


=head2 comment

  Args       : none
  Example    : none
  Description: Additional getter for attribute 'description'.
               The object's description.
               Only include for BioPerl interface compliance, please use
               $self->description().
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub comment {
  my $self = shift;
  return $self->description();
}


=head2 primary_id_linkable

  Arg [1]    : (optional) Boolean $arg - value to set
  Example    : none
  Description: Getter/setter for attribute 'primary_id_linkable'.
  Returntype : String
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

  Arg [1]    : (optional) Boolean $arg - value to set
  Example    : none
  Description: Getter/setter for attribute 'display_id_linkable'.
  Returntype : String
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
  Priority   : Getter/setter for attribute 'priority'. Note this
               is the priority from the external_db table.
  Returntype : String
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
  Description: Getter/setter for attribute 'db_display_name'.
               The preferred display name for the external database. Has
               "Projected " prepended if info_type='PROJECTION'.
  Returntype : String
  Exceptions : none
  Caller     : general

=cut

sub db_display_name {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    $self->{db_display_name} = $arg;
  }

  my $name;
  if ($self->{info_type} && $self->{info_type} eq "PROJECTION") {
    $name = "Projected " . $self->{db_display_name};
  } else {
    $name =  $self->{db_display_name};
  }

  return $name;
}


=head2 info_type

  Arg [1]    : String $info_type
  Example    : none
  Description: Getter/setter for attribute 'info_type'.
  Returntype : String
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
  Description: Getter/setter for attribute 'info_text'.
  Returntype : String
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

=head2 linkage_annotation

  Arg [1]    : String $object_xref_text
  Example    : none
  Description: Getter/setter for attribute 'linkage_annotation'.
  Returntype : String
  Exceptions : none
  Caller     : general

=cut

sub linkage_annotation {
  my ( $self, $arg ) = @_;

  $self->{linkage_annotation} = $arg if defined $arg;
  
  return $self->{linkage_annotation};
}


=head2 type

  Arg [1]    : String $type
  Example    : none
  Description: Getter/setter for attribute 'type'.
  Returntype : String
  Exceptions : none
  Caller     : general

=cut

sub type {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    $self->{type} = $arg;
  }
  return $self->{type};
}

=head2 secondary_db_name

  Arg [1]    : String $secondary_db_name
  Example    : none
  Description: Getter/setter for attribute 'secondary_db_name'.
  Returnsecondary_db_name : String
  Exceptions : none
  Caller     : general

=cut

sub secondary_db_name {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    $self->{secondary_db_name} = $arg;
  }
  return $self->{secondary_db_name};
}


=head2 secondary_db_table

  Arg [1]    : String $secondary_db_table
  Example    : none
  Description: Getter/setter for attribute 'secondary_db_table'.
  Returnsecondary_db_table : String
  Exceptions : none
  Caller     : general

=cut

sub secondary_db_table {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    $self->{secondary_db_table} = $arg;
  }
  return $self->{secondary_db_table};
}


=head2 add_synonym

  Arg [1]    : String $arg - synonym to add
  Example    : none
  Description: Add a synonym for the external object.
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
  Example    : my @synonyms = @{ $db_entry->get_all_synonyms };
  Description: Get a list of synonyms known for this object.
               Synonyms are lazy-loaded if required.
  Returntype : listref of strings. May be empty.
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_synonyms {

  my $self = shift;

  # lazy-load synonyms if required
  if (!$self->{synonyms} && $self->{adaptor}) {
    $self->{synonyms} = $self->{adaptor}->fetch_all_synonyms($self->dbID());
  }

  return $self->{synonyms};
}


=head2 flush_synonyms

  Args       : none
  Example    : none
  Description: Remove all synonyms from this object.
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

  Arg [1]    : (optional) String $arg - value to set
  Example    : none
  Description: Getter/setter for attribute 'status'.
  Returntype : String
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

