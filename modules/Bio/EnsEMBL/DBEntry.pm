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

Bio::EnsEMBL::DBEntry -
Object representing an external reference (xref)

=head1 DESCRIPTION

This object holds information about external references (xrefs) to
Ensembl objects.

=head1 METHODS

=cut

package Bio::EnsEMBL::DBEntry;

use strict;
use warnings;
no warnings qw(uninitialized);

use Bio::EnsEMBL::Storable;
use Bio::Annotation::DBLink;

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(deprecate);

our @ISA = qw(Bio::EnsEMBL::Storable Bio::Annotation::DBLink);

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
               Due to 'PRIORITY',
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
       $priority,
       $db_display_name, $info_type, $info_text, $type,
       $secondary_db_name, $secondary_db_table, $link_annotation, $analysis) =
    rearrange ( ['ADAPTOR','DBID','PRIMARY_ID','VERSION',
                 'DBNAME','RELEASE','DISPLAY_ID','DESCRIPTION',
		 'PRIORITY',
		 'DB_DISPLAY_NAME', 'INFO_TYPE', 'INFO_TEXT', 'TYPE',
                 'SECONDARY_DB_NAME', 'SECONDARY_DB_TABLE', 'LINKAGE_ANNOTATION', 'ANALYSIS'], @args );

  $self->adaptor($adaptor);
  $self->{'dbID'}    = $dbID;

  if( defined $primary_id ) { $self->primary_id( $primary_id ) }
  if( defined $version ) { $self->version( $version ) } else
    { $self->version( 0 ); }
  if( defined $dbname ) { $self->dbname( $dbname ) }
  if( defined $release) { $self->release( $release ) }
  if( defined $display_id) { $self->display_id( $display_id ) }
  if( defined $description) { $self->description($description) }
  if( defined $priority) { $self->priority($priority) }
  if( defined $db_display_name) { $self->db_display_name($db_display_name) }
  if( defined $info_type) { $self->info_type($info_type) }
  if( defined $info_text) { $self->info_text($info_text) }
  if( defined $type) { $self->type($type) }
  if( defined $secondary_db_name) { $self->secondary_db_name($secondary_db_name) }
  if( defined $secondary_db_table) { $self->secondary_db_table($secondary_db_table) }

  $self->linkage_annotation($link_annotation) if defined $link_annotation;
  $self->analysis($analysis) if defined $analysis;


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

=head2 db_version

  Arg [1]    : (optional) String $arg - value to set
  Example    : none
  Description: Alias for release(). The release/version of the external DB
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub db_version {
  my ( $self, $arg ) = @_;
  return $self->release($arg);
}



=head2 description

  Arg [1]    : (optional) String $arg - value to set
  Example    : none
  Description: Getter/setter for attribute 'description'.
               The object's description, from the xref table
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub description {
  my ( $self, $arg ) = @_;

  if ( defined($arg) ) { $self->{'description'} = $arg }

  return $self->{description};
}

=head2 analysis

  Arg [1]    : Bio::EnsEMBL::Analysis $analysis
  Example    : none
  Description: get/set for attribute analysis
  Returntype : Bio::EnsEMBL::Analysis
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub analysis {
   my $self = shift;
  $self->{analysis} = shift if( @_ );
  return $self->{analysis};
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


=head2 priority

  Arg [1]    : int $priority
  Example    : none
  Priority   : Getter/setter for attribute 'priority'.
               The external database priority.
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
               Defines the method by which the linked object was
               connected to this xref. Available types are:
               'CHECKSUM', 'COORDINATE_OVERLAP', 'DEPENDENT',
               'DIRECT', 'INFERRED_PAIR', 'MISC', 'NONE', 'PROBE',
               'PROJECTION', 'SEQUENCE_MATCH' AND 'UNMAPPED'.
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
               Additional information recorded during the xref
               association. Intended to be used as metadata.
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
               The object xref 'linkage annotation'.
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
               The external database type.
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
  Description: Getter/setter for attribute 'secondary_db_name'.
               The external database 'secondary' database name.
  Returntype : String
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
  Description: Getter/setter for attribute 'secondary_db_table'.
               The external database 'secondary' database table.
  Returns    : String
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

  Args      : none
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
  if (!$self->{synonyms} && $self->adaptor()) {
    $self->{synonyms} = $self->adaptor()->fetch_all_synonyms($self->dbID());
  }

  return $self->{synonyms};
}


=head2 get_all_dependents

  Args[1]    : (optional) Bio::EnsEMBL::Gene, Transcript or Translation object
  Example    : my @dependents = @{ $db_entry->get_all_dependents };
  Description: Get a list of DBEntrys that are depenednet on the DBEntry.
               if an ensembl gene transcript or translation is given then only
               the ones on that object will be given
  Returntype : listref of DBEntrys. May be empty.
  Exceptions : none
  Caller     : general
  Status     : UnStable

=cut

sub get_all_dependents {
  my $self = shift;
  my $ensembl_object = shift;

  return  $self->adaptor()->get_all_dependents($self->dbID(), $ensembl_object);
}

=head2 get_all_masters

  Args[1]    : (optional) Bio::EnsEMBL::Gene, Transcript or Translation object
  Example    : my @masters = @{ $db_entry->get_all_masters };
  Description: Get a list of DBEntrys that are the masters of the DBEntry.
               if an ensembl gene transcript or translation is given then only
               the ones on that object will be given.
  Returntype : listref of DBEntrys. May be empty.
  Exceptions : none
  Caller     : general
  Status     : UnStable

=cut

sub get_all_masters {
  my $self = shift;
  my $ensembl_object = shift;

  return  $self->adaptor()->get_all_masters($self->dbID(), $ensembl_object);
}


=head2 flush_synonyms

  Args       : none
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
  Description: Getter/setter for attribute 'status'.
               The external database status.
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

1;

