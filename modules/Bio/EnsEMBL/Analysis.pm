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

Bio::EnsEMBL::Analysis.pm - Stores details of an analysis run

=head1 SYNOPSIS

    my $obj = new Bio::EnsEMBL::Analysis(
      -id              => $id,
      -logic_name      => 'SWIRBlast',
      -db              => $db,
      -db_version      => $db_version,
      -db_file         => $db_file,
      -program         => $program,
      -program_version => $program_version,
      -program_file    => $program_file,
      -gff_source      => $gff_source,
      -gff_feature     => $gff_feature,
      -module          => $module,
      -module_version  => $module_version,
      -parameters      => $parameters,
      -created         => $created,
      -description     => 'some warm words about this analysis',
      -display_label   => 'UNIprot alignment',
      -displayable     => '1',
      -web_data        => 'web metadata info'
    );

=head1 DESCRIPTION

Object to store details of an analysis run.

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Storable;

use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

use base 'Bio::EnsEMBL::Storable';


=head2 new

  Arg [..]   :  Takes a set of named arguments
  Example    : $analysis = new Bio::EnsEMBL::Analysis::Analysis(
                                -id              => $id,
                                -logic_name      => 'SWIRBlast',
                                -db              => $db,
                                -db_version      => $db_version,
                                -db_file         => $db_file,
                                -program         => $program,
                                -program_version => $program_version,
                                -program_file    => $program_file,
                                -gff_source      => $gff_source,
                                -gff_feature     => $gff_feature,
                                -module          => $module,
                                -module_version  => $module_version,
                                -parameters      => $parameters,
                                -created         => $created );
  Description: Creates a new Analysis object
  Returntype : Bio::EnsEMBL::Analysis
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub new {
  my($class,@args) = @_;

  my $self = bless {},$class;

  my ($id, $adaptor, $db, $db_version, $db_file, $program, $program_version,
      $program_file, $gff_source, $gff_feature, $module, $module_version,
      $parameters, $created, $logic_name, $description, $display_label,
      $displayable, $web_data) =

  rearrange([qw(ID
    ADAPTOR
    DB
    DB_VERSION
    DB_FILE
    PROGRAM
    PROGRAM_VERSION
    PROGRAM_FILE
    GFF_SOURCE
    GFF_FEATURE
    MODULE
    MODULE_VERSION
    PARAMETERS
    CREATED
    LOGIC_NAME
    DESCRIPTION
    DISPLAY_LABEL
    DISPLAYABLE
    WEB_DATA
  )],@args);

  $displayable ||= 0;

  $self->dbID             ($id);
  $self->adaptor        ($adaptor);
  $self->db             ($db);
  $self->db_version     ($db_version);
  $self->db_file        ($db_file);
  $self->program        ($program);
  $self->program_version($program_version);
  $self->program_file   ($program_file);
  $self->module         ($module);
  $self->module_version ($module_version);
  $self->gff_source     ($gff_source);
  $self->gff_feature    ($gff_feature);
  $self->parameters     ($parameters);
  $self->created        ($created);
  $self->logic_name ( $logic_name );
  $self->description( $description );
  $self->display_label( $display_label );
  $self->displayable( $displayable );
  $self->web_data       ( $web_data );
  return $self; # success - we hope!
}


=head2 db

  Arg [1]    : string $db
  Description: get/set for the attribute db
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub db {
    my ($self,$arg) = @_;

    if (defined($arg)) {
      $self->{_db} = $arg;
    }

    return $self->{_db};
}


=head2 db_version

  Arg [1]    : string $db_version
  Description: get/set for attribute db_version
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub db_version {
    my ($self,$arg) = @_;

    if (defined($arg)) {
      $self->{_db_version} = $arg;
    }

    return $self->{_db_version};
}


=head2 db_file

  Arg [1]    : string $db_file
  Description: get/set for attribute db_file
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub db_file {
    my ($self,$arg) = @_;

    if (defined($arg)) {
      $self->{_db_file} = $arg;
    }

    return $self->{_db_file};
}



=head2 program

  Arg [1]    : string $program
  Description: get/set for attribute program
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub program {
    my ($self,$arg) = @_;

    if (defined($arg)) {
      $self->{_program} = $arg;
    }

    return $self->{_program};
}


=head2 program_version

  Arg [1]    : string $program_version
  Description: get/set for attribute program_version
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub program_version {
    my ($self,$arg) = @_;

    if (defined($arg)) {
      $self->{_program_version} = $arg;
    }

    return $self->{_program_version};
}


=head2 program_file

  Arg [1]    : string $program_file
  Description: get/set for attribute program_file
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub program_file {
    my ($self,$arg) = @_;

    if (defined($arg)) {
      $self->{_program_file} = $arg;
    }

    return $self->{_program_file};
}


=head2 module

  Arg [1]    : string $module
  Description: get/set for attribute module. Usually a RunnableDB perl 
               module that executes this analysis job. 
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub module {
    my ($self,$arg) = @_;

    if (defined($arg)) {
      $self->{_module} = $arg;
    }

    return $self->{_module};
}


=head2 module_version

  Arg [1]    : string $module_version
  Description: get/set for attribute module_version
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub module_version {
    my ($self,$arg) = @_;

    if (defined($arg)) {
      $self->{_module_version} = $arg;
    }

    return $self->{_module_version};
}


=head2 gff_source

  Arg [1]    : string $gff_source
  Description: get/set for attribute gff_source
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub gff_source {
    my ($self,$arg) = @_;

    if (defined($arg)) {
      $self->{_gff_source} = $arg;
    }

    return $self->{_gff_source};
}


=head2 gff_feature

  Arg [1]    : string $gff_feature
  Description: get/set for attribute gff_feature
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub gff_feature {
    my ($self,$arg) = @_;

    if (defined($arg)) {
      $self->{_gff_feature} = $arg;
    }

    return $self->{_gff_feature};
}


=head2 parameters

  Arg [1]    : string $parameters
  Description: get/set for attribute parameters. This should be evaluated
               by the module if given or the program that is specified.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub parameters {
    my ($self,$arg) = @_;

    if (defined($arg)) {
      $self->{_parameters} = $arg;
    }

    return $self->{_parameters};
}


=head2 created

  Arg [1]    : string $created
  Description: get/set for attribute created time.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub created {
    my ($self,$arg) = @_;

    if (defined($arg)) {
      $self->{_created} = $arg;
    }

    return $self->{_created};
}


=head2 logic_name

  Arg [1]    : string $logic_name
  Description: Get/set method for the logic_name, the name under 
               which this typical analysis is known.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub logic_name {
  my ($self, $arg ) = @_;
  ( defined $arg ) &&
    ($self->{_logic_name} = $arg);
  $self->{_logic_name};
}


=head2 has_database

  Args       : none
  Description: tests if the db attribute is set, returns 1 if so,
               0 if not.
  Returntype : int 0,1
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub has_database{
   my ($self,@args) = @_;

   if( defined $self->db ){ return 1; }
   return 0;
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
    my ($self,$arg) = @_;

    if (defined($arg)) {
      $self->{_description} = $arg;
    }

    return $self->{_description};
}


=head2 display_label

  Arg [1]    : string $display_label
  Description: get/set for attribute display_label
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub display_label {
    my ($self,$arg) = @_;

    if (defined($arg)) {
      $self->{_display_label} = $arg;
    }

    return $self->{_display_label};
}

=head2 displayable

  Arg [1]    : string $displayable
  Description: get/set for attribute displayable
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub displayable {
    my ($self,$arg) = @_;

    if (defined($arg)) {
      $self->{_displayable} = $arg;
    }

    return $self->{_displayable};
}


=head2 web_data

  Arg [1]    : string $web_data
  Description: get/set for attribute web_data
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
# See also _objFromHashref in AnalysisAdaptor for context.
sub web_data {
    my ($self,$arg) = @_;

    if (defined($arg)) {
      $self->{_web_data} = $arg;
    }

    return $self->{_web_data};
}

=head2 compare

  Arg  1     : Bio::EnsEMBL::Analysis $ana
               The analysis to compare to
  Description: returns 1 if this analysis is special case of given analysis
               returns 0 if they are equal
	           returns -1 if they are completely different
  Returntype : int -1,0,1
  Exceptions : none
  Caller     : unknown
  Status     : Stable

=cut

sub compare {
  my ($self, $ana ) = @_;

  throw("Object is not a Bio::EnsEMBL::Analysis") 
    unless $ana->isa("Bio::EnsEMBL::Analysis");

  my $detail = 0;

  foreach my $methodName ( 'program', 'program_version', 'program_file',
    'db','db_version','db_file','gff_source','gff_feature', 'module',
    'module_version', 'parameters','logic_name' ) {
    if( defined $self->$methodName() && ! $ana->can($methodName )) {
      $detail = 1;
    }
    if( defined $self->$methodName() && ! defined $ana->$methodName() ) {
      $detail = 1;
    }
    # if given anal is different from this, defined or not, then its different
    if( defined($ana->$methodName()) && defined($self->$methodName()) &&
        ( $self->$methodName() ne $ana->$methodName() )) {
      return -1;
    }
  }
  if( $detail == 1 ) { return 1 };
  return 0;
}


1;
















