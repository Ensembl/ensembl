#
# Object for storing sequence analysis details
#
# Cared for by Michele Clamp  <michele@sanger.ac.uk>
#
# Copyright Michele Clamp
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Analysis::Analysis.pm - Stores details of an analysis run

=head1 SYNOPSIS

    my $obj    = new Bio::EnsEMBL::Analysis::Analysis(-db              => $db,
						      -db_version      => $db_version,
						      -program         => $program,
						      -program_version => $program_version,
						      -gff_source      => $gff_source,
						      -gff_feature     => $gff_feature,
						      )

=head1 DESCRIPTION

Object to store details of an analysis run

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Analysis::Analysis;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object;

use Bio::Root::Object;

# Inherits from the base bioperl object
@ISA = qw(Bio::Root::Object);

# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;
  
  my $make = $self->SUPER::_initialize;


  return $self; # success - we hope!
}


=head2 id

  Title   : id
  Usage   : $self->id
  Function: Get/set method for the id
  Returns : int
  Args    : int

=cut

sub id {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_id} = $arg;
    }

    return $self->{_id};
}


=head2 db

  Title   : db
  Usage   : $self->db
  Function: Get/set method for database
  Returns : String
  Args    : String

=cut

sub db {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_db} = $arg;
    }

    return $self->{_db};
}


=head2 db_version

  Title   : db_version
  Usage   : $self->db_version
  Function: Get/set method for the database version number
  Returns : int
  Args    : int

=cut

sub db_version {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_db_version} = $arg;
    }

    return $self->{_db_version};
}


=head2 program

  Title   : program
  Usage   : $self->program
  Function: Get/set method for the program name
  Returns : String
  Args    : String

=cut

sub program {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_program} = $arg;
    }

    return $self->{_program};
}


=head2 program_version

  Title   : program_version
  Usage   : $self->program_version
  Function: Get/set method for the program version number
  Returns : int
  Args    : int

=cut

sub program_version {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_program_version} = $arg;
    }

    return $self->{_program_version};
}


=head2 gff_source

  Title   : gff_source
  Usage   : $self->gff_source
  Function: Get/set method for the gff_source tag
  Returns : String
  Args    : String

=cut

sub gff_source {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_gff_source} = $arg;
    }

    return $self->{_gff_source};
}

=head2 gff_feature

  Title   : gff_feature
  Usage   : $self->gff_feature
  Function: Get/set method for the gff_feature tag
  Returns : String
  Args    : String

=cut

sub gff_feature {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_gff_feature} = $arg;
    }

    return $self->{_gff_feature};
}

1;
