#
# Object for storing details of an analysis job
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

Bio::EnsEMBL::Analysis::Job

=head1 SYNOPSIS

=head1 DESCRIPTION

Object to store details of an analysis run

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Analysis::Job;

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

  my ($id,$db,$db_version,$program,$program_version,$gff_source,$gff_feature) = 

#      $self->_rearrange([qw(
#			    )],@args);

  return $make; # success - we hope!
}


=head2 queue

  Title   : queue
  Usage   : $self->queue
  Function: Get/set method for the LSF queue name
  Returns : String
  Args    : String

=cut

sub queue {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_queue} = $arg;
    }

    return $self->{_queue};
}


=head2 input

  Title   : input
  Usage   : $self->input
  Function: Get/set method for the input objects
            This should be overwritten by derived
            classes to check/manipulate the input
            objects into the correct format.
  Returns : \@ Objects
  Args    : \@ Objects

=cut

sub input {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_input} = $arg;
    }

    return $self->{_input};
}


=head2 submit

  Title   : submit
  Usage   : $self->submit
  Function: Submits the job to the specified LSF queue
  Returns : 
  Args    : 

=cut

sub submit {
    my ($self) = @_;

}


=head2 freeze

  Title   : freeze
  Usage   : $self->freeze
  Function: Freezes the object into a string
  Returns : String
  Args    : None

=cut

sub freeze {
    my ($self) = @_;

}


=head2 submission_checks

  Title   : submission_checks
  Usage   : $self->submission_checks
  Function: After submission to the LSF queue when 
            the wrapper script is run - these are
            the checks to run (on binaries,databases etc)
            before the job is run.
  Returns : String
  Args    : None

=cut

sub submission_checks {
    my ($self) = @_;

    return $self->{_submission_checks};
}


=head2 add_check

  Title   : add_check
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
