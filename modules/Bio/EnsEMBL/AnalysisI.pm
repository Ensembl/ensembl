
#
# BioPerl module for Bio::EnsEMBL::AnalysisI
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::AnalysisI - Abstract interface for an analysis

=head1 SYNOPSIS

  # analysis objects are attached to EnsEMBL::SeqFeatureI objects

  $analysis = $sf->analysis();
  # has the following methods
  $analysis->db;
  $analysis->db_version;
  $analysis->program;
  $analysis->program_version;
  
=head1 DESCRIPTION

Abstract analysis object, indicating what methods one can call on
analysis objects.


=head1 AUTHOR - Ewan Birney

Email birney@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::AnalysisI;
use vars qw($AUTOLOAD @ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object

use Bio::EnsEMBL::Root;

@ISA = qw ( Bio::EnsEMBL::Root );


=head2 db

 Title   : db
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub db{
   my ($self) = @_;

   $self->throw("No database method provided on analysis implementing object");
}

=head2 db

 Title   : db_version
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub db_version{
   my ($self,@args) = @_;

   $self->throw("no database_version method on analysis implementing object");
}

=head2 program

 Title   : program
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub program{
   my ($self,@args) = @_;

   $self->throw("no program method on analysis object");
}

=head2 program_version

 Title   : program_version
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub program_version{
   my ($self,@args) = @_;

   $self->throw("No program version on analysis implementing object");
}

=head2 has_database

 Title   : has_database
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub has_database {
   my ($self,@args) = @_;

   $self->throw("No has_database analysis implementing object");

}


1;







