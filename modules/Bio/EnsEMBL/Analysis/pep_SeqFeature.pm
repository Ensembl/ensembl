
#
# BioPerl module for Bio::SeqFeature::Homol
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Analysis::pep_SeqFeature

=head1 SYNOPSIS

=head1 DESCRIPTION

Extends the bioperl Bio::SeqFeature::Homol to cope with 
fractional start/end points in peptides.

Creation:

    my $pepf = new Bio::EnsEMBL::Analysis::pep_SeqFeature(-start  => $start,
						-end    => $end,
						-strand => $strand);

       $pepf->start_frac(2);
       $pepf->end_frac  (3);

Manipulation:

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Analysis::pep_SeqFeature;

use vars qw($AUTOLOAD @ISA);
use strict;

# Object preamble - inherits from Bio::SeqFeature::Homol

use Bio::SeqFeature::Homol;


@ISA = qw(Bio::SeqFeature::Homol);

# new() is inherited from Bio::SeqFeature::Generic

sub _initialize {
  my($self,@args) = @_;
  my $make = $self->SUPER::_initialize(@args);

  # Defaults
  $self->start_frac(1);
  $self->end_frac  (3);

  return $make;
}

=head1 Methods unique to exon

=head2 start_frac

 Title   : start_frac
 Usage   : $pepf->start_frac(3)
 Function: Get/set method for the start fraction of a residue
 Example : 
 Returns : 1,2,3
 Args    : 1,2,3


=cut

sub start_frac {
    my ($self,$arg) = @_;

    if ($arg ne "") {
	$self->throw("ERROR: only 1,2,3 allowed for end_frac : $arg") unless ($arg == 1 || $arg == 2 || $arg == 3);

	$self->{_start_frac} = $arg;
    }
    
    return $self->{_start_frac};
}

=head2 end_frac

 Title   : end_frac
 Usage   : $pepf->end_frac(3)
 Function: Get/set method for the end fraction of a residue
 Example : 
 Returns : 1,2,3
 Args    : 1,2,3


=cut

sub end_frac {
    my ($self,$arg) = @_;

    if ($arg ne "") {
    $self->throw("ERROR: only 1,2,3 allowed for end_frac : $arg") unless ($arg == 1 || $arg == 2 || $arg == 3);


	$self->{_end_frac} = $arg;
    }
    
    return $self->{_end_frac};
}




