#
# Object for submitting jobs to and querying the LSF queue
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

Bio::EnsEMBL::MappedExon

=head1 SYNOPSIS

=head1 DESCRIPTION

Extension of a normal exon that knows whether its DNA 
sequence has changed from the previous version.

Usage is the same for Bio::EnsEMBL::Exon except we have one
new method 

  $exon->has_identical_sequence(0)  # Sequence not identical to last version

  or

  $exon->has_identical_sequence(1)  # Sequence is the same as the last version

  and

  my $changed = $exon->has_identical_sequence


=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::MappedExon;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object;

use Bio::EnsEMBL::Exon;

@ISA = qw(Bio::EnsEMBL::Exon Bio::Root::Object);

sub _initialize {
    my ($self,@args) = @_;

    my $make = $self->SUPER::_initialize(@args);

    return $make; # success - we hope!
}

=head2 has_identical_sequence

  Title   : has_identical_sequence
  Usage   : my $changed = $exon->has_identical_sequence
  Function: Get/set flag for whether the exon\'s sequence
            has changed.
  Returns : 0,1
  Args    : 0,1

=cut


sub has_identical_sequence {
    my ($self,$arg) = @_;

    if (defined($arg)) {

	if ($arg != 0 && $arg != 1) {
	    $self->throw("Argument to has_identical_sequence should be 0,1 [$arg]");
	}
	$self->{_identical_sequence} = $arg;
    }

    return $self->{_identical_sequence};

}

1;



