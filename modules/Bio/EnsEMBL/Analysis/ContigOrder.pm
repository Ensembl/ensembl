
#
# BioPerl module for ContigOrder
#
# Cared for by Tim Hubbard <th@sanger.ac.uk>
#
# Copyright Tim Hubbard
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

ContigOrder - DESCRIPTION of Object

Not OO.  Subroutine calls transcript2contigOrder and contigOrder.

contigOrder takes a list of contigpair relationships (which can be in
either direction and orientation), such as contig1.f->contig3.r,
contig2.f->contig1.r etc. and returns a set of contigcluster objects,
such that each contains a best order, such as:

    contig3.f:contig2.f:contig1.r

The ordering will be ambigious: initially a montycarlo algorithm will
be used to find the optimal ordering which minimises the skipping
number of contigpair relationships which skip over a contig.

The evidence for this list can come from multiple sources: initially
from transcript objects, but could also be backend sequence matches.

transcript2contigOrder is a routine to build contigpair relationships
from the set of contigs of each exon in a transcript object.

=head1 SYNOPSIS

Give standard usage here

    &contigOrder(\@contigpairs,
		 \%contigclusters,
		 \$mess);

    &transcript2contigOrder(\@transcripts,
			    \@contigpairs,
			    \$mess);

=head1 DESCRIPTION

Describe the object here

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package EnsEMBL::ContigOrder;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;


@ISA = qw(Bio::Root::Object);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;

# set stuff in self from @args
 return $make; # success - we hope!
}
