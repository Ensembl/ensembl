
#
# BioPerl module for Exon2Gene
#
# Cared for by Tim Hubbard <th@sanger.ac.uk>
#
# Copyright Tim Hubbard
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Exon2Gene - DESCRIPTION of Object

Not OO.  Single subroutine call exon2gene.

Given a list of exon pairings (x followed by y) builds a list of all
transcripts implied by the list (uses all starts, and follows all
branches).  Issues warnings on finding circular references and
terminates transcript construction at that point.  Clusters
transcripts into gene objects.

Later: given an existing list of transcripts and genes, extends
existing transcripts, where extra exons have been added.  Allocates
new transcript and gene identifiers where no existing object found.

=head1 SYNOPSIS

Give standard usage here

    &exon2gene(\@exons,\%pairs,
	       \%transcripts,\%genes,
	       \$mess);

Where

@exons is an array of exon objects

%pairs is a hash of pairwise relationships.

Unless you want objects for these too, for now I propose using the
format in my old script, which is:
$pairs{HExxxxxxxxx}="HExxxxxxxxy:abce;HExxxxxxxxz:abcd", which
indicates the relationships HExxxxxxxxx->HExxxxxxxxy and
HExxxxxxxxx->HExxxxxxxxz

%transcripts is a hash of transcript objects, indexed by unique id
(HTxxxxxxx).  This is what I previously called a 'gene'.

%genes is a hash of gene objects, indexed by unique id (HGxxxxxxx).
This is what I previously called a 'gene cluster'.

warning messages from the routine get written to $mess

Notes:

@exons is not really required - currently I check \%pairs against
@exons to see if there are inconsistencies, however if routine is to
be used in an incremental mode, this become irrelevant.

Incremental mode:

Initially routine would just build complete set of %transcripts,
%genes each time it is called.  Later it should modify contents of
existing %transcripts, %genes.  In this case either %pairs can be just
new pair relations or complete set.

Required:

Definition of transcript and gene objects

=head1 DESCRIPTION

Describe the object here

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Analysis::Exon2Gene;
use vars qw(@ISA);
use strict;


# since this is not currently OO following part is not used at present

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


sub exon2gene{
    my($raexons,$rhpairs,$rhtranscripts,$rhgenes,$logobj)=@_;
}
