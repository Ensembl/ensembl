# EnsEMBL module for Bio::EnsEMBL::Utils::Sequence
#
#

=head1 NAME

Bio::EnsEMBL::Utils::Sequence - Utility functions for sequences

=head1 SYNOPSIS

  use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);

  my $seq = 'ACTTTAAAGGCTATCCCAATATG';

  print "my sequence = $seq\n";

  reverse_comp(\$seq);

  print "my reverse comp = $seq\n";


=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::Utils::Sequence;

use Exporter;

use vars qw(@ISA @EXPORT_OK);

@ISA = qw(Exporter);

@EXPORT_OK = qw(&reverse_comp);


=head2 reverse_comp

  Arg [1]    : reference to a string $seqref
  Example    : use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);

               $seq = 'ACCTGAA';
               reverse_comp(\$seq);
               print $seq;

  Description: Does an in place reverse compliment of a passed in string
               reference.  The string is passed by reference
               rather than by value for memory efficiency.
  Returntype : none
  Exceptions : none
  Caller     : SequenceAdaptor, SliceAdaptor

=cut

sub reverse_comp {
  my $seqref = shift;

  $$seqref = reverse( $$seqref );
  $$seqref =~
    tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;

  return;
}


1;
