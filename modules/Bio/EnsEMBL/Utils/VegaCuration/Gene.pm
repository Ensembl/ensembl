=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 METHODS

=cut

package Bio::EnsEMBL::Utils::VegaCuration::Gene;

use strict;
use warnings;
use vars qw(@ISA);

use Bio::EnsEMBL::Utils::ConversionSupport;

@ISA = qw(Bio::EnsEMBL::Utils::ConversionSupport);


=head2 find_gaps

   Args       : arrayref of B::E::Transcripts
   Example    : my $gaps = find_gaps($all_transcripts)
   Description: identifies regions of a gene that are not covered by any transcript
   Returntype : int
   Exceptions : none
   Caller     : internal

=cut

sub find_gaps {
  my $self = shift;
  my ($all_transcripts) = @_;
  my $gaps = 0;
  my @sorted_transcripts = sort {$a->start <=> $b->start || $b->end <=> $a->end} @{$all_transcripts};
  if ( my $first_transcript = shift @sorted_transcripts ) {
    my $pos = $first_transcript->end;
    foreach my $transcript (@sorted_transcripts) {
      next if ($transcript->end < $pos );
      if ($transcript->start < $pos && $transcript->end > $pos ) {
	$pos = $transcript->end;			
	next;
      }
      elsif ($transcript->end > $pos) {
	$gaps++;
	$pos = $transcript->end;
      }
    }
  }
  return $gaps;		
}
