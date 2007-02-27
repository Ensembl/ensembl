package Bio::EnsEMBL::Utils::VegaCuration::Gene;

=head1 NAME

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 LICENCE

This code is distributed under an Apache style licence:
Please see http://www.ensembl.org/code_licence.html for details

=head1 AUTHOR

Steve Trevanion <st3@sanger.ac.uk>

=head1 CONTACT

Post questions to the EnsEMBL development list ensembl-dev@ebi.ac.uk

=cut

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
	my $first_transcript = shift @sorted_transcripts;
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
	return $gaps;		
}
