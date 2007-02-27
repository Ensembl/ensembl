package Bio::EnsEMBL::Utils::VegaCuration::Translation;

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

use Bio::EnsEMBL::Utils::VegaCuration::Transcript;

@ISA = qw(Bio::EnsEMBL::Utils::VegaCuration::Transcript);

=head2 check_CDS_end_remarks

   Args       : B::E::Transcript
   Example    : my $results = $support->check_CDS_end_remarks($transcript)
   Description: identifies incorrect 'CDS end...' transcript remarks
   Returntype : hashref

=cut

sub check_CDS_start_end_remarks {
	my $self = shift;
	my $trans = shift;

	# info for checking
	my @remarks = @{$trans->get_all_Attributes('remark')};
	my $coding_end   = $trans->cdna_coding_end;
	my $coding_start = $trans->cdna_coding_start;
	my $trans_end    = $trans->length;
	my $trans_seq    = $trans->seq->seq;
	my $stop_codon   = substr($trans_seq, $coding_end-3, 3);
	my $start_codon  = substr($trans_seq, $coding_start-1, 3);

	#hasref to return results
	my $results;

	# check consistency of 'CDS end not found' and the translation end
	if (grep {$_ eq $stop_codon} qw(TGA TAA TAG)) {
		if (    ($coding_end != $trans_end)
			 && (grep {$_->value eq 'CDS end not found'} @remarks) ) {
			$results->{'END_EXTRA'} = 1;
		}
	}
	else {
		unless (grep {$_->value eq 'CDS end not found'} @remarks) {
			$results->{'END_MISSING'} = $stop_codon;
		}
	}

	# check consistency of 'mRNA start not found' tags and the translation start
	if ($start_codon eq 'ATG') {
		if (    ($coding_start != 1)
			 && (grep {$_->value eq 'CDS start not found'} @remarks) ) {
			$results->{'START_EXTRA'} = 1;
		}
	}
	else {
		unless (grep {$_->value eq 'CDS start not found'} @remarks) {
			$results->{'START_MISSING'} = $start_codon;
		}
	}

	return $results;
}
