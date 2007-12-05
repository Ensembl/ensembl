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

	#hashref to return results
	my $results;

	#extra CDS end not found remarks
	if (grep {$_->value eq 'CDS end not found'} @remarks) {
		if (   ($coding_end != $trans_end) 
            && ( grep {$_ eq $stop_codon} qw(TGA TAA TAG) ) ) {
			$results->{'END_EXTRA'} = 1;
		}
	}

	#missing CDS end not found remark
	if ( $coding_end == $trans_end ) {
		if (! grep {$_->value eq 'CDS end not found'} @remarks) {
			if (grep {$_ eq $stop_codon} qw(TGA TAA TAG)) {
				$results->{'END_MISSING_2'} = 1;
			}
			else {
				$results->{'END_MISSING_1'} = $stop_codon;
			}
		}
	}


	#extra CDS start not found remark
	if (grep {$_->value eq 'CDS start not found'} @remarks) {
		if (   ($coding_start != 1) 
			&& ($start_codon eq 'ATG') ) {
			$results->{'START_EXTRA'} = 1;
		}
	}

	#missing CDS start not found remark
	if ( $coding_start == 1) {
		if ( ! grep {$_->value eq 'CDS start not found'} @remarks) {
			if ($start_codon eq 'ATG') {
				$results->{'START_MISSING_2'} = 1;
			} else {
				$results->{'START_MISSING_1'} = $start_codon;
			}
		}
	}

	return $results;
}

=head2 check_CDS_end_remarks_loutre

   Args       : B::E::Transcript
   Example    : my $results = $support->check_CDS_end_remarks($transcript)
   Description: identifies incorrect 'CDS end...' transcript attribs
   Returntype : hashref

=cut

sub check_CDS_start_end_remarks_loutre {
	my $self = shift;
	my $trans = shift;

	# info for checking
	my @stops = qw(TGA TAA TAG);
	my %attributes;
	foreach my $attribute (@{$trans->get_all_Attributes()}) {
		$attributes{$attribute->code} = $attribute;
	}
	my $coding_end   = $trans->cdna_coding_end;
	my $coding_start = $trans->cdna_coding_start;
	my $trans_end    = $trans->length;
	my $trans_seq    = $trans->seq->seq;
	my $stop_codon   = substr($trans_seq, $coding_end-3, 3);
	my $start_codon  = substr($trans_seq, $coding_start-1, 3);

	#hashref to return results
	my $results;

	#extra CDS end not found remarks
	if ( ($attributes{'cds_end_NF'}->value == 1)
			 && ($coding_end != $trans_end) 
				 && ( grep {$_ eq $stop_codon} @stops) ) {
		$results->{'END_EXTRA'} = 1;
	}
	#missing CDS end not found remark
	if ( $coding_end == $trans_end ) {
		if ($attributes{'cds_end_NF'}->value == 0 ) {
			if (grep {$_ eq $stop_codon} @stops) {
				$results->{'END_MISSING_2'} = 1;
			}
			else {
				$results->{'END_MISSING_1'} = $stop_codon;
			}
		}
	}
	#extra CDS start not found remark
	if ( ($attributes{'cds_start_NF'}->value == 1 )
			 && ($coding_start != 1)
				 && ($start_codon eq 'ATG') ) {
		$results->{'START_EXTRA'} = 1;
	}
	#missing CDS start not found remark
	if ( $coding_start == 1) {
		if ( $attributes{'cds_start_NF'}->value == 0 ) {
			if ($start_codon eq 'ATG') {
				$results->{'START_MISSING_1'} = 1;
			} else {
				$results->{'START_MISSING_2'} = $start_codon;
			}
		}
	}
	return $results;
}
