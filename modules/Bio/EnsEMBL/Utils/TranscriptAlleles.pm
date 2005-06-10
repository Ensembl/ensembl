#
# Ensembl module for TranscriptAlleles
#
# Copyright (c) 2005 Ensembl
#
# You may distribute this module under the same terms as perl itself
#

=head1 NAME

TranscriptAlleles - A utility class used to obtain information about the
relationships between a transcript and Alleles

=head1 SYNOPSIS

  use Bio::EnsEMBL::Utils::TranscriptAlleles;


  # get the peptide variations caused by a set of Alleles

  %variations = %{Bio::EnsEMBL::Utils::TranscriptAlleles::get_all_peptide_variations($transcript, $alleles)};


=head1 DESCRIPTION

This is a utility class which can be used to find consequence type of an AlleleFeature in a
transcript, and to determine the amino acid changes caused by the AlleleFeature in the Transcript

=head1 CONTACT

Email questions to the ensembl developer mailing list <ensembl-dev@ebi.ac.uk>

=head1 METHODS

=cut
package Bio::EnsEMBL::Utils::TranscriptAlleles;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Variation::ConsequenceType;
use vars qw(@ISA @EXPORT_OK);

@ISA = qw(Exporter);

@EXPORT_OK = qw(&get_all_ConsequenceType &type_variation);

=head2 get_all_ConsequenceType

  Arg [1]    : $transcript the transcript to obtain the peptide variations for
  Arg [2]    : $alleles listref of AlleleFeatures
  Example    : $consequence_types = get_all_ConsequenceType($transcript, \@alleles);
               foreach my $ct (@{$consequence_types}){
                  print "Allele : ", $ct->allele_string, " has a consequence type of :",$ct->type;
                  print " and is affecting the transcript with ",@{$ct->aa_alleles}, "in position ", 
		              $ct->aa_start,"-", $ct->aa_end if (defined $ct->aa_alleles);
		  print "\n";
	      }
  Description: Takes a list of AlleleFeatures and a Transcritpt, and return a list
                     of ConsequenceType of the alleles in the given Transcript
  Returntype : listref of Bio::EnsEMBL::Variation::ConsequenceType
  Exceptions : none
  Caller     : general 

=cut

sub get_all_ConsequenceType {
  my $transcript = shift;
  my $alleles = shift;

  if(!ref($transcript) || !$transcript->isa('Bio::EnsEMBL::Transcript')) {
    throw('Bio::EnsEMBL::Transcript argument is required.');
  }

  if(!ref($alleles) eq 'ARRAY') {
    throw('Reference to a list of Bio::EnsEMBL::Variation::AlleleFeature objects is required');
  }


  my @alleles_ordered = sort { $a->start <=> $b->start} @$alleles; #sort the alleles by the genomic position
  my @same_codon; #contains up to 3 allele features, that are in the same codon
  my @out; #array containing the consequence types of the alleles in the transcript
  foreach my $allele (@alleles_ordered) {
    #get consequence type of the AlleleFeature
    my $consequence_type = Bio::EnsEMBL::Variation::ConsequenceType->new($transcript->dbID(),'',$allele->start,$allele->end,$allele->strand,[$allele->allele_string]);
    #calculate the consequence type of the Allele
    my $ref_consequences = type_variation($transcript,$consequence_type);
    if ($allele->start != $allele->end){
	#do not calculate for indels effects of 2 or more in same codon
	push @out, @{$ref_consequences};
	next;
    }

    my $new_consequence = shift @{$ref_consequences};
    if (!defined $new_consequence->aa_start){
	push @out, $new_consequence;
	next;
    }
    #for alleles with aa effect, find out if they are in the same codon
    if (!defined $same_codon[0] || $same_codon[0]->aa_start == $new_consequence->aa_start){
	push @same_codon,$new_consequence;
    }
    else{
	#if there is more than one element in the same_codon array, calculate the effect of the codon
	if (@same_codon > 1){
	    calculate_same_codon(\@same_codon);
	}
	push @out, @same_codon;
	@same_codon = ();
	push @same_codon, $new_consequence; #push the element not in the same codon
    }
  }
  #add last consequence_type
  if (@same_codon == 1){
      push @out, @same_codon;
  }
  elsif (@same_codon > 1){
    calculate_same_codon(\@same_codon);
    push @out, @same_codon;
  }

  return \@out;
}

# recalculates the effect of 2 or 3 SNPs in the same codon
sub calculate_same_codon{
    my $same_codon = shift;
    my $new_codon;
    my $old_aa;
    my $codon_table = Bio::Tools::CodonTable->new;
    if (@{$same_codon} == 3){
	#if there are 3 alleles in the same codon
	map {$new_codon .= @{$_->alleles};$old_aa = $_->aa_alleles()->[0]} @{$same_codon};
    }
    else{
	#if there are 2 alleles affecting the same codon
	my $first_pos = ($same_codon->[0]->cdna_start -1) % 3; #position of the first allele in the codon
	my $second_pos = ($same_codon->[1]->cdna_start -1)% 3; #position of the second allele in the codon
	if ($first_pos == 0){
	    #codon starts with first allele
	    $new_codon = $same_codon->[0]->alleles->[0]; #first base in the codon
	    if ($second_pos == 1){
		$new_codon .= $same_codon->[1]->alleles->[0]; #second base in the codon
		$new_codon .= substr($same_codon->[1]->codon,2,1); #third base in the codon
	    }
	    else{
		$new_codon .= substr($same_codon->[1]->codon,1,1); #second base in the codon
		$new_codon .= $same_codon->[1]->alleles->[0]; #third base in the codon
	    }
	}
	else{
	    #alleles are in position 1 and 2 in the codon
	    $new_codon = substr($same_codon->[1]->codon,0,1); #first base in the codon
	    $new_codon .= $same_codon->[0]->alleles->[0]; #second base in the codon
	    $new_codon .= $same_codon->[1]->alleles->[0]; #third base in the codon
	}
	$old_aa = $same_codon->[0]->aa_alleles->[0];	
    }
    #calculate the new_aa
    my $new_aa = $codon_table->translate($new_codon);
    #and update the aa_alleles field
    map {$_->aa_alleles([$old_aa,$new_aa])} @{$same_codon};

}
#
# Classifies a variation which is in the vicinity of a transcript
#
sub type_variation {
  my $tr  = shift;
  my $var = shift;

  my $UPSTREAM = 5000;
  my $DOWNSTREAM = 5000;

  if (!$var->isa('Bio::EnsEMBL::Variation::ConsequenceType')){
      throw("Not possible to calculate the consequence type for ",ref($var)," : Bio::EnsEMBL::Variation::ConsequenceType object expected");      
  }

  # ARNE: need the following:
  # $var and $tr need to have the same slice attached
  # if not, its difficult to compare the variations and the 
  # transcript. If one or both are on a strain slice the 
  # coordinate transformations DONT WORK CORRECTLY

  # if( $tr->slice() != $var->slice() ) {

#     # problem: Are the coordinates comparable?
#     # throw here, or warn and return, otherwise its too difficult
#     if( $tr->slice()->is_not_strainSlice() &&
# 	$var->slice()->is_not_strainSlice() ) {
#       # coordintes are transferable
#     } else {
#       # coordinates are not transferable
#     }
#   } else {
#     # yes they are comparable
#   }

  #necessary to compare the variation with the transcript to determine if it is in the
  #vicinity or not
  if (($var->start < $tr->start - $UPSTREAM) || ($var->end  > $tr->end + $DOWNSTREAM)){
      #since the variation is more than UPSTREAM and DOWNSTREAM of the transcript, there is no effect in the transcript
      return [];
  }
  
  my $tm = $tr->get_TranscriptMapper();
  my @coords = $tm->genomic2cdna($var->start,
                                 $var->end,
                                 $var->strand);

  # Handle simple cases where the variation is not split into parts.
  # Call method recursively with component parts in complicated case.
  # E.g. a single multi-base variation may be both intronic and coding


  if(@coords > 1) {
    my @out;

    foreach my $c (@coords) {
      my %new_var = %{$var};
      $new_var{'end'} = $var->start + $c->length() - 1;
      $var->start( $new_var{'end'} + 1);      
      push @out, @{type_variation($tr, bless \%new_var, ref($var))};
    }

    return \@out;
  }

  my $c = $coords[0];
  if($c->isa('Bio::EnsEMBL::Mapper::Gap')) {

    # check if the variation is completely outside the transcript:

    if($var->end < $tr->start()) {
      $var->type( ($tr->strand() == 1) ? 'UPSTREAM' : 'DOWNSTREAM' );
      return [$var];
    }
    if($var->start > $tr->end()) {
      $var->type( ($tr->strand() == 1) ? 'DOWNSTREAM' : 'UPSTREAM' );
      return [$var];
    }

    # variation must be intronic since mapped to cdna gap, but is within
    # transcript
    $var->type('INTRONIC');
    return [$var];
  }

  $var->cdna_start( $c->start() );
  $var->cdna_end( $c->end() );

  @coords = $tm->genomic2cds($var->start, $var->end,$var->strand);

  if(@coords > 1) {
    my @out;

    foreach my $c (@coords) {
      my %new_var = %{$var};
      $new_var{'end'} = $var->start + $c->length() - 1;
      $var->start( $new_var{'end'} + 1);
      push @out, @{type_variation($tr, bless \%new_var, ref($var))};
    }
    return \@out;
  }

  $c = $coords[0];

  if($c->isa('Bio::EnsEMBL::Mapper::Gap')) {
    # mapped successfully to CDNA but not to CDS, must be UTR

    if($var->end < $tr->coding_region_start()) {
      $var->type( ($tr->strand() == 1) ? '5PRIME_UTR' : '3PRIME_UTR' );
    }
    elsif($var->start > $tr->coding_region_end()) {
      $var->type( ($tr->strand() == 1) ? '3PRIME_UTR' : '5PRIME_UTR');
    }
    else {
      throw('Unexpected: CDNA variation which is not in CDS is not in UTR');
    }
    return [$var];
  }

  $var->cds_start( $c->start() );
  $var->cds_end( $c->end() );

  @coords = $tm->genomic2pep($var->start, $var->end, $var->strand);


  if(@coords != 1 || $coords[0]->isa('Bio::EnsEMBL::Mapper::Gap')) {
    throw("Unexpected: Could map to CDS but not to peptide coordinates.");
  }

  my @coords_amp = $tm->genomic2pep($var->start -2, $var->end +2, $var->strand);

  if (@coords_amp != @coords){
      $var->type( 'SPLICE_SITE' );
  }

  $c = $coords[0];

  $var->aa_start( $c->start());
  $var->aa_end( $c->end());

  apply_aa_change($tr, $var);

  return [$var];
}

#
# Determines the effect of a coding variation on the peptide sequence
#

sub apply_aa_change {
  my $tr = shift;
  my $var = shift;

  #my $peptide = $tr->translate->seq();
  #to consider stop codon as well
  my $mrna = $tr->translateable_seq();
  my $mrna_seqobj = Bio::Seq->new( -seq        => $mrna,
				   -moltype    => "dna",
				   -alphabet   => "dna");

  my ($attrib) = @{$tr->slice()->get_all_Attributes('codon_table')}; #for mithocondrial dna it is necessary to change the table

  my $codon_table;
  $codon_table = $attrib->value() if($attrib);
  $codon_table ||= 1; # default vertebrate codon table 

  my $peptide = $mrna_seqobj->translate(undef,undef,undef,$codon_table)->seq;
  my $len = $var->aa_end - $var->aa_start + 1;
  my $old_aa = substr($peptide, $var->aa_start -1 , $len);
  my $codon_cds_start = $var->aa_start * 3 - 2;
  my $codon_cds_end   = $var->aa_end   * 3;
  my $codon_len = $codon_cds_end - $codon_cds_start + 1;

  my @alleles = @{$var->alleles};
  
#  shift(@alleles); # ignore reference allele

  my $var_len = $var->cds_end - $var->cds_start + 1;

  my @aa_alleles = ($old_aa);

  foreach my $a (@alleles) {
    $a =~ s/\-//;
    my $cds = $tr->translateable_seq();
    if($var_len != length($a)) {
      if(abs(length($a) - $var_len) % 3) {
        # frameshifting variation, do not set peptide_allele string
        # since too complicated and could be very long
        $var->type('FRAMESHIFT_CODING');
        return;
      }

      if($codon_len == 0) { # insertion
        $aa_alleles[0] = '-';
        $old_aa    = '-';
      }
    }

    my $new_aa;

    if(length($a)) {
      substr($cds, $var->cds_start-1, $var_len) = $a;
      my $codon_str = substr($cds, $codon_cds_start-1, $codon_len + abs(length($a)-$var_len));
      
      $var->codon($codon_str); #add the codon to the ConsequenceType object
      my $codon_seq = Bio::Seq->new(-seq      => $codon_str,
                                    -moltype  => 'dna',
                                    -alphabet => 'dna');


      $new_aa = $codon_seq->translate(undef,undef,undef,$codon_table)->seq();
    } else {
      $new_aa = '-'; # deletion
    }

    if(uc($new_aa) ne uc($old_aa)) {
      push @aa_alleles, $new_aa;
      if ($new_aa =~ /\*/) {
	$var->type('STOP_GAINED');
      }
      elsif ($old_aa =~ /\*/) {
	$var->type('STOP_LOST');
      }
    }
  }

  if(@aa_alleles > 1) {
    if (! defined $var->type) {
      $var->type('NON_SYNONYMOUS_CODING');
    }
  } else {
    $var->type('SYNONYMOUS_CODING');
  }

  $var->aa_alleles(\@aa_alleles);
}


1;
