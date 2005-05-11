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
use vars qw(@ISA @EXPORT_OK);

@ISA = qw(Exporter);

@EXPORT_OK = qw(&type_variation);

use Data::Dumper;

=head2 get_all_peptide_variations

  Arg [1]    : $transcript the transcript to obtain the peptide variations for
  Arg [2]    : $alleles listref of AlleleFeatures
  Example    : $pep_hash = get_all_peptide_variations($transcript, \@alleles);
  Description: Takes a list of coding alleles on this transcript in
               which are in cdna coordinates and returns a hash with peptide
               coordinate keys and listrefs of alternative amino acids as
               values.  

               Note that the peptide encoded by the reference sequence is
               also present in the results and that duplicate peptides
               (e.g. resulting from synonomous mutations) are discarded.
               It is possible to have greated than two peptides variations
               at a given location given adjacent or overlapping snps.
               Insertion/deletion variations are ignored by this method.
               Example of a data structure that could be returned:
               {  1  => ['I', 'M'],
                 10  => ['I', 'T'],
                 37  => ['N', 'D'],
                 56  => ['G', 'E'],
                 118 => ['R', 'K'],
                 159 => ['D', 'E'],
                 167 => ['Q', 'R'],
                 173 => ['H', 'Q'] }
  Returntype : hashref
  Exceptions : none
  Caller     : general

=cut

sub get_all_peptide_variations {
  my $transcript = shift;
  my $alleles = shift;

  if(!ref($transcript) || !$transcript->isa('Bio::EnsEMBL::Transcript')) {
    throw('Bio::EnsEMBL::Transcript argument is required.');
  }

  if(!ref($alleles) eq 'ARRAY') {
    throw('Reference to a list of Bio::EnsEMBL::Variation::AlleleFeature objects is required');
  }

  my $codon_table = Bio::Tools::CodonTable->new;
  my $codon_length = 3;
  my $cdna = $transcript->spliced_seq;
  print "cdna seq $cdna\n";
  my $variant_alleles;
  my $translation_start = $transcript->cdna_coding_start;
  print "translation start ",$translation_start,"\n";
  foreach my $allele (@$alleles) {
    #skip variations not on a single base
    next if ($allele->start != $allele->end);
    print "allele in coords ",$allele->start,"-",$allele->end,"\n";
    #convert the AlleleFeature in cdna coordinates
    my $new_allele = convert_to_cdna($transcript,$allele);
    print "allele in cdna oords ",$new_allele->start,"-",$new_allele->end,"\n";
    next if(!defined $new_allele); #the AlleleFeature it is not in a coding region
    my $start = $new_allele->start;
    my $strand = $new_allele->strand;

    #calculate offset of the nucleotide from codon start (0|1|2)
    my $codon_pos = ($start - $translation_start) % $codon_length;
    #calculate the peptide coordinate of the allele
    my $peptide = ($start - $translation_start +
		   ($codon_length - $codon_pos)) / $codon_length;

    #retrieve the codon
    my $codon = substr($cdna, $start - $codon_pos-1, $codon_length);
    print "start $start and codon $codon\n";
    #store each alternative allele by its location in the peptide
    my $allele_string = lc($new_allele->allele_string);
    print "allele $allele_string\n";
    next if $allele_string eq '-';       #skip deletions
    next if CORE::length($allele_string) != 1; #skip insertions

    
    if($strand == -1) {
	#complement the allele if the snp is on the reverse strand
	$allele_string =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
    }

    #create a data structure of variant alleles sorted by both their
    #peptide position and their position within the peptides codon
    $variant_alleles ||= {};
    if(exists $variant_alleles->{$peptide}) {
        my $alleles_arr = $variant_alleles->{$peptide}->[1];
        push @{$alleles_arr->[$codon_pos]}, $allele_string;
    } else {
        #create a list of 3 lists (one list for each codon position)
        my $alleles_arr = [[],[],[]];
        push @{$alleles_arr->[$codon_pos]}, $allele_string;
        $variant_alleles->{$peptide} = [$codon, $alleles_arr];
    }
}
  
  my %out;
  #now generate all possible codons for each peptide and translate them
  foreach my $peptide (keys %$variant_alleles) {
    my ($codon, $alleles) = @{$variant_alleles->{$peptide}};

    #need to push original nucleotides onto each position
    #so that all possible combinations can be generated
    push @{$alleles->[0]}, substr($codon,0,1);
    push @{$alleles->[1]}, substr($codon,1,1);
    push @{$alleles->[2]}, substr($codon,2,1);

    my %alt_amino_acids;
    foreach my $a1 (@{$alleles->[0]}) {
      substr($codon, 0, 1) = $a1;
      foreach my $a2 (@{$alleles->[1]}) {
        substr($codon, 1, 1) = $a2;
        foreach my $a3 (@{$alleles->[2]}) {
          substr($codon, 2, 1) = $a3;
          my $aa = $codon_table->translate($codon);
	  $aa = '*' if ($codon_table->is_ter_codon($codon));
          #print "$codon translation is $aa\n";
          $alt_amino_acids{$aa} = 1;
        }
      }
    }

    my @aas = keys %alt_amino_acids;
    $out{$peptide} = \@aas;
  }

  return \%out;
}

#creates a new Allele Feature in cdna coordinates

sub convert_to_cdna{
    my ($transcript,$allele) = @_;
    my $new_allele; #new alleleFeature in cdna coordinates
    print "before converting transcript start ", $transcript->cdna_coding_start,"-",$transcript->cdna_coding_end,"\n";
    my $slice = $transcript->slice();
    my $sa = $slice->adaptor();

    $slice = $sa->fetch_by_Feature($transcript);

    $transcript = $transcript->transfer($slice);
    print "transcript start ", $transcript->cdna_coding_start,"-",$transcript->cdna_coding_end,"\n";

    my @coords = $transcript->genomic2cdna($allele->start,
                                 $allele->end,
                                 $allele->strand);
    print Dumper(@coords);    
    exit;
    if((@coords == 1) && ($coords[0]->isa('Bio::EnsEMBL::Mapper::Coordinate'))){
	#copy the AlleleFeature and convert to cdna coords...
	%$new_allele = %$allele;
	bless $new_allele, ref $allele;
	$new_allele->start($coords[0]->start);
	$new_allele->end($coords[0]->end);	
    }
    
    return $new_allele;
}

#
# Classifies a variation which is in the vicinity of a transcript
#
sub type_variation {
  my $tr  = shift;
  my $var = shift;

  if (!$var->isa('Bio::EnsEMBL::Variation::ConsequenceType')){
      throw("Not possible to calculate the consequence type for," ref($var),": Bio::EnsEMBL::Variation::ConsequenceType object expected");
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
