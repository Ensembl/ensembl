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

TranscriptSNPs - A utility class used to obtain information about the
relationships between a transcript and SNPs

=head1 SYNOPSIS

  use Bio::EnsEMBL::Utils::TranscriptSNPs;

  # get and type all snps in the region of the transcript

  %snps = %{
    Bio::EnsEMBL::Utils::TranscriptSNPs::get_all_SNPs( $transcript,
      $flanking ) };

  # get all snps overlapping the transcript in cdna coordinates

  %snps =
    %{ Bio::EnsEMBL::Utils::TranscriptSNPs::get_all_cdna_SNPs(
      $transcript) };

  # get the peptide variations caused by a set of SNPs

  %variations = %{
    Bio::EnsEMBL::Utils::TranscriptSNPs::get_all_peptide_variations(
      $transcript, $snps ) };

=head1 DESCRIPTION

This is a utility class which can be used to get snps associated with a
transcript, and to determine the amino acid changes caused by the SNPs

=head1 METHODS

=cut

package Bio::EnsEMBL::Utils::TranscriptSNPs;

use strict;
use warnings;
no warnings 'uninitialized';



use Bio::EnsEMBL::Utils::Exception qw(throw warning);


=head2 get_all_peptide_variations

  Arg [1]    : $transcript the transcript to obtain the peptide variations for
  Arg [2]    : $snps listref of coding snps in cdna coordinates
  Example    : $pep_hash = get_all_peptide_variations($transcript, \@snps);
  Description: Takes a list of coding snps on this transcript in
               which are in cdna coordinates and returns a hash with peptide
               coordinate keys and listrefs of alternative amino acids as
               values.  The SNPs must additionally have a strand of 1 for the
               sake of simplicity.  Normally these could be generated using the
               get_all_cdna_SNPs method.

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
  my $snps = shift;

  if(!ref($transcript) || !$transcript->isa('Bio::EnsEMBL::Transcript')) {
    throw('Bio::EnsEMBL::Transcript argument is required.');
  }

  if(!ref($snps) eq 'ARRAY') {
    throw('Reference to a list of Bio::EnsEMBL::SNP objects is required');
  }

  my $codon_table = Bio::Tools::CodonTable->new;
  my $codon_length = 3;
  my $cdna = $transcript->spliced_seq;

  my $variant_alleles;
  my $translation_start = $transcript->cdna_coding_start;
  foreach my $snp (@$snps) {
    #skip variations not on a single base
    next if ($snp->start != $snp->end);

    my $start = $snp->start;
    my $strand = $snp->strand;

    #calculate offset of the nucleotide from codon start (0|1|2)
    my $codon_pos = ($start - $translation_start) % $codon_length;

    #calculate the peptide coordinate of the snp
    my $peptide = ($start - $translation_start +
		   ($codon_length - $codon_pos)) / $codon_length;
	
	# skip this SNP if it falls in a partial codon
	next if $start - $codon_pos + $codon_length > length($cdna);

    #retrieve the codon
    my $codon = substr($cdna, $start - $codon_pos-1, $codon_length);

    #store each alternative allele by its location in the peptide
    my @alleles = split(/\/|\|/, lc($snp->allele_string));
    #my @alleles = split(/\/|\|/, lc($snp->alleles));

    foreach my $allele (@alleles) {
      next if $allele eq '-';       #skip deletions
      next if CORE::length($allele) != 1; #skip insertions

      #get_all_cdna_SNPs always gives strand of 1 now
      #if($strand == -1) {
      #  #complement the allele if the snp is on the reverse strand
      #  $allele =~ 
      #  tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
      #}

      #create a data structure of variant alleles sorted by both their
      #peptide position and their position within the peptides codon
      $variant_alleles ||= {};
      if(exists $variant_alleles->{$peptide}) {
        my $alleles_arr = $variant_alleles->{$peptide}->[1];
        push @{$alleles_arr->[$codon_pos]}, $allele;
      } else {
        #create a list of 3 lists (one list for each codon position)
        my $alleles_arr = [[],[],[]];
        push @{$alleles_arr->[$codon_pos]}, $allele;
        $variant_alleles->{$peptide} = [$codon, $alleles_arr];
      }
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



=head2 get_all_SNPs

  Arg [1]    : Bio::EnsEMBL::Transcript $transcript

  Arg [2]    : (optional) int $flanking
               The number of basepairs of transcript flanking sequence to
               retrieve snps from (default 0)
  Arg [3]    : $source type of database source (dbSNP, Glovar)
  Example    : $snp_hashref = get_all_transcript_SNPs($transcript)
  Description: Retrieves all snps found within the region of the 
               provided transcript
               The snps are returned in a hash with keys corresponding
               to the region the snp was found in.  Possible keys are:
               'three prime UTR', 'five prime UTR', 'coding', 'intronic',
               'three prime flanking', 'five prime flanking'
               If no flanking argument is provided no flanking snps will be
               obtained.
               The listrefs which are the values of the returned hash
               contain snps in coordinates of the transcript region 
               (i.e. first base = first base of the first exon on the
               postive strand - flanking bases + 1)

               Multiple base variations and inserts/deletes are discarded
               by this method and not used.

  Returntype : hasref with string keys and listrefs of Bio::EnsEMBL::SNPs for
               values
  Exceptions : none
  Caller     : general

=cut

sub get_all_SNPs {
  my $transcript = shift;
  my $flanking = shift || 0;
  my $source = shift;

  if(!ref($transcript) || !$transcript->isa('Bio::EnsEMBL::Transcript')) {
    throw('Bio::EnsEMBL::Transcript argument required.');
  }

  my $slice = $transcript->slice();

  if(!$slice) {
    warning("Cannot obtain SNPs for transcript without attached Slice.");
    return {};
  }

  my $sa = $slice->adaptor();

  if(!$sa) {
    warning('Cannot obtain SNPs for transcript unless attached slice ' .
            'has attached adaptor');
    return {};
  }

  my %snp_hash;

  # retrieve slice in the region of the transcript
  $slice = $sa->fetch_by_Feature($transcript, $flanking );

  # copy transcript, to work in coord system we are interested in
  $transcript = $transcript->transfer( $slice );

  # get all snps in the transcript region
  my $snps;
  if ($source eq 'glovar') {
    $snps = $slice->get_all_ExternalFeatures('GlovarSNP');
  }
  elsif ($source eq 'variation') {
    $snps = $slice->get_all_VariationFeatures;
  }
  else {
    $snps = $slice->get_all_SNPs;   # dont need once use new snp api (i think)
  }

  my $trans_start  = $flanking + 1;
  my $trans_end    = $slice->length - $flanking;
  my $trans_strand = $transcript->get_all_Exons->[0]->strand;

  # classify each snp
  foreach my $snp (@$snps) {
    my $key;

    if(($trans_strand == 1 && $snp->end < $trans_start) ||
       ($trans_strand == -1 && $snp->start > $trans_end)) {
      # this snp is upstream from the transcript
      $key = 'five prime flanking';
    }

    elsif(($trans_strand == 1 && $snp->start > $trans_end) ||
	  ($trans_strand == -1 && $snp->start < $trans_start)) {
      # this snp is downstream from the transcript
      $key = 'three prime flanking';
    }

    else {
      #snp is inside transcript region check if it overlaps an exon
      foreach my $e (@{$transcript->get_all_Exons}) {
        if($snp->end >= $e->start && $snp->start <= $e->end) {
          # this snp is in an exon
	  
          if(($trans_strand == 1 && 
              $snp->end < $transcript->coding_region_start) ||
             ($trans_strand == -1 && 
              $snp->start > $transcript->coding_region_end)) {
            # this snp is in the 5' UTR
            $key = 'five prime UTR';
          }

          elsif(($trans_strand == 1 && 
                 $snp->start > $transcript->coding_region_end)||
                ($trans_strand == -1 && 
                 $snp->end < $transcript->coding_region_start)) {
            # this snp is in the 3' UTR
            $key = 'three prime UTR';
          }

          else {
            # snp is coding
            $key = 'coding';
          }
          last;
        }
      }
      unless($key) {
        # snp was not in an exon and is therefore intronic
        $key = 'intronic';
      }
    }

    unless($key) {
      #warning('SNP could not be mapped. In/Dels not supported yet...');
      next;
    }

    if(exists $snp_hash{$key}) {
      push @{$snp_hash{$key}}, $snp;
    } 
    else {
      $snp_hash{$key} = [$snp];
    }
  }

  return \%snp_hash;
}



=head2 get_all_cdna_SNPs

  Arg [1]    : Bio::EnsEMBL::Transcript $transcript
  Arg [2]    : $source type of database source (dbSNP, Glovar)
  Example    : $cdna_snp_hasref = $transcript->get_all_cdna_SNPs;
  Description: Retrieves all snps found within exons of the provided
               transcript.
               The snps are returned in a hash with three keys corresponding
               to the region the snp was found in.  Valid keys are:
               'three prime UTR', 'five prime UTR', 'coding'
               The listrefs which are the values of the returned hash
               contain snps in CDNA coordinates.

               Multiple base variations and insertions/deletions are not
               used by this function and are discarded.
  Returntype : hasref with string keys and listrefs of Bio::EnsEMBL::SNPs for 
               values
  Exceptions : none
  Caller     : general

=cut

sub get_all_cdna_SNPs {
  my ($transcript, $source) = @_;

  #retrieve all of the snps from this transcript
  my $all_snps = get_all_SNPs($transcript, 0, $source);
  my %snp_hash;

  my @cdna_types = ('three prime UTR', 'five prime UTR','coding');

  my $slice = $transcript->slice();
  my $sa    = $slice->adaptor();

  $slice = $sa->fetch_by_Feature($transcript);

  # copy transcript in order to work in coord system of interest
  $transcript = $transcript->transfer($slice);

  foreach my $type (@cdna_types) {
    $snp_hash{$type} = [];
    foreach my $snp (@{$all_snps->{$type}}) {
      my @coords = $transcript->genomic2cdna($snp->start, $snp->end,
                                             $snp->strand);

      #skip snps that don't map cleanly (possibly an indel...)
      if(scalar(@coords) != 1) {
        #warning("snp of type $type does not map cleanly\n");
        next;
      }

      my ($coord) = @coords;

      unless($coord->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
        #warning("snp of type $type maps to gap\n");
        next;
      }

      my $alleles;
      my $ambicode;

      # get alleles and ambig_code (with fallback to old snp API)
      $alleles = $snp->allele_string || $snp->{'alleles'};
      $ambicode = $snp->ambig_code || $snp->{'_ambiguity_code'};

      #we arbitrarily put the SNP on the +ve strand because it is easier to
      #work with in the webcode
      if($coord->strand == -1) {
        $alleles =~
         tr/acgthvmrdbkynwsACGTDBKYHVMRNWS\//tgcadbkyhvmrnwsTGCAHVMRDBKYNWS\//;
        $ambicode =~
         tr/acgthvmrdbkynwsACGTDBKYHVMRNWS\//tgcadbkyhvmrnwsTGCAHVMRDBKYNWS\//;
      }
      #copy the snp and convert to cdna coords...
      my $new_snp;
      %$new_snp = %$snp;
      bless $new_snp, ref $snp;
      $new_snp->start($coord->start);
      $new_snp->end($coord->end);
      $new_snp->strand(1);
      $new_snp->allele_string($alleles);
      $new_snp->ambig_code($ambicode);
      push @{$snp_hash{$type}}, $new_snp;
    }
  }

  return \%snp_hash;
}




1;
