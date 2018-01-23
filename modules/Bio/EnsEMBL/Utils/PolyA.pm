=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Utils::PolyA

=head1 SYNOPSIS

  my $seq;    # a Bio::Seq object
  my $polyA = Bio::EnsEMBL::Utils::PolyA->new();

  # returns a new Bio::Seq object with the trimmed sequence
  my $trimmed_seq = $polyA->clip($seq);

  # cat put Ns in the place of the polyA/polyT tail
  my $masked_seq = $polyA->mask($seq);

  # can put in lower case the polyA/polyT using any flag:
  my $softmasked_seq = $poly->mask( $seq, 'soft' );

=head1 DESCRIPTION

  It reads a Bio::Seq object, it first finds out whether it has a
  polyA or a polyT and then performs one operation in the seq string:
  clipping, masking or softmasking.  It then returns a new Bio::Seq
  object with the new sequence.

=head1 METHODS

=cut

package Bio::EnsEMBL::Utils::PolyA;

use Bio::Seq;
use vars qw(@ISA);

use strict;


=head2 new

=cut

sub new{
  my ($class) = @_;
   if (ref($class)){
    $class = ref($class);
  }
  my $self = {};
  bless($self,$class);
  
  return $self;
}


############################################################

sub clip{
  my ($self, $bioseq) = @_;

  # print STDERR "past a $bioseq\n";
  my $seq = $bioseq->seq;
  $self->_clip(1);
  $self->_mask(0);
  $self->_softmask(0);
  my $new_seq = $self->_find_polyA($seq);
  my $cdna = Bio::Seq->new();
  if (length($new_seq) > 0){
    $cdna->seq($new_seq);
  } 
  else {
    print "While clipping the the polyA tail, sequence ".$bioseq->display_id." totally disappeared.\n";
    print "Returning undef\n";
    return undef;
  }
  $cdna->display_id( $bioseq->display_id );
  $cdna->desc( $bioseq->desc );
  
  return $cdna;
}

############################################################

sub mask{
  my ($self, $bioseq, $flag ) = @_;
  
  my $seq = $bioseq->seq;
  $self->_clip(0);
  if ( $flag ){
    $self->_mask(0);
    $self->_softmask(1);
  }
  else{
    $self->_mask(1);
    $self->_softmask(0);
  }
  my $new_seq = $self->_find_polyA($seq);
  my $cdna = new Bio::Seq;
  $cdna->seq($new_seq);
  $cdna->display_id( $bioseq->display_id );
  $cdna->desc( $bioseq->desc );
  
  return $cdna;
}

############################################################




############################################################

sub _find_polyA{
  my ($self, $seq) = @_;
  my $new_seq;
  my $length = length($seq);
  
  # is it a polyA or polyT?
  my $check_polyT = substr( $seq, 0, 6 );
  
  my $check_polyA = substr( $seq, -6 );
  
  my $t_count = $check_polyT =~ tr/Tt//;
  my $a_count = $check_polyA =~ tr/Aa//;
  
  #### polyA ####
  if ( $a_count >= 5 && $a_count > $t_count ){
    
    # we calculate the number of bases we want to chop
    my $length_to_mask = 0;
    
    # we start with 3 bases
    my ($piece, $count ) = (3,0);
    
    # count also the number of Ns, consider the Ns as potential As
    my $n_count = 0;

    # take 3 by 3 bases from the end
    while( $length_to_mask < $length ){
      my $chunk  = substr( $seq, ($length - ($length_to_mask + 3)), $piece);
      $count   = $chunk =~ tr/Aa//;
      $n_count = $chunk =~ tr/Nn//;
      if ( ($count + $n_count) >= 2*( $piece )/3 ){
  $length_to_mask += 3;
      }
      else{
  last;
      }
    }
    
    if ( $length_to_mask > 0 ){
      # do not mask the last base if it is not an A:
      my $last_base        = substr( $seq, ( $length - $length_to_mask    ), 1);
      my $previous_to_last = substr( $seq, ( $length - $length_to_mask - 1), 1);
      if ( !( $last_base eq 'A' || $last_base eq 'a') ){
  $length_to_mask--;
      }
      elsif( $previous_to_last eq 'A' || $previous_to_last eq 'a' ){
  $length_to_mask++;
      }
      my $clipped_seq = substr( $seq, 0, $length - $length_to_mask );
      my $mask;
      if ( $self->_clip ){
  $mask = "";
      }
      elsif( $self->_mask ){
  $mask    = "N" x ($length_to_mask);
      }
      elsif ( $self->_softmask ){
  $mask = lc substr( $seq, ( $length - $length_to_mask ) );
      }
      $new_seq =  $clipped_seq . $mask;
    }
    else{
      $new_seq = $seq;
    }  
  }
  #### polyT ####
  elsif( $t_count >=5 && $t_count > $a_count ){
    
    # calculate the number of bases to chop
    my $length_to_mask = -3;
    
    # we start with 3 bases:
    my ($piece, $count) = (3,3);
    
    # count also the number of Ns, consider the Ns as potential As
    my $n_count = 0;
    
    # take 3 by 3 bases from the beginning
    while ( $length_to_mask < $length ){
      my $chunk = substr( $seq, $length_to_mask + 3, $piece );
      #print STDERR "length to mask: $length_to_mask\n";
      #print "chunk: $chunk\n";
      $count = $chunk =~ tr/Tt//;
      $n_count = $chunk =~ tr/Nn//;
      if ( ($count+$n_count)  >= 2*( $piece )/3 ){
  $length_to_mask +=3;
      }
      else{
  last;
  
      }
    }
    if ( $length_to_mask >= 0 ){
      # do not chop the last base if it is not a A:
      #print STDERR "clipping sequence $seq\n";
      my $last_base        = substr( $seq, ( $length_to_mask + 3 - 1 ), 1 );
      my $previous_to_last = substr( $seq, ( $length_to_mask + 3     ), 1 );
      if ( !( $last_base eq 'T' || $last_base eq 't' ) ){
  $length_to_mask--;
      }
      elsif( $previous_to_last eq 'T' || $previous_to_last eq 't' ){
  $length_to_mask++;
      }
      my $clipped_seq = substr( $seq, $length_to_mask + 3);
      my $mask;
      if ( $self->_clip ){
  $mask = "";
      }
      elsif( $self->_mask ){
  $mask        = "N" x ($length_to_mask+3);
      }
      elsif ($self->_softmask){
  $mask = lc substr( $seq, 0, ($length_to_mask + 3) );
      }
      $new_seq     = $mask.$clipped_seq;
    }
    else{
      $new_seq = $seq;
    }
  }
  else{
    # we cannot be sure of what it is
    # do not clip
    $new_seq = $seq;
  }

  return $new_seq;
}

############################################################

sub _mask{
  my ($self,$mask) = @_;
  if (defined($mask)){
    $self->{_mask} = $mask;
  }
  return $self->{_mask};
}

############################################################

sub _clip{
  my ($self,$clip) = @_;
  if (defined($clip)){
    $self->{_clip} = $clip;
  }
  return $self->{_clip};
}

############################################################

sub _softmask{
  my ($self,$softmask) = @_;
  if (defined($softmask)){
    $self->{_softmask} = $softmask;
  }
  return $self->{_softmask};
}

############################################################


sub has_polyA_track{
  my ($self, $seq) = @_;
  my $new_seq;
  my $length = length($seq);
  
  # is it a polyA or polyT?
  my $check_polyT = substr( $seq, 0, 10 );
  
  my $check_polyA = substr( $seq, -10 );
  
  print STDERR "polyA: $check_polyA\n";

  my $t_count = $check_polyT =~ tr/Tt//;
  my $a_count = $check_polyA =~ tr/Aa//;
  
  ## testing with this short cut
  if ( $a_count >=7 || $t_count >=7 ){
    return 1;
  }
  else{
    return 0;
  }
}


################
1;

