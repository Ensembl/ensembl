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

Bio::EnsEMBL::PaddedSlice

=head1 DESCRIPTION

Used when dumping Slices which represet a portion of the sequence region
they map to e.g. the first section of human Y. The code will return N
as sequence if an attempt is made to retrieve sequence not covered by the
Slice given. This makes the code very memory efficient if sequence dumping
is carried out using C<subseq()> calls.

=head1 METHODS

=cut

package Bio::EnsEMBL::PaddedSlice;

use Bio::EnsEMBL::Utils::Argument qw/rearrange/;
use Bio::EnsEMBL::Utils::Scalar qw/assert_ref assert_strand/;
use base qw/Bio::EnsEMBL::Utils::Proxy/;

=head2 new()

  Arg [SLICE] : The Slice to proxy  
  Example     : my $newobj = Bio::EnsEMBL::PaddedSlice->new($myobj);
  Description : Provides a new instance of a padded slice
  Returntype  : Bio::EnsEMBL::PaddedSlice
  Exceptions  : None 
  Caller      : public
  Status      : -

=cut

sub new {
  my ($class, @args) = @_;
  my ($slice) = rearrange([qw/slice/], @args);
  return $class->SUPER::new($slice);
}

=head2 start()

  Example     : $slice->start();
  Description : Always returns 1 since all padded slices start at 1
  Returntype  : Int
  Exceptions  : None 
  Caller      : public
  Status      : -

=cut

sub start {
  my ($self) = @_;
  return 1;
}

=head2 end()
  
  Example     : $slice->end();
  Description : Always returns the backing slice sequence region length
  Returntype  : Int
  Exceptions  : None 
  Caller      : public
  Status      : -

=cut

sub end {
  my ($self) = @_;
  return $self->seq_region_length();
}

=head2 length()
  
  Example     : $slice->length();
  Description : Delegates to C<end()>
  Returntype  : Int
  Exceptions  : None 
  Caller      : public
  Status      : -

=cut

sub length {
  my ($self) = @_;
  return $self->end();
}

=head2 seq()
  
  Example     : my $seq = $slice->seq()
  Description : Returns the entire sequence of the backing slice but padded
                with N's at the beginning and the end of the slice where
                applicable
  Returntype  : Scalar string
  Exceptions  : None 
  Caller      : public
  Status      : -

=cut

sub seq {
  my ($self) = @_;
  my $parent_slice = $self->__proxy();
  my $pad_start = 'N' x ( $parent_slice->start() - 1 );
  my $pad_end   = 'N' x ( $parent_slice->seq_region_length() - $parent_slice->end() );
  my $seq = $parent_slice->seq();
  return $pad_start . $seq . $pad_end;
}

=head2 subseq()
  
  Arg [1]     : Int; start position of the subslice
  Arg [2]     : Int; end position of the subslice
  Arg [3]     : Int; strand of the subslice  
  Example     : my $subseq = $slice->subseq(1, 1_000_000);
  Description : Returns a portion of the sequence padded with N's if required 
  Returntype  : Scalar string
  Exceptions  : None 
  Caller      : public
  Status      : -

=cut

sub subseq {
  my ( $self, $start, $end, $strand ) = @_;

  if ( $end+1 < $start ) {
    throw("End coord + 1 is less than start coord");
  }
  
  return '' if( $start == $end + 1);
  
  $strand = 1 unless(defined $strand);
  assert_strand($strand, 'strand');
  
  my $parent_slice = $self->__proxy();
    
  #Coords relative to the SeqRegion i.e. huge
  my $parent_start = $parent_slice->start();
  my $parent_end = $parent_slice->end();
  
  #Return if we were upstream of overlap
  if($start < $parent_start && $end < $parent_start) {
    return N x (( $end - $start )+1);
  }
  #Return if we were downstream of overlap
  if($start > $parent_end && $end > $parent_end) {
    return N x (( $end - $start )+1);
  }
  
  my $prefix  = '';
  my $suffix = '';
  my $subslice_start = ($start - $parent_start)+1;
  my $subslice_end = ($end - $parent_start) + 1;
  if($start < $parent_start) {
    $prefix = N x ($parent_start - $start);
    $subslice_start = 1;
  }
  if($end > $parent_end) {
    $suffix = N x ($end - $parent_end);
    $subslice_end = (($parent_end - $parent_start)+1);
  }
  
  my $subseq = $parent_slice->subseq($subslice_start, $subslice_end, $strand);
  
  return $prefix . $subseq . $suffix;
}

=head2 subseq()
  
  Arg [1]     : Int; start position of the subslice
  Arg [2]     : Int; end position of the subslice
  Arg [3]     : Int; strand of the subslice  
  Example     : my $subseq = $slice->subseq(1, 1_000_000);
  Description : Returns a portion of the sequence padded with N's if required 
  Returntype  : Scalar string
  Exceptions  : None 
  Caller      : public
  Status      : -

=cut

sub sub_Slice {
  die "Unsupported";
}


=head2 __resolver()

  Description : Delegates all non-overriden actions onto the backing slice 
  Returntype  : CodeRef
  Exceptions  : None 
  Caller      : public
  Status      : -

=cut

sub __resolver {
  my ($self, $package_name, $method) = @_;
  return sub {
    my ($local_self, @args) = @_;
    return $local_self->__proxy()->$method(@args);
  };
}


1;
