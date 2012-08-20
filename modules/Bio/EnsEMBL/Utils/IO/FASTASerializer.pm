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

Bio::EnsEMBL::Utils::IO::FASTASerializer

=head1 SYNOPSIS

  my $serializer = Bio::EnsEMBL::Utils::IO::FASTASerializer->new($filehandle);
  $serializer->chunk_factor(1000);
  $serializer->line_width(60);
  $serializer->print_Seq($slice);

  $serializer = Bio::EnsEMBL::Utils::IO::FASTASerializer->new($filehandle,
      sub {
        my $slice = shift;
        return "Custom header";
      }
  );

=head1 DESCRIPTION

  Replacement for SeqDumper, making better use of shared code. Outputs FASTA
  format with optional custom header and formatting parameters. Set line_width
  and chunk_factor to dictate buffer size depending on application. A 60kb
  buffer is used by default with a line width of 60 characters.

  Custom headers are set by supplying an anonymous subroutine to new(). Custom
  header code must accept a Slice or Bio::PrimarySeqI compliant object as
  argument and return a string.

  The custom header method can be overridden later through set_custom_header()
  but this is not normally necessary.

=cut

package Bio::EnsEMBL::Utils::IO::FASTASerializer;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception;
use Bio::EnsEMBL::Utils::Scalar qw/assert_ref check_ref/;

use base qw(Bio::EnsEMBL::Utils::IO::Serializer);

=head2 new

  Arg [1]    : Filehandle (optional)
  Arg [2]    : CODEREF subroutine for writing custom headers
  Arg [3]    : [optional] Chunking size (integer)
  Arg [4]    : [optional] Line width (integer)
  Example    : $dumper = Bio::EnsEMBL::Utils::IO::FASTASerializer->new($filehandle,$header_function,1000,60);
  Description: Constructor
               Allows the specification of a custom function for rendering
               header lines.
  Returntype : Bio::EnsEMBL::Utils::IO::FASTASerializer;
  Exceptions : none
  Caller     : general

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  my $filehandle = shift;
  my $header_function = shift;
  my $chunk_factor = shift;
  my $line_width = shift;

  my $self = $class->SUPER::new($filehandle);

  $self->{'header_function'} = $header_function;
  $self->{'line_width'} = ($line_width)? $line_width : 60;
  $self->{'chunk_factor'} = ($chunk_factor)? $chunk_factor : 1000;
  # gives a 60kb buffer by default, increase for higher database and disk efficiency.

  if ( defined($self->{'header_function'}) ) {
    if (ref($self->{'header_function'}) ne "CODE") {
      throw("Custom header function must be an anonymous subroutine when instantiating FASTASerializer");}
  }
  else {
    $self->{'header_function'} = sub {
      my $slice = shift;

      if(check_ref($slice, 'Bio::EnsEMBL::Slice')) {
        my $id     = $slice->seq_region_name;
        my $seqtype  = 'dna';
        my $idtype   = $slice->coord_system->name;
        my $location = $slice->name;

        return "$id $seqtype:$idtype $location";
      }
      else {
        # must be a Bio::Seq , or we're doomed

        return $slice->display_id;
      }
    };

  }

  return $self;
}

=head2 print_metadata

  Arg [1]    : Bio::EnsEMBL::Slice
  Description: Printing header lines into FASTA files. Usually handled
               internally to the serializer.
  Returntype : None
  Caller     : print_Seq
=cut

sub print_metadata {
  my $self = shift;
  my $slice = shift;
  my $fh = $self->{'filehandle'};
  my $function = $self->header_function();
  my $metadata = $function->($slice);
  print $fh '>'.$metadata."\n" or throw "Error writing to file handle: $!";
}

=head2 print_Seq

  Arg [1]    : Bio::EnsEMBL::Slice or other Bio::PrimarySeqI compliant object

  Description: Serializes the slice into FASTA format. Buffering is used
               While other Bioperl PrimarySeqI implementations can be used,
               a custom header function will be required to accommodate it.

  Returntype : None

=cut

sub print_Seq {
  my $self = shift;
  my $slice = shift;
  my $fh = $self->{'filehandle'};

  $self->print_metadata($slice);
  my $width = $self->{line_width};

  # set buffer size
  my $chunk_size = $self->{'chunk_factor'} * $width;

  my $start = 1;
  my $end = $slice->length();

  #chunk the sequence to conserve memory, and print

  my $here = $start;

  while($here <= $end) {
    my $there = $here + $chunk_size - 1;
    $there = $end if($there > $end);
    my $seq = $slice->subseq($here, $there);
    $seq =~ s/(.{1,$width})/$1\n/g;
    print $fh $seq or throw "Error writing to file handle: $!";
    $here = $there + 1;
  }

  if ($slice->length > 0) {$self->{'achieved_something'} = 1;}

}

=head2 line_width

  Arg [1]    : Integer e.g. 60 or 80
  Description: Set and get FASTA format line width. Default is 60
  Returntype : Integer

=cut

sub line_width {
  my $self = shift;
  my $line_width = shift;
  if ($line_width) { $self->{'line_width'} = $line_width };
  return $self->{'line_width'}
}

=head2 chunk_factor
  
  Arg [1]    : Integer e.g. 1000
  Description: Set and get the multiplier used to dictate buffer size
               Chunk factor x line width = buffer size in bases.
  Returntype : Integer

=cut

sub chunk_factor {
  my $self = shift;
  my $chunk_factor = shift;
  if ($chunk_factor) { $self->{'chunk_factor'} = $chunk_factor};
  return $self->{'chunk_factor'}
}

=head2 set_custom_header

  Arg [1]    : CODE reference
  Description: Set the custom header function. Normally this is done at
               construction time, but can be overridden here.
  Example    : $serializer->set_custom_header( sub { return 'New header'});
  Returntype :

=cut

sub set_custom_header {
  my ($self, $new_header_function) = @_;
  $self->header_function($new_header_function);
  return;
}

=head2 header_function

  Arg [1]     : CODE reference (optional)
  Description : Getter/setter for the custom header code
  Example     : $serializer->header_function( sub { return 'New header'});
  Returntype  : CODE

=cut

sub header_function {
  my ($self, $header_function) = @_;
  if($header_function) {
    assert_ref($header_function, 'CODE', 'header_function');
    $self->{header_function} = $header_function;
  }
  return $self->{header_function};
}

1;
