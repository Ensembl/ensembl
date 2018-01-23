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

Bio::EnsEMBL:SeqEdit - A class representing a post transcriptional edit to a
sequence.

=head1 SYNOPSIS

  use Bio::EnsEMBL::SeqEdit;
  use Bio::EnsEMBL::Attribute;

  # construct a SeqEdit object using a Transcript attribute

  ($attribute) = @{ $translation->get_all_Attributes('_rna_edit') };

  $seq_edit = Bio::EnsEMBL::SeqEdit( -ATTRIB => $attribute );

  print $seq_edit->start(),   "\n";
  print $seq_edit->end(),     "\n";
  print $seq_edit->alt_seq(), "\n";

  # apply the edit to some sequence
  $seq = $transcript->spliced_seq();
  print "Before modifiction: $seq\n";

  $seq_edit->apply_edit( \$seq );
  print "After modification: $seq\n";

  # construct an attribute object from a SeqEdit and add it to a
  # translation

  $seq_edit = Bio::EnsEMBL::SeqEdit->new(
    -CODE    => '_selenocysteine',
    -NAME    => 'Selenocysteine',
    -DESC    => 'Selenocysteine',
    -START   => 10,
    -END     => 10,
    -ALT_SEQ => 'U'
  );

  $attribute = $seq_edit->get_Attribute();
  $translation->add_Attributes($attribute);

=head1 DESCRIPTION

This is a class used to represent post transcriptional
modifications to sequences.  SeqEdit objects are stored as ordinary
Bio::EnsEMBL::Attributes with a parseable value and can be used to
represent RNA editing, selenocysteines etc.

Also see B<Bio::EnsEMBL::Attribute>

=head1 METHODS

=cut

package Bio::EnsEMBL::SeqEdit;

use strict;
use warnings;

use Bio::EnsEMBL::Attribute;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw);


=head2 new

  Arg [-ATTRIB] : Bio::EnsEMBL::Attribute
                  Constructs a new SeqEdit from an Attribute.
                  Can only be provided if no other constructor arguments
                  are provided.
  Arg [-START]       : The start position of the edit.
  Arg [-END]         : The end position of the edit.
  Arg [-ALT_SEQ]     : The alternate sequence
  Arg [-CODE]        : A code for this SeqEdit
  Arg [-NAME]        : A name for this SeqEdit
  Arg [-DESCRIPTION] : Arg passed to superclass constructor
  Example    : my $sea = Bio::EnsEMBL::SeqEdit->new(-ATTRIB => $attrib);
               my $sea = Bio::EnsEMBL::SeqEdit->new
                             (-START => 10,
                              -END   => 12,
                              -ALT_SEQ => 'ACG',
                              -CODE    => '_rna_edit',
                              -NAME    => 'RNA Edit',
                              -DESCRIPTION => 'RNA edit');
  Description: Constructs a SeqEdit representing a single edit to a
               sequence, such as an rna modification or a selenocysteine.
  Returntype : Bio::EnsEMBL::SeqEdit
  Exceptions : throws if attribute set and other args aswell
               throws if start and end not set correctly of attribure not set
  Caller     : general
  Status     : Stable

=cut

sub new {
  my $class = shift;

  my ($attrib, $start, $end, $alt_seq, $name, $desc, $code) =
    rearrange([qw(ATTRIB START END ALT_SEQ NAME DESCRIPTION CODE)], @_);

  my $self;

  if($attrib) {
    if(defined($start) || defined($end) || defined($alt_seq) ||
       defined($name)  || defined($desc) || defined($code)) {
      throw("Cannot specify -ATTRIB argument with additional arguments.");
    }

    if(!ref($attrib) || !$attrib->isa('Bio::EnsEMBL::Attribute')) {
      throw('Bio::EnsEMBL::Attribute argument expected.');
    }

    ($start, $end, $alt_seq) = split(/\s+/, $attrib->value());

    if($start !~ /\d+/ || $end !~ /\d+/) {
      throw('Could not parse value of attribute: '.$attrib->value());
    }

    $name  = $attrib->name();
    $code  = $attrib->code();
    $desc = $attrib->description();


  }

  if(defined($end) && defined($start) && $start > $end+1) {
    throw("start must be less than or equal to end + 1");
  }

  if(defined($start) && $start < 1) {
    throw("start must be greater than or equal to 1");
  }

  if(defined($end) && $end < 0) {
    throw("end must be greater than or equal to 0");
  }

  $alt_seq ||= '';

  return bless {'start'        => $start,
                'end'          => $end,
                'alt_seq'      => $alt_seq,
                'description'  => $desc,
                'name'         => $name,
                'code'         => $code}, $class;
}



=head2 start

  Arg [1]    : (optional) int $start - the new start position
  Example    : $start = $se_attrib->start();
  Description: Getter/Setter for the start position of the region replaced
               by the alt_seq.

               Coordinates are inclusive and one-based, which means that
               inserts are unusually represented by a start 1bp higher than
               the end.

               E.g. start = 1, end = 1 is a replacement of the first base but 
               start = 1, end = 0 is an insert BEFORE the first base.
  Returntype : int
  Exceptions : none
  Caller     : Transcript, Translation
  Status     : Stable

=cut

sub start {
  my $self = shift;

  if(@_) {
    my $start = shift;
    if(defined($start) && $start < 1) {
      throw("start must be greater than or equal to 1");
    }
    $self->{'start'} = $start;
  }

  return $self->{'start'};
}


=head2 end

  Arg [1]    : (optional) int $end - the new end position
  Example    : $end = $se_attrib->end();
  Description: Getter/Setter for the end position of the region replaced
               by the alt_seq.

               Coordinates are inclusive and one-based, which means that
               inserts are unusually represented by a start 1bp higher than
               the end.

               E.g. start = 1, end = 1 is a replacement of the first base but
               start = 1, end = 0 is an insert BEFORE the first base.
  Returntype : int
  Exceptions : throws if end  <= 0
  Caller     : Transcript, Translation
  Status     : Stable

=cut

sub end {
  my $self = shift;

  if(@_) {
    my $end = shift;
    if(defined($end) && $end < 0) {
      throw("end must be greater than or equal to 0");
    }
    $self->{'end'} = $end;
  }

  return $self->{'end'};
}


=head2 alt_seq

  Arg [1]    : (optional) string $alt_seq
  Example    : my $alt_seq = $se_attrib->alt_seq();
  Description: Getter/Setter for the replacement sequence used by this edit.
               The sequence may either be a string of amino acids or
               nucleotides depending on the context in which this edit is
               used.

               In the case of a deletion the replacement sequence is an empty
               string.
  Returntype : string
  Exceptions : none
  Caller     : Transcript, Translation
  Status     : Stable

=cut

sub alt_seq {
  my $self = shift;
  $self->{'alt_seq'} = shift || '' if(@_);
  return $self->{'alt_seq'};
}


=head2 length_diff

  Arg [1]    : none
  Example    : my $diff = $sea->length_diff();
  Description: Returns the difference in length caused by applying this
               edit to a sequence.  This may be be negative (deletion),
               positive (insertion) or 0 (replacement).

               If either start or end are not defined 0 is returned.
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub length_diff {
  my $self = shift;

  return 0 if(!defined($self->{'end'}) || !defined($self->{'start'}));

  return length($self->{'alt_seq'}) - ($self->{'end'} - $self->{'start'} + 1);
}



=head2 name

  Arg [1]    : (optional) string $name
  Example    : my $name = $seqedit->name();
  Description: Getter/Setter for the name of this SeqEdit
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub name {
  my $self = shift;
  $self->{'name'} = shift if(@_);
  return $self->{'name'};
}




=head2 code

  Arg [1]    : (optional) string $code
  Example    : my $code = $seqedit->code();
  Description: Getter/Setter for the code of this SeqEdit
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub code {
  my $self = shift;
  $self->{'code'} = shift if(@_);
  return $self->{'code'};
}



=head2 description

  Arg [1]    : (optional) string $desc
  Example    : my $desc = $seqedit->description();
  Description: Getter/Setter for the description of this SeqEdit
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub description {
  my $self = shift;
  $self->{'description'} = shift if(@_);
  return $self->{'description'};
}



=head2 get_Attribute

  Arg [1]    : none
  Example    : my $attrib = $seqedit->get_Attribute();
               $transcript->add_Attributes($attrib);
  Description: Converts a SeqEdit object into an Attribute object.  This
               allows the SeqEdit to be stored as any other attribute in the
               ensembl database.  The start/end and alt_seq properties
               should be set before calling this method.
  Returntype : Bio::EnsEMBL::Attribute
  Exceptions : warning if start/end or alt_seq properties are not defined
  Caller     : general
  Status     : Stable

=cut

sub get_Attribute {
  my $self = shift;

  my $start = $self->start();
  my $end  = $self->end();
  my $alt_seq = $self->alt_seq();

  my $value;

  if(defined($start) && defined($end) && defined($alt_seq)) {
    $value = join(' ', $start, $end, $alt_seq);
  } else {
    warning('Attribute value cannot be created unless start, end and alt_seq' .
            'properties are defined');
    $value = '';
  }

  return Bio::EnsEMBL::Attribute->new(-CODE  => $self->code(),
                                      -VALUE => $value,
                                      -NAME  => $self->name(),
                                      -DESCRIPTION => $self->description());
}


=head2 apply_edit

  Arg [1]    : reference to string $seqref
  Example    : $sequence = 'ACTGAATATTTAAGGCA';
               $seqedit->apply_edit(\$sequence);
               print $sequence, "\n";
  Description: Applies this edit directly to a sequence which is
               passed by reference.  The coordinates of this SeqEdit
               are assumed to be relative to the start of the sequence
               argument.
               If either the start or end of this SeqEdit are not defined
               this function will not do anything to the passed sequence.
  Returntype : reference to the same sequence that was passed in
  Exceptions : none
  Caller     : Transcript, Translation
  Status     : Stable

=cut

sub apply_edit {
  my $self   = shift;
  my $seqref = shift;

  if(ref($seqref) ne 'SCALAR') {
    throw("Reference to scalar argument expected");
  }

  if(!defined($self->{'start'}) || !defined($self->{'end'})) {
    return $seqref;
  }

  my $len = $self->{'end'} - $self->{'start'} + 1;
  substr($$seqref, $self->{'start'} - 1, $len) = $self->{'alt_seq'};

  return $seqref;
}


1;
