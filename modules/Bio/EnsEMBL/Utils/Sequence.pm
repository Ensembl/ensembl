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

Bio::EnsEMBL::Utils::Sequence - Utility functions for sequences

=head1 SYNOPSIS

  use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp expand);

  my $seq = 'ACTTTAAAGGCTATCCCAATATG';

  print "my sequence = $seq\n";

  reverse_comp( \$seq );

  print "my reverse comp = $seq\n";

  my $compressed_seq = '(AC)3';

  print "my expanded seq is = expand($compressed_seq)";

=head1 METHODS

=cut


package Bio::EnsEMBL::Utils::Sequence;

use strict;
use warnings;

use Exporter;

use vars qw(@ISA @EXPORT_OK);

@ISA = qw(Exporter);

@EXPORT_OK = qw(&reverse_comp &expand);


=head2 reverse_comp

  Arg [1]    : reference to a string $seqref
  Example    : use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);

               $seq = 'ACCTGAA';
               reverse_comp(\$seq);
               print $seq;

  Description: Does an in place reverse compliment of a passed in string
               reference.  The string is passed by reference
               rather than by value for memory efficiency.
  Returntype : none
  Exceptions : none
  Caller     : SequenceAdaptor, SliceAdaptor

=cut

sub reverse_comp {
  my $seqref = shift;

  $$seqref = reverse( $$seqref );
  $$seqref =~
    tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;

  return;
}

=head2 expand

  Arg [1]    : reference to a string $seqref
  Example    : use Bio::EnsEMBL::Utils::Sequence qw(expand);

               $seq = '(AC)3';
               expand(\$seq);
               print $seq;
              

  Description: Expands a genomic sequence. The string is passed by reference
               rather than by value for memory efficiency.
  Returntype : none
  Exceptions : none
  Caller     : SequenceAdaptor, SliceAdaptor

=cut

sub expand {
  my $seq_ref = shift;
     $$seq_ref =~ s/(\w*)\((\w+)\)(\d+)/$1.$2 x $3/eg  if ($$seq_ref =~ /\(/);#expressions with parenthesis, expand the alleles			    
  return;
}


1;
