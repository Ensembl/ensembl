# Creator: Arne Stabenau <stabenau@ebi.ac.uk>
# Date of creation: 07.04.2000
# Last modified : 16.02.2003 by James Smith
#
# Copyright EMBL-EBI/WTSI 2000
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Chromosome

=head1 SYNOPSIS

=head1 DESCRIPTION

The chromosome object is deprecated.  Use Bio::EnsEMBL::Slice instead.

=head1 CONTACT

Post questions to the EnsEMBL developer mailing list: <ensembl-dev@ebi.ac.uk>

=cut


package Bio::EnsEMBL::Chromosome;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Utils::Exception qw(deprecate);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

use Bio::EnsEMBL::Slice;

@ISA = qw( Bio::EnsEMBL::Slice );


=head2 new

  Description: The chromosome class is deprecated. Bio::EnsEMBL::Slice class
               should be used instead.

=cut

sub new {
  my $caller = shift;

  deprecate("The Bio::EnsEMBL::Chromosome class is deprecated." .
            "The Bio::EnsEMBL::Slice class should be used instead");

  my ($adaptor, $chr_name) = rearrange(['ADAPTOR', 'CHR_NAME'], @_);

  if($chr_name) {
    if($adaptor) {
      my $csa = $adaptor->db->get_CoordSystemAdaptor();
      my ($top_cs) = @{$csa->fetch_all()};
      my $chr = $adaptor->fetch_by_region($top_cs->name(),$chr_name,
                                         undef,undef,undef,$top_cs->version());
      bless $chr, 'Bio::EnsEMBL::Chromosome';
      return $chr;
    } else {
      return $caller->SUPER::new(@_,
                                 -ADAPTOR => $adaptor,
                                 -SEQ_REGION_NAME => $chr_name);
    }
  }

  my $self = $caller->SUPER::new(@_);
  if($self->adaptor) {
    $self->adaptor($self->adaptor->db->get_ChromosomeAdaptor);
  }

  return $self;
}


#by name actually mean seq_region_name
sub name {
  my $self = shift;
  return $self->seq_region_name(@_);
}



=head2 stats

  Arg [1]    : 
  Example    : $obj->stats()
  Description: returns hashref of additional statistics stored in Chromosome table
  Returntype : hashref
  Exceptions : none
  Caller     : mapview

=cut

sub stats {
   my $self = shift;
   return %{$self->{'stats'} || {}};
}

=head2 stat

  Arg [1]    : Name of attribute to set
  Arg [2]    : Optional value to set to...
  Example    : $obj->stat( $value)
  Description: get/set for a statistic...
  Returntype : number
  Exceptions : none
  Caller     : mapview

=cut

sub stat {
  my( $self, $key, $value ) = @_;
         $self->{'stats'}{$key} = $value if defined $value;
  return $self->{'stats'}{$key};
}


## Deprecated calls - these should now use "stat"
sub xref_genes    { return $_[0]->stat('xref_genes'); }
sub known_genes   { return $_[0]->stat('known_genes'); }
sub unknown_genes { return $_[0]->stat('unknown_genes'); }
sub snps          { return $_[0]->stat('snps'); }

1;

