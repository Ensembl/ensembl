#
#
# BioPerl module for SimilarityXref.pl
#
# Cared for by Emmanuel Mongin <mongin@ebi.ac.uk>
#
# Copyright Emmanuel Mongin
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::IdentiyXref

=head1 SYNOPSIS

my $xref = Bio::EnsEMBL::IdentityXref->new;

=head1 CONTACT

Post questions to the ensembl development list: <ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...
package Bio::EnsEMBL::IdentityXref;
use vars qw(@ISA $AUTOLOAD);
use strict;


@ISA = qw( Bio::EnsEMBL::DBEntry );

=head2 new
  
  See Bio::EnsEMBL::DBEntry::new

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    my ($query_identity, $target_identity) = $self->_rearrange(
        [qw(QUERY_IDENTITY TARGET_IDENTITY)], @args);
    if(defined $query_identity) { $self->query_identity($query_identity); }
    if(defined $target_identity) { $self->target_identity($target_identity); }
    return $self;
}

=head2 query_identity

  Arg [1]    : (optional) string $value
  Example    : $query_identity = $id_xref->query_identity;
  Description: Getter/Setter for query identity
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub query_identity{
   my $obj = shift;
   if( @_ ) {
       my $value = shift;
       $obj->{'query_identity'} = $value;
   }
   return $obj->{'query_identity'};

}


=head2 target_identity

  Arg [1]    : (optional) string $value
  Example    : $target_identity = $id_xref->target_identity;
  Description: Getter/Setter for query identity
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub target_identity{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'target_identity'} = $value;
    }
    return $obj->{'target_identity'};

}

1;
