
#
# EnsEMBL module for GoXref.pl
#
# Cared for by Arne Stabenau <stabenau@ebi.ac.uk>
#
# Copyright EnsEMBL 2000-2003
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::GoXref

=head1 SYNOPSIS

my $goxref = Bio::EnsEMBL::GoXref->new;
$goxref->linkage_type()
  => "computational"
  => "experimental"
  => "curated"

=head1 CONTACT

Post questions to the ensembl development list: <ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...
package Bio::EnsEMBL::GoXref;
use vars qw(@ISA $AUTOLOAD);
use strict;


@ISA = qw( Bio::EnsEMBL::DBEntry );



=head2 linkage_type

  Arg [1]    : (optional) string $value
    allowed are "computational", "experimental", "curated"
  Example    : none
  Description: Getter/Setter for linkage_type
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub linkage_type {
   my $obj = shift;
   if( @_ ) {
       my $value = shift;
       $obj->{'linkage_type'} = $value;
   }
   return $obj->{'linkage_type'};

}

1;
