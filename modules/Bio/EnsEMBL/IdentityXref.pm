
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

SimilarityXref.pl - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...
package Bio::EnsEMBL::IdentityXref;
use vars qw(@ISA $AUTOLOAD);
use strict;


@ISA = qw( Bio::EnsEMBL::DBEntry );

=head2 query_identity

 Title   : query_identity
 Usage   : $obj->query_identity($newval)
 Function: 
 Returns : value of query_identity
 Args    : newvalue (optional)


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

 Title   : target_identity
 Usage   : $obj->target_identity($newval)
 Function: 
 Returns : value of target_identity
 Args    : newvalue (optional)


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
