# A factory for ContextMapper objects
# 
# Author: Arne Stabenau
#
# Copyright EMBL/EBI 2001
#
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL:ContextMapper::Factory

=head1 SYNOPSIS

Everything in here is static. 

=head1 DESCRIPTION

The ContextMapper::Factory produces context mapper objects on
request. It contains the EnsEMBL knowledge on which modules can
translate coordinates and when and how can you get them.

Eg. It knows that from RawContig to Chromosome context and back you
use the GoldenPathMapper. It may even cache this one, as its
frequently used.

=head1 CONTACT

  Arne Stabenau: stabenau@ebi.ac.uk
  Ewan Birney  : birney@ebi.ac.uk


=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::ContextMapper::Factory;


=head2 get_Mapper

 Title   : get_Mapper
 Usage   : $factory->get_Mapper( undef, "RawContig", undef, "Chromosome" )
 Function: Get a ContextMapper. 
 Example : 
 Returns : you get undef if you cant do the map
  Args    : 1. name of source context
               this is mostly dbID for give type of object
            2. Source context type.
               This is a fixed set of types.
            3.4. same for destination.
            each of them can be the context object instead.
=cut

sub get_Mapper {
  my ($class, $sourceName, $sourceType, 
      $destinationName, $destinationType ) = @_;
  # some argument wizardry
  my $contextMapper = undef;
  return $contextMapper;
}


1;
