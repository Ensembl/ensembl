#
# Ensembl module for Bio::EnsEMBL::MarkerFeature
#
# Cared for by Arne Stabenau<stabenau@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::MarkerFeature - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=head1 AUTHOR - Arne Stabenau

This modules is part of the Ensembl project http://www.ensembl.org

Email stabenau@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::MarkerFeature;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::EnsEMBL::Root

use Bio::EnsEMBL::SimpleFeature;

@ISA = qw(Bio::EnsEMBL::SimpleFeature);

# new() is inherieted from SeqFeature


sub marker_name {
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'marker_name'} = $value;
    }
    return $self->{'marker_name'};

}


1;
