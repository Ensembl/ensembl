#
# Ensembl module for Bio::EnsEMBL::Lite::ChromosomeAdaptor
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Lite::ChromosomeAdaptor - DESCRIPTION of Object

=head1 SYNOPSIS

$chromosome_adaptor = $db_adaptor->get_ChromosomeAdaptor();
$chromosome = $chromosome_adaptor->fetch_by_chr_name('12');

=head1 DESCRIPTION

This is a database adaptor used to retrieve chromosome objects from a database.

=head1 AUTHOR - Ewan Birney

This modules is part of the Ensembl project http://www.ensembl.org

Email birney@ebi.ac.uk


=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Lite::ChromosomeAdaptor;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::EnsEMBL::Root

use Bio::EnsEMBL::DBSQL::ChromosomeAdaptor;
use Bio::EnsEMBL::Chromosome;

@ISA = qw(Bio::EnsEMBL::DBSQL::ChromosomeAdaptor);

# new inherited from DBSQL::ChromosomeAdaptor


=head2 _create_object_from_hashref

  Args       : hash ref containing a row of the chromosome table
  Example    : $self->_create_object_from_hashref
  Description: Creates object from hash reference
  Returntype : Bio::Chromosome object
  Exceptions : none
  Caller     : general

=cut

sub _create_object_from_hashref {
  my( $self,$h ) =@_;
  local $_;
  my $chr = new Bio::EnsEMBL::Chromosome(
    -adaptor  => $self,        -dbID   => $h->{'chromosome_id'},
    -chr_name => $h->{'name'}, -length => $h->{'length'},
  );
  $chr->stat( $_, $h->{$_} ) foreach
    ( grep { $_ ne 'chromosome_id' && $_ ne 'name' && $_ ne 'length' } keys %$h);
  return $self->{'_chr_cache'     }->{$h->{'chromosome_id'}} = 
         $self->{'_chr_name_cache'}->{$h->{'name'}         } = $chr ;
}

1;
