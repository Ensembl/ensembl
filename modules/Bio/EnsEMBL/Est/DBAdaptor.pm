
=Head1 NAME - Bio::EnsEMBL::Est::DBAdaptor

=head1 SYNOPSIS

    $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -user   => 'root',
        -dbname => 'pog',
        -host   => 'caldy',
        -driver => 'mysql',
        );

    
=head1 DESCRIPTION

This adaptor is specialized in
the way, that it just returns a DnaAlignFeatureAdaptor that can be connected
to a different DB than the core Database, making it possible to keep the Ests in a 
different db. 

The extra est - dna_align_feature database only needs to contain the dna_align_feature
table as in the core db. all other data is used from core DB, whose DBAdaptor has to be
attached to the DBAdaptor in this directory.

This might not be the final solution to the problem ....


=head1 CONTACT

stabenau@ebi.ac.uk

=head1 APPENDIX
=cut


# Let the code begin...

package Bio::EnsEMBL::Est::DBAdaptor;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DBSQL::DBConnection;

@ISA = qw(Bio::EnsEMBL::DBSQL::DBConnection);

# following the only reason why this DBAdaptor exists
sub get_DnaAlignFeatureAdaptor {
  my $self = shift;
  
  $self->_get_adaptor( "Bio::EnsEMBL::DBSQL::DnaAlignFeatureAdaptor" );
}

# following adaptors make this db work like a core DB for hopefully all
# relevant cases ...
# by this feature->db->get_XXXAdaptor calls still work and give the right adaptors

sub get_RawContigAdaptor {
  my $self = shift;

  $self->core_DBAdaptor()->get_RawContigAdaptor();
}

sub get_AnalysisAdaptor {
  my $self = shift;

  $self->core_DBAdaptor()->get_AnalysisAdaptor();
}

sub get_AssemblyMapperAdaptor {
  my $self = shift;

  $self->core_DBAdaptor()->get_AssemblyMapperAdaptor();
}

sub core_DBAdaptor {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    $self->{'_core_DBAdaptor'} = $arg;
  }

  return $self->{'_core_DBAdaptor'};
}

1;
