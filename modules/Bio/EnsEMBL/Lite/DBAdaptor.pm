
=Head1 NAME - Bio::EnsEMBL::Lite::DBAdaptor

=head1 SYNOPSIS

    $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -user   => 'root',
        -dbname => 'pog',
        -host   => 'caldy',
        -driver => 'mysql',
        );


    $db->get_GeneAdaptor

=head1 DESCRIPTION

This DBAdaptor provides the link to the lite database. It allows the rapid creation of drawable objects, mainly of gene objects. The created Gene objects are not entirely complete, but they are connected to the DBSQL::GeneAdaptor to allow lazy loading of missing data. 

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Lite::DBAdaptor;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DBSQL::DBConnection;

@ISA = qw(Bio::EnsEMBL::DBSQL::DBConnection);


#Override constructor inherited by Bio::EnsEMBL::DBSQL::DBConnection
sub new {
  my($class, @args) = @_;

  #call superclass constructor
  my $self = $class->SUPER::new(@args);
  
  return $self;
}

sub get_GeneAdaptor {
  my $self = shift;

    return 
      $self->_get_adaptor("Bio::EnsEMBL::Lite::GeneAdaptor");
}

sub get_LandmarkMarkerAdaptor {
  my $self = shift;

    return 
      $self->_get_adaptor("Bio::EnsEMBL::Lite::LandmarkMarkerAdaptor");
}

sub get_RepeatFeatureAdaptor {
  my $self = shift;

  return $self->_get_adaptor("Bio::EnsEMBL::Lite::RepeatFeatureAdaptor");
}

sub get_DensityAdaptor {
  my $self = shift;

  return $self->_get_adaptor("Bio::EnsEMBL::Lite::DensityPlot::DensityAdaptor");
}

sub get_ChromosomeAdaptor {
  my $self = shift;

  return $self->_get_adaptor("Bio::EnsEMBL::Lite::ChromosomeAdaptor");
}


sub get_SNPAdaptor {
  my $self = shift;

  return $self->_get_adaptor("Bio::EnsEMBL::Lite::SNPAdaptor");
}


1;
