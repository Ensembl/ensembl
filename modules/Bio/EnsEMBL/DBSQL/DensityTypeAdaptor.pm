#
# EnsEMBL module for Bio::EnsEMBL::DBSQL::DensityTypeAdaptor
#
#
# Copyright EMBL/EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code




#mysql> desc density_type;
#+-----------------+---------------------+------+-----+---------+----------------+
#| Field           | Type                | Null | Key | Default | Extra          |
#+-----------------+---------------------+------+-----+---------+----------------+
#| density_type_id | int(11)             |      | PRI | NULL    | auto_increment |
#| analysis_id     | int(11)             |      | MUL | 0       |                |
#| block_size      | int(11)             |      |     | 0       |                |
#| value_type      | enum('sum','ratio') |      |     | sum     |                |
#+-----------------+---------------------+------+-----+---------+----------------+

=head1 NAME

Bio::EnsEMBL::DBSQL::DensityTypeAdaptor

=head1 SYNOPSIS

my $density_type_adaptor = $database_adaptor->get_DensityTypeAdaptor();
@density_types = @{$density_type_adaptor->fetch_all_by_Slice($slice)};

=head1 DESCRIPTION

Simple Feature Adaptor - database access for simple features

=head1 AUTHOR - Ewan Birney

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::DensityTypeAdaptor;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
#use Bio::EnsEMBL::DensityType;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor);


=head2 store

  Arg [1]    : list of Bio::EnsEMBL::DensityType @dt
               the density types to store in the database
  Example    : $density_type->store(1234, @density_types);
  Description: Stores a list of density type objects in the database
  Returntype : none
  Exceptions : thrown if @dt is not defined
               or if any elements of @dt are not Bio::EnsEMBL::DensityType 
  Caller     : general

=cut
#mysql> desc density_type;
#+-----------------+---------------------+------+-----+---------+----------------+
#| Field           | Type                | Null | Key | Default | Extra          |
#+-----------------+---------------------+------+-----+---------+----------------+
#| density_type_id | int(11)             |      | PRI | NULL    | auto_increment |
#| analysis_id     | int(11)             |      | MUL | 0       |                |
#| block_size      | int(11)             |      |     | 0       |                |
#| value_type      | enum('sum','ratio') |      |     | sum     |                |
#+-----------------+---------------------+------+-----+---------+----------------+

sub store{
  my ($self,@dt) = @_;

  if( scalar(@dt) == 0 ) {
    throw("Must call store with list of Density Types");
  }

  my $sth = $self->prepare
    ("INSERT INTO density_type (analysis_id,".
                                  "block_size, value_type) ". 
    "VALUES (?,?,?)");

  my $db = $self->db();
  my $analysis_adaptor = $db->get_AnalysisAdaptor();

 FEATURE: foreach my $dt ( @dt ) {
    
    if( !ref $dt || !$dt->isa("Bio::EnsEMBL::DensityType") ) {
      throw("Density Type must be an Ensembl DensityType, " .
            "not a [".ref($dt)."]");
    }
    
    if($dt->is_stored($db)) {
      warning("DensityType [".$dt->density_type_id."] is already stored" .
              " in this database.");
      next FEATURE;
    }
    
    if(!defined($dt->analysis)) {
      throw("An analysis must be attached to the density type to be stored.");
    }
    
    #store the analysis if it has not been stored yet
    if(!$dt->analysis->is_stored($db)) {
      throw("Analysis must already exist???? But does not??");
    }
	

    $sth->execute($dt->analysis->dbID,
		  $dt->{'block_size'}, $dt->{'value_type'});


    # next two lines are to set the density type as stored
    $dt->dbID($sth->{'mysql_insertid'});
    $dt->adaptor($self);
  }
}

    
    
   

