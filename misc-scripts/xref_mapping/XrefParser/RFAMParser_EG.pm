package XrefParser::RFAMParser_EG;
#RFAM Parser which loads the registry from the staging servers for EnsemblGenomes databases. 

use strict;
use warnings;
use Carp;
use DBI;

use base qw( XrefParser::BaseParser);
use base qw( XrefParser::RFAMParser);
use Bio::EnsEMBL::Registry;

sub load_registry {

  my ($self, $registry) = @_;

  $registry->load_registry_from_multiple_dbs( 
      {
        '-host'     => 'mysql-eg-staging-1.ebi.ac.uk',
	'-port'     => 4160,
        '-user'     => 'ensro',
      },
      {
        '-host'     => 'mysql-eg-staging-2.ebi.ac.uk',
	'-port'     => 4275,
        '-user'     => 'ensro',
      },
 
  );
  return $registry;
}


1;
