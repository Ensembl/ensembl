package XrefParser::RFAMParser_E;
#RFAM Parser which loads the registry from the staging servers for EnsEMBL databases. 

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
        '-host'    => 'ens-staging1',
        '-user'    => 'ensro',
      },
      {
        '-host'     => 'ens-staging2',
        '-user'     => 'ensro',
      },
  );

  return $registry;
}


1;
