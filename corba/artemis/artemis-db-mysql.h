
#ifndef ARTEMIS_DB_HEADER
#define ARTEMIS_DB_HEADER

#include "artemis-mysql-impl.h"

Ensembl_artemis_DB new_EA_Database(PortableServer_POA poa,MYSQL * connection,int verbose,CORBA_Environment * ev);

#endif
