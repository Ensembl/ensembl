
#ifndef ARTEMIS_MYSQL_IMPL_HEADER
#define ARTEMIS_MYSQL_IMPL_HEADER

#include <mysql/mysql.h>
#include "artemis.h"

#include "artemis-exon-impl.h"


Ensembl_artemis_Entry new_Ensembl_artemis_Entry(PortableServer_POA poa, MYSQL * c,char * c_id,CORBA_Environment * ev);

#endif
