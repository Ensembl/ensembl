
#ifndef ARTEMIS_MYSQL_IMPL_HEADER
#define ARTEMIS_MYSQL_IMPL_HEADER

#include <mysql/mysql.h>
#include "artemis.h"
#include "simpleobjectmanager.h"

#include "artemis-exon-impl.h"
#define  ARTEMIS_MAX_SEQLENGTH 100000

Ensembl_artemis_Entry new_Ensembl_artemis_Entry(PortableServer_POA poa, MYSQL * c,const char * c_id,SimpleObjectManagerAdaptor soma,CORBA_Environment * ev);

#define RETHROW(ev, val ) if ((ev)->_major != CORBA_NO_EXCEPTION ) return val;
#define RETHROW_VOID(ev) RETHROW(ev,)

#endif
