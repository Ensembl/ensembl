
#ifndef ENSEMBL_SEQ_SERVER
#define ENSEMBL_SEQ_SERVER


#include "bioseq.h"
#include <mysql/mysql.h>
#include <sys/time.h>


BioSource_SeqDB new_EnsEMBL_BioSource_SeqDB(PortableServer_POA poa, 
					    MYSQL * connection,
					    CORBA_Environment * ev);

BioSource_Seq   new_EnsEMBL_BioSource_Seq(PortableServer_POA poa,
					  MYSQL * connection,
					  char * dna_database_id,
					  CORBA_Environment * ev);

#endif


