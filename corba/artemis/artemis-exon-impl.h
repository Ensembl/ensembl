
#ifndef ARTEMIS_EXON_IMPL_HEADER
#define ARTEMIS_EXON_IMPL_HEADER

#include "artemis.h"
#include "artemis-qual-helper.h"
#include "simpleobjectmanager.h"

Ensembl_artemis_Feature new_EA_Exon_Feature(PortableServer_POA poa,const char * id,const char * created,const char * modified,int start,int end,int strand,int phase,SimpleObjectManagerAdaptor soma,CORBA_Environment *ev);

#endif

