#include "artemis-exon-impl.h"
#include <stdio.h>

/*** App-specific servant structures ***/

typedef struct
{
  POA_Ensembl_artemis_Feature servant;
  PortableServer_POA poa;
  int start;
  int end;
  int strand;
  char * created;
  char * modified;
  char * id;
  int start_phase;
}
impl_POA_Ensembl_artemis_Feature;

/*** Implementation stub prototypes ***/

static void
impl_Ensembl_artemis_Feature__destroy(impl_POA_Ensembl_artemis_Feature *
				      servant, CORBA_Environment * ev);
static CORBA_char
   *impl_Ensembl_artemis_Feature_getKey(impl_POA_Ensembl_artemis_Feature *
					servant, CORBA_Environment * ev);

static CORBA_char
   *impl_Ensembl_artemis_Feature_getLocation(impl_POA_Ensembl_artemis_Feature
					     * servant,

					     CORBA_Environment * ev);

static Ensembl_artemis_QualifierList
   *impl_Ensembl_artemis_Feature_getQualifiers
   (impl_POA_Ensembl_artemis_Feature * servant, CORBA_Environment * ev);

/*** epv structures ***/

static PortableServer_ServantBase__epv impl_Ensembl_artemis_Feature_base_epv = {
   NULL,			/* _private data */
   NULL,			/* finalize routine */
   NULL,			/* default_POA routine */
};
static POA_Ensembl_artemis_Feature__epv impl_Ensembl_artemis_Feature_epv = {
   NULL,			/* _private */
   (gpointer) & impl_Ensembl_artemis_Feature_getKey,

   (gpointer) & impl_Ensembl_artemis_Feature_getLocation,

   (gpointer) & impl_Ensembl_artemis_Feature_getQualifiers,

};

/*** vepv structures ***/

static POA_Ensembl_artemis_Feature__vepv impl_Ensembl_artemis_Feature_vepv = {
   &impl_Ensembl_artemis_Feature_base_epv,
   &impl_Ensembl_artemis_Feature_epv,
};

/*** Stub implementations ***/

static Ensembl_artemis_Feature
impl_Ensembl_artemis_Feature__create(PortableServer_POA poa,
				     CORBA_Environment * ev)
{
   Ensembl_artemis_Feature retval;
   impl_POA_Ensembl_artemis_Feature *newservant;
   PortableServer_ObjectId *objid;

   newservant = g_new0(impl_POA_Ensembl_artemis_Feature, 1);
   newservant->servant.vepv = &impl_Ensembl_artemis_Feature_vepv;
   newservant->poa = poa;
   POA_Ensembl_artemis_Feature__init((PortableServer_Servant) newservant, ev);
   objid = PortableServer_POA_activate_object(poa, newservant, ev);
   CORBA_free(objid);
   retval = PortableServer_POA_servant_to_reference(poa, newservant, ev);

   return retval;
}

Ensembl_artemis_Feature new_EA_Exon_Feature(PortableServer_POA poa,const char * id,const char * created,const char * modified,int start,int end,int strand,int phase,CORBA_Environment *ev)
{
   Ensembl_artemis_Feature retval;
   impl_POA_Ensembl_artemis_Feature *newservant;
   PortableServer_ObjectId *objid;

   newservant = g_new0(impl_POA_Ensembl_artemis_Feature, 1);
   newservant->servant.vepv = &impl_Ensembl_artemis_Feature_vepv;
   newservant->poa = poa;
   newservant->id  = g_strdup(id);
   newservant->created = g_strdup(created);
   newservant->modified = g_strdup(modified);
   newservant->start = start;
   newservant->end = end;
   newservant->start_phase = phase;
   newservant->strand = strand;

   POA_Ensembl_artemis_Feature__init((PortableServer_Servant) newservant, ev);
   objid = PortableServer_POA_activate_object(poa, newservant, ev);
   CORBA_free(objid);
   retval = PortableServer_POA_servant_to_reference(poa, newservant, ev);

   return retval;
}
  
static void
impl_Ensembl_artemis_Feature__destroy(impl_POA_Ensembl_artemis_Feature *
				      servant, CORBA_Environment * ev)
{
   PortableServer_ObjectId *objid;

   objid = PortableServer_POA_servant_to_id(servant->poa, servant, ev);
   PortableServer_POA_deactivate_object(servant->poa, objid, ev);
   CORBA_free(objid);

   POA_Ensembl_artemis_Feature__fini((PortableServer_Servant) servant, ev);
   g_free(servant);
}

static CORBA_char *
impl_Ensembl_artemis_Feature_getKey(impl_POA_Ensembl_artemis_Feature *
				    servant, CORBA_Environment * ev)
{
   CORBA_char *retval;
   fprintf(stderr,"Getting exon key...\n");
   retval = CORBA_string_dup("exon");
   return retval;
}

static CORBA_char *
impl_Ensembl_artemis_Feature_getLocation(impl_POA_Ensembl_artemis_Feature *
					 servant, CORBA_Environment * ev)
{
   CORBA_char *retval;
   char buffer[1024];
   fprintf(stderr,"Getting exon location...\n");
   if( servant->strand == 1 ) {
     sprintf(buffer,"%d..%d",servant->start,servant->end);
   } else {
     sprintf(buffer,"complement(%d..%d)",servant->start,servant->end);
   }
   retval = CORBA_string_dup(buffer);
   return retval;
}

static Ensembl_artemis_QualifierList *
impl_Ensembl_artemis_Feature_getQualifiers(impl_POA_Ensembl_artemis_Feature *
					   servant, CORBA_Environment * ev)
{
  Ensembl_artemis_QualifierList *retval;
   
  retval = CORBA_sequence_Ensembl_artemis_Qualifier__alloc();
  retval->_buffer = (Ensembl_artemis_Qualifier *) calloc (3,sizeof(Ensembl_artemis_Qualifier));
  
  retval->_buffer[0] = new_EA_Qualifier("created",servant->created);
  retval->_buffer[1] = new_EA_Qualifier("modified",servant->modified);
  retval->_buffer[2] = new_EA_Qualifier("exon_id",servant->id);
  retval->_length = 3;
  retval->_maximum = 3;
  /*retval->_buffer[3] = new_EA_Qualifier("exon_id",servant->id);*/
  
  return retval;
}



