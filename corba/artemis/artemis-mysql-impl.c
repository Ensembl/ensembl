#include "artemis-mysql-impl.h"
#include "timetolease.h" 
#include <stdio.h>
#include <stdlib.h>


/*
 * This file makes implementations of Artemis Entries and sequence objects
 * using a mysql database system. It is where most of the clever stuff occurs
 * complex SQL queries for example to build up transcript and feature type
 * objects.
 */

/*** App-specific servant structures ***/
typedef struct {
  POA_Ensembl_artemis_Sequence servant;
  PortableServer_POA poa;
  MYSQL * connection; /* database connection */
  char * contig_id;   /* need this to retrieve the data */
  long length;        /* we store this in memory as it is cheap */
  SimpleObjectManagerAdaptor soma;
} impl_POA_Ensembl_artemis_Sequence;

/*
 * This is an implementation only of Transcripts, not of anything else
 * Other features are implementated elsewhere
 */
 
typedef struct {
  POA_Ensembl_artemis_Feature servant;
  PortableServer_POA poa;
  MYSQL * connection;   /* database connection */
  char * transcript_id; /* transcript id for this feature... */
  char * location;
  SimpleObjectManagerAdaptor soma;
} impl_POA_Ensembl_artemis_Feature;

typedef struct {
  POA_Ensembl_artemis_Entry servant;
  PortableServer_POA poa;
  MYSQL * connection;
  char * contig_id;
  long transcript_number;

  /*this is so we can cache transcript/exons*/
  char ** transcript_ids;
  char ** exon_ids;
  int has_done_transcript;

  SimpleObjectManagerAdaptor soma;
} impl_POA_Ensembl_artemis_Entry;

/*** Implementation stub prototypes ***/
static void
impl_Ensembl_artemis_Sequence__destroy(impl_POA_Ensembl_artemis_Sequence *
				       servant, CORBA_Environment * ev);
static CORBA_char
   *impl_Ensembl_artemis_Sequence_getSubSequence
   (impl_POA_Ensembl_artemis_Sequence * servant, CORBA_long start,
    CORBA_long end, CORBA_Environment * ev);

static CORBA_long
impl_Ensembl_artemis_Sequence_length(impl_POA_Ensembl_artemis_Sequence *
				     servant, CORBA_Environment * ev);

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

static void impl_Ensembl_artemis_Entry__destroy(impl_POA_Ensembl_artemis_Entry
						* servant,

						CORBA_Environment * ev);
static CORBA_char
   *impl_Ensembl_artemis_Entry_getName(impl_POA_Ensembl_artemis_Entry *
				       servant, CORBA_Environment * ev);

static CORBA_long
impl_Ensembl_artemis_Entry_getFeatureCount(impl_POA_Ensembl_artemis_Entry *
					   servant, CORBA_Environment * ev);

static Ensembl_artemis_FeatureList
   *impl_Ensembl_artemis_Entry_getAllFeatures(impl_POA_Ensembl_artemis_Entry *
					      servant,

					      CORBA_Environment * ev);

static Ensembl_artemis_Sequence
impl_Ensembl_artemis_Entry_getSequence(impl_POA_Ensembl_artemis_Entry *
				       servant, CORBA_Environment * ev);

/*** epv structures ***/
static PortableServer_ServantBase__epv impl_Ensembl_artemis_Sequence_base_epv
   = {
   NULL,			/* _private data */
   NULL,			/* finalize routine */
   NULL,			/* default_POA routine */
};
static POA_Ensembl_artemis_Sequence__epv impl_Ensembl_artemis_Sequence_epv = {
   NULL,			/* _private */
   (gpointer) & impl_Ensembl_artemis_Sequence_getSubSequence,

   (gpointer) & impl_Ensembl_artemis_Sequence_length,

};

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

static PortableServer_ServantBase__epv impl_Ensembl_artemis_Entry_base_epv = {
   NULL,			/* _private data */
   NULL,			/* finalize routine */
   NULL,			/* default_POA routine */
};
static POA_Ensembl_artemis_Entry__epv impl_Ensembl_artemis_Entry_epv = {
   NULL,			/* _private */
   (gpointer) & impl_Ensembl_artemis_Entry_getName,

   (gpointer) & impl_Ensembl_artemis_Entry_getFeatureCount,

   (gpointer) & impl_Ensembl_artemis_Entry_getAllFeatures,

   (gpointer) & impl_Ensembl_artemis_Entry_getSequence,

};

/*** vepv structures ***/
static POA_Ensembl_artemis_Sequence__vepv impl_Ensembl_artemis_Sequence_vepv = {
   &impl_Ensembl_artemis_Sequence_base_epv,
   &impl_Ensembl_artemis_Sequence_epv,
};

static POA_Ensembl_artemis_Feature__vepv impl_Ensembl_artemis_Feature_vepv = {
   &impl_Ensembl_artemis_Feature_base_epv,
   &impl_Ensembl_artemis_Feature_epv,
};

static POA_Ensembl_artemis_Entry__vepv impl_Ensembl_artemis_Entry_vepv = {
   &impl_Ensembl_artemis_Entry_base_epv,
   &impl_Ensembl_artemis_Entry_epv,
};

/*** Stub implementations ***/
static Ensembl_artemis_Sequence
impl_Ensembl_artemis_Sequence__create(PortableServer_POA poa,
				      CORBA_Environment * ev)
{
   Ensembl_artemis_Sequence retval;
   impl_POA_Ensembl_artemis_Sequence *newservant;
   PortableServer_ObjectId *objid;

   newservant = g_new0(impl_POA_Ensembl_artemis_Sequence, 1);
   newservant->servant.vepv = &impl_Ensembl_artemis_Sequence_vepv;
   newservant->poa = poa;
   POA_Ensembl_artemis_Sequence__init((PortableServer_Servant) newservant,
				      ev);
   objid = PortableServer_POA_activate_object(poa, newservant, ev);
   CORBA_free(objid);
   retval = PortableServer_POA_servant_to_reference(poa, newservant, ev);

   return retval;
}
  

impl_POA_Ensembl_artemis_Sequence * new_impl_EA_Sequence(PortableServer_POA poa,MYSQL * c,char * c_id,SimpleObjectManagerAdaptor soma)
{
   impl_POA_Ensembl_artemis_Sequence *newservant;

   newservant = g_new0(impl_POA_Ensembl_artemis_Sequence, 1);
   newservant->servant.vepv = &impl_Ensembl_artemis_Sequence_vepv;
   newservant->poa = poa;
   newservant->connection =c ;
   newservant->soma = soma;
   newservant->contig_id = c_id;

   return newservant;
}

int remove_Sequence_func(gpointer data)
{
  impl_POA_Ensembl_artemis_Sequence * servant;
  servant = (impl_POA_Ensembl_artemis_Sequence*) data;
  SimpleObjectManagerAdaptor_log_message(&servant->soma,G_LOG_LEVEL_MESSAGE,"Removing Sequence %s",servant->contig_id);
  impl_Ensembl_artemis_Sequence__destroy(servant,servant->soma.ev);
  return 0;
}

Ensembl_artemis_Sequence new_Ensembl_artemis_Sequence(PortableServer_POA poa,MYSQL *c,char * c_id,SimpleObjectManagerAdaptor soma,CORBA_Environment * ev)
{
  Ensembl_artemis_Sequence retval;
  impl_POA_Ensembl_artemis_Sequence *newservant;
  PortableServer_ObjectId *objid;
  
  /* MySQL stuff */
  char sqlbuffer[1024];
  MYSQL_RES * result;
  MYSQL_ROW row;
  int state;
  int no;
  SimpleObjectManagerAdaptor * ret;
  PortableServer_Servant tempservant;

  if( c == NULL ) {
    SimpleObjectManagerAdaptor_log_message(&soma,G_LOG_LEVEL_ERROR,"Passed in NULL connection. Cannot build sequence object on null connection");
    CORBA_exception_set(ev, CORBA_USER_EXCEPTION, NULL,NULL);
    return;
  }

  if( c_id == NULL ) {
    SimpleObjectManagerAdaptor_log_message(&soma,G_LOG_LEVEL_ERROR,"Passed in NULL c_id. Cannot build sequence object on null contig name");
    CORBA_exception_set_system(ev,ex_CORBA_UNKNOWN,CORBA_COMPLETED_NO);
    return;
  }

  /* see if we have this object already... */
  ret = SimpleObjectManagerAdaptor_find_object(&soma,"EASequence",c_id);
  if( ret != NULL ) {
    SimpleObjectManagerAdaptor_log_message(ret,G_LOG_LEVEL_MESSAGE,"Reactivating sequence object");
    retval = SimpleObjectManagerAdaptor_reactivate(ret);
    return retval;
  } else {

    /*
     * Check in exists in database *now* not later!
     */
    sprintf(sqlbuffer,"SELECT length from contig where id = '%s'",c_id);
    state = mysql_query(c,sqlbuffer);
    result = mysql_store_result(c);
    
    no = mysql_num_rows(result);
    if( no == 0 ) {
      SimpleObjectManagerAdaptor_log_message(&soma,G_LOG_LEVEL_ERROR,"No sequence of this name %s",c_id);
      CORBA_exception_set_system(ev,ex_CORBA_UNKNOWN,CORBA_COMPLETED_NO);
      return;
    }
        
    row = mysql_fetch_row(result);
    newservant = new_impl_EA_Sequence(poa,c,g_strdup(c_id),soma);
    newservant->length = atol(row[0]);
    
    SimpleObjectManagerAdaptor_log_message(&soma,G_LOG_LEVEL_DEBUG,"Made sequence %s with length",c_id,newservant->length);
    mysql_free_result(result);
    
    
    /*
     * ok. We can rock and roll now 
     */
    POA_Ensembl_artemis_Sequence__init((PortableServer_Servant) newservant, ev);
    objid = PortableServer_POA_activate_object(poa, newservant, ev);
    CORBA_free(objid);
    retval = PortableServer_POA_servant_to_reference(poa, newservant, ev);
    SimpleObjectManagerAdaptor_activate(&newservant->soma,"EASequence",c_id,retval,(gpointer)newservant,remove_Sequence_func);
    return retval;
  }
}


static void
impl_Ensembl_artemis_Sequence__destroy(impl_POA_Ensembl_artemis_Sequence *
				       servant, CORBA_Environment * ev)
{
   PortableServer_ObjectId *objid;

   objid = PortableServer_POA_servant_to_id(servant->poa, servant, ev);
   PortableServer_POA_deactivate_object(servant->poa, objid, ev);
   CORBA_free(objid);
   g_free(servant->contig_id);
   POA_Ensembl_artemis_Sequence__fini((PortableServer_Servant) servant, ev);
   g_free(servant);
}

static CORBA_char *
impl_Ensembl_artemis_Sequence_getSubSequence(impl_POA_Ensembl_artemis_Sequence
					     * servant, CORBA_long start,
					     CORBA_long end,
					     CORBA_Environment * ev)
{
   CORBA_char *retval;
   char sqlbuffer[1024];
   char * temp;
   MYSQL_RES * result;
   MYSQL_ROW row;
   int state;
   char * str;

   sprintf(sqlbuffer,"SELECT sequence from dna where contig = '%s'",servant->contig_id);
   state = mysql_query(servant->connection,sqlbuffer);
   result = mysql_store_result(servant->connection);
   /* FIXME: error trapping here */
   
   row = mysql_fetch_row(result);

   str = row[0];
   /** check for DNA sequences in here? **/
 
   temp = calloc(end-start+2,sizeof(char));
   strncpy(temp,str+start-1,end-start+1);
   temp[end-start+1] = '\0';

   retval = CORBA_string_dup(temp);
   free(temp);
   mysql_free_result(result);
   return retval;
}

static CORBA_long
impl_Ensembl_artemis_Sequence_length(impl_POA_Ensembl_artemis_Sequence *
				     servant, CORBA_Environment * ev)
{
   CORBA_long retval;
   retval = servant->length;
   return retval;
}

int remove_Feature_func(gpointer data)
{
  impl_POA_Ensembl_artemis_Feature * servant;
  servant = (impl_POA_Ensembl_artemis_Feature*) data;
  SimpleObjectManagerAdaptor_log_message(&servant->soma,0,"Removing Feature");
  impl_Ensembl_artemis_Feature__destroy(servant,servant->soma.ev);
  return 0;
}

impl_POA_Ensembl_artemis_Feature * new_EA_Feature(PortableServer_POA poa,char * location,char * transcript_id,SimpleObjectManagerAdaptor soma)
{
   impl_POA_Ensembl_artemis_Feature *newservant;

   newservant = g_new0(impl_POA_Ensembl_artemis_Feature, 1);
   newservant->servant.vepv = &impl_Ensembl_artemis_Feature_vepv;
   newservant->poa = poa;
   newservant->location = location;
   newservant->transcript_id = transcript_id;
   newservant->soma = soma;

   return newservant;
}

Ensembl_artemis_Feature new_Ensembl_artemis_Feature(PortableServer_POA poa,char * location,char * transcript_id,SimpleObjectManagerAdaptor soma,CORBA_Environment * ev) 
{
   Ensembl_artemis_Feature retval;
   impl_POA_Ensembl_artemis_Feature *newservant;
   PortableServer_ObjectId *objid;

   newservant = new_EA_Feature(poa,location,transcript_id,soma);
   POA_Ensembl_artemis_Feature__init((PortableServer_Servant) newservant, ev);
   objid = PortableServer_POA_activate_object(poa, newservant, ev);
   CORBA_free(objid);
   retval = PortableServer_POA_servant_to_reference(poa, newservant, ev);
   SimpleObjectManagerAdaptor_activate(&newservant->soma,"EAFeatureT",transcript_id,retval,(gpointer)newservant,remove_Feature_func);
   return retval;
}


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

static void
impl_Ensembl_artemis_Feature__destroy(impl_POA_Ensembl_artemis_Feature *
				      servant, CORBA_Environment * ev)
{
   PortableServer_ObjectId *objid;

   objid = PortableServer_POA_servant_to_id(servant->poa, servant, ev);
   PortableServer_POA_deactivate_object(servant->poa, objid, ev);
   CORBA_free(objid);
   fprintf(stderr,"Removing trans id\n");
   if( servant->transcript_id != NULL ) 
     g_free(servant->transcript_id);
   fprintf(stderr,"Removing loc id\n");
   if( servant->location != NULL )
     g_free(servant->location);
   POA_Ensembl_artemis_Feature__fini((PortableServer_Servant) servant, ev);
   fprintf(stderr,"Removing servant\n");
   g_free(servant);
}

static CORBA_char *
impl_Ensembl_artemis_Feature_getKey(impl_POA_Ensembl_artemis_Feature *
				    servant, CORBA_Environment * ev)
{
   CORBA_char *retval;
   retval = CORBA_string_dup("CDS");
   return retval;
}

static CORBA_char *
impl_Ensembl_artemis_Feature_getLocation(impl_POA_Ensembl_artemis_Feature *
					 servant, CORBA_Environment * ev)
{
   CORBA_char *retval;
   retval = CORBA_string_dup(servant->location);
   return retval;
}

static Ensembl_artemis_QualifierList *
impl_Ensembl_artemis_Feature_getQualifiers(impl_POA_Ensembl_artemis_Feature *
					   servant, CORBA_Environment * ev)
{
   Ensembl_artemis_QualifierList *retval;

   retval = CORBA_sequence_Ensembl_artemis_Qualifier__alloc();
   retval->_buffer = CORBA_sequence_Ensembl_artemis_Qualifier_allocbuf (1);
   retval->_length = 1;
   retval->_maximum = 1;

   retval->_buffer[0] = new_EA_Qualifier("transcript_id",servant->transcript_id);
   CORBA_sequence_set_release(retval,1);
   return retval;
}

static Ensembl_artemis_Entry
impl_Ensembl_artemis_Entry__create(PortableServer_POA poa,
				   CORBA_Environment * ev)
{
   Ensembl_artemis_Entry retval;
   impl_POA_Ensembl_artemis_Entry *newservant;
   PortableServer_ObjectId *objid;

   newservant = g_new0(impl_POA_Ensembl_artemis_Entry, 1);
   newservant->servant.vepv = &impl_Ensembl_artemis_Entry_vepv;
   newservant->poa = poa;
   POA_Ensembl_artemis_Entry__init((PortableServer_Servant) newservant, ev);
   objid = PortableServer_POA_activate_object(poa, newservant, ev);
   CORBA_free(objid);
   retval = PortableServer_POA_servant_to_reference(poa, newservant, ev);

   return retval;
}

impl_POA_Ensembl_artemis_Entry * new_EA_Entry(PortableServer_POA poa, MYSQL * c,const char * c_id,SimpleObjectManagerAdaptor soma,long trans_number)
{
  impl_POA_Ensembl_artemis_Entry * newservant;

  newservant = g_new0(impl_POA_Ensembl_artemis_Entry, 1);
  newservant->servant.vepv = &impl_Ensembl_artemis_Entry_vepv;
  newservant->poa = poa;

  g_assert(c);
  g_assert(c_id);

  newservant->connection = c;
  newservant->contig_id  = g_strdup(c_id);
  newservant->transcript_number = trans_number;
  newservant->soma = soma;
  newservant->has_done_transcript = 0;
  newservant->transcript_ids = NULL;
  newservant->exon_ids = NULL;
  return newservant;
}

int remove_Entry_func(gpointer data)
{
  char ** ids;
  impl_POA_Ensembl_artemis_Entry * servant;
  servant = (impl_POA_Ensembl_artemis_Entry*) data;
  SimpleObjectManagerAdaptor_log_message(&servant->soma,0,"Removing Entry %s",servant->contig_id);

  fprintf(stderr,"Jumping into destroy\n");
  impl_Ensembl_artemis_Entry__destroy(servant,servant->soma.ev);
  return 0;
}

Ensembl_artemis_Entry new_Ensembl_artemis_Entry(PortableServer_POA poa, MYSQL * c,const char * c_id,SimpleObjectManagerAdaptor soma,CORBA_Environment * ev)
{
   Ensembl_artemis_Entry retval;
   impl_POA_Ensembl_artemis_Entry *newservant;
   PortableServer_ObjectId *objid;

   /* MySQL stuff */
   char sqlbuffer[1024];
   MYSQL_RES * result;
   MYSQL_ROW row;
   int state;
   int no;
   
   SimpleObjectManagerAdaptor * ret;

   if( c == NULL ) {
     SimpleObjectManagerAdaptor_log_message(&soma,G_LOG_LEVEL_ERROR,"Passed in NULL connection. Cannot build sequence object on null connection");
     CORBA_exception_set_system(ev,ex_CORBA_UNKNOWN,CORBA_COMPLETED_NO);
     return;
   }
   
   if( c_id == NULL ) {
     SimpleObjectManagerAdaptor_log_message(&soma,G_LOG_LEVEL_ERROR,"Passed in NULL c_id. Cannot build sequence object on null contig name");
     CORBA_exception_set_system(ev,ex_CORBA_UNKNOWN,CORBA_COMPLETED_NO);
     return;
   }

   ret = SimpleObjectManagerAdaptor_find_object(&soma,"EAEntry",c_id);
   if( ret != NULL ) {
     SimpleObjectManagerAdaptor_log_message(ret,G_LOG_LEVEL_MESSAGE,"Reactivating entry object");
     return SimpleObjectManagerAdaptor_reactivate(ret);
   } else {
     /* check it does exist */
     sprintf(sqlbuffer,"SELECT 1 from contig where id = '%s'",c_id);
     state = mysql_query(c,sqlbuffer);
     result = mysql_store_result(c);
     no = mysql_num_rows(result);
     mysql_free_result(result);
     
     if( no == 0 ) {
       SimpleObjectManagerAdaptor_log_message(&soma,G_LOG_LEVEL_ERROR,"No entry of this name %s",c_id);
       CORBA_exception_set_system(ev,ex_CORBA_UNKNOWN,CORBA_COMPLETED_NO);
       return;
     }
     
    newservant = new_EA_Entry(poa,c,c_id,soma,0);
    POA_Ensembl_artemis_Entry__init((PortableServer_Servant) newservant, ev);
    objid = PortableServer_POA_activate_object(poa, newservant, ev);
    CORBA_free(objid);
    retval = PortableServer_POA_servant_to_reference(poa, newservant, ev);
    SimpleObjectManagerAdaptor_activate(&newservant->soma,"EAEntry",c_id,retval,(gpointer)newservant,remove_Entry_func);
    return retval;
  }
}

static void
impl_Ensembl_artemis_Entry__destroy(impl_POA_Ensembl_artemis_Entry * servant,
				    CORBA_Environment * ev)
{
   PortableServer_ObjectId *objid;
   char ** ids;

   objid = PortableServer_POA_servant_to_id(servant->poa, servant, ev);
   PortableServer_POA_deactivate_object(servant->poa, objid, ev);
   CORBA_free(objid);

   if( servant->contig_id != NULL ) {
     g_free(servant->contig_id);
   }
   if( servant->transcript_ids != NULL ) {
     for(ids = servant->transcript_ids;*ids != NULL;ids++) {
       fprintf(stderr,"Removing %s\n",*ids);
       g_free(*ids);
     } 
   }
   g_free(servant->transcript_ids);
   servant->transcript_ids = NULL;
   if( servant->exon_ids != NULL ) {
     
     for(ids = servant->exon_ids;*ids != NULL;ids++) {
       fprintf(stderr,"Removing %s %d\n",*ids,(int)*ids);
       g_free(*ids);
     } 
   }
   g_free(servant->exon_ids);
   servant->exon_ids = NULL;

   POA_Ensembl_artemis_Entry__fini((PortableServer_Servant) servant, ev);
   fprintf(stderr,"Going to free servant\n");
   g_free(servant);
}

static CORBA_char *
impl_Ensembl_artemis_Entry_getName(impl_POA_Ensembl_artemis_Entry * servant,
				   CORBA_Environment * ev)
{
   CORBA_char *retval;
   retval = CORBA_string_dup(servant->contig_id);
   return retval;
}

static CORBA_long
impl_Ensembl_artemis_Entry_getFeatureCount(impl_POA_Ensembl_artemis_Entry *
					   servant, CORBA_Environment * ev)
{
   CORBA_long retval;
   retval = servant->transcript_number;
   return retval;
}

static Ensembl_artemis_FeatureList *
impl_Ensembl_artemis_Entry_getAllFeatures(impl_POA_Ensembl_artemis_Entry *
					  servant, CORBA_Environment * ev)
{
  Ensembl_artemis_Feature temp_buffer[1024];
  Ensembl_artemis_FeatureList *retval;

  impl_POA_Ensembl_artemis_Feature * feats;

  /* cache stuff */
  char * temp_trans[128];
  char * temp_exon[1024];
  int t,e,l;
  
  /* MySQL stuff */
  char sqlbuffer[1024];
  MYSQL_RES * result;
  MYSQL_ROW row;
  
  char loc[2048];
  char * trans_id;
  int i;
  MYSQL * c;
  int state;
  int touched;
  int no;
  int seen = 0;
  char ** id = NULL;
  
  SimpleObjectManagerAdaptor * ret;
  c = servant->connection;
  
  /*S_ASSERT(servant->soma,c);*/

  /* check to see if we have made this before */
  fprintf(stderr,"Got into get all entries... %d\n",servant->has_done_transcript);
  i =0;
  t = 0; e= 0;
  if( servant->has_done_transcript == 1 ){
    SimpleObjectManagerAdaptor_log_message(&servant->soma,G_LOG_LEVEL_DEBUG,"Attempting cache");

    for(id = servant->transcript_ids;*id != NULL;id++) {
      if( (ret=SimpleObjectManagerAdaptor_find_object(&servant->soma,"EAFeatureT",*id)) != NULL ) {
	SimpleObjectManagerAdaptor_log_message(&servant->soma,G_LOG_LEVEL_DEBUG,"Cacheing %s",*id);
	temp_buffer[i++] = SimpleObjectManagerAdaptor_reactivate(ret);
	} else { 
	  SimpleObjectManagerAdaptor_log_message(&servant->soma,G_LOG_LEVEL_DEBUG,"Failed cache on %s",*id);
	  break;
	}
    }
    if( *id == NULL ) {
      /* see if we can do the exons */
      fprintf(stderr,"Doing the exon loop\n");

      for(id = servant->exon_ids;*id != NULL;id++) {
	if( (ret=SimpleObjectManagerAdaptor_find_object(&servant->soma,"EAFeatureE",*id)) != NULL ) {
	  temp_buffer[i++] = SimpleObjectManagerAdaptor_reactivate(ret);
	  SimpleObjectManagerAdaptor_log_message(&servant->soma,G_LOG_LEVEL_DEBUG,"Cacheing %s",*id);
	} else {
	  SimpleObjectManagerAdaptor_log_message(&servant->soma,G_LOG_LEVEL_DEBUG,"Failed cache on %s",*id);
	  break;
	}
      }/* each id */
    } /* over exons */
  }
  
  /* id == NULL means we have not got into caching. 
   * *id == NULL means we were successful
   */
  if( id != NULL && *id == NULL ) {
    SimpleObjectManagerAdaptor_log_message(&servant->soma,G_LOG_LEVEL_MESSAGE,"Reactivated feature objects");
  } else {
    /** we have to make the SQL query **/
    i =0; /* reset buffer */
    sprintf(sqlbuffer,"SELECT p1.transcript,p1.rank,p2.start,p2.end,p2.strand,p2.id,p2.created,p2.modified from exon_transcript as p1,exon as p2 where p2.contig = '%s' and p1.exon = p2.id order by p1.transcript,p1.rank",servant->contig_id);
    SimpleObjectManagerAdaptor_log_message(&servant->soma,G_LOG_LEVEL_DEBUG,"Going to issue [%s] for select statement over transcripts",sqlbuffer);
    
   
    state = mysql_query(c,sqlbuffer);
    if( state == -1 ) {
      retval = CORBA_sequence_Ensembl_artemis_Feature__alloc();
      retval->_buffer = NULL;
      retval->_length = 0;
      retval->_maximum = 0;
      return retval;
    }
   
    result = mysql_store_result(c);
    
    trans_id = NULL;
    seen = 0;
    i = 0;
    
   
    if( mysql_num_rows(result) == 0 ) {
      /* return empty list */
      retval = CORBA_sequence_Ensembl_artemis_Feature__alloc();
      retval->_buffer = NULL;
      retval->_length = 0;
      retval->_maximum = 0;
      
      return retval;
    }
    
    while ( row = mysql_fetch_row(result) ) {
      
      /* transcript id */
      if( trans_id == NULL || strcmp(trans_id,row[0]) != 0 ) {
	/* put away old feature object */
	if( seen != 0 ) {
	  strcat(loc,")");
	  temp_buffer[i++] = new_Ensembl_artemis_Feature(servant->poa,g_strdup(loc),trans_id,servant->soma,ev);
	  temp_trans[t++] = g_strdup(trans_id);
	} 
	
	/* get ready for the next transcript */
	trans_id = g_strdup(row[0]);
	loc[0] = '\0';
	strcat(loc,"join(");
	touched = 0;
	seen =1;
      }
      
      SimpleObjectManagerAdaptor_log_message(&servant->soma,G_LOG_LEVEL_DEBUG,"Making a new exon");
      
      temp_buffer[i++] = new_EA_Exon_Feature(servant->poa,row[5],row[6],row[7],atol(row[2]),atol(row[3]),strcmp(row[4],"1") == 0 ? 1 : -1,0,servant->soma,ev);
      temp_exon[e++] = g_strdup(row[5]);
      fprintf(stderr,"Allocated %s to temp %d\n",temp_exon[e-1],(int)temp_exon[e-1]);
      if( touched == 1 ) {
	strcat(loc,",");
      } else {
	touched = 1;
      }
      
      
      if( strcmp(row[4],"1") == 0) {
	strcat(loc,row[2]);
	strcat(loc,"..");
	strcat(loc,row[3]);
      } else {
	strcat(loc,"complement(");
	strcat(loc,row[2]);
	strcat(loc,"..");
	strcat(loc,row[3]);
	strcat(loc,")");
      }
      
    }
    mysql_free_result(result);

    /* last feature */
    if( seen != 0 ) {
      strcat(loc,")");
      temp_buffer[i++] = new_Ensembl_artemis_Feature(servant->poa,g_strdup(loc),trans_id,servant->soma,ev);
      temp_trans[t++] = g_strdup(trans_id);
    }
    
    servant->transcript_ids = g_new(char*,t+2);
    for(l=0;l<t;l++)
      servant->transcript_ids[l] = temp_trans[l];
    servant->transcript_ids[l] = NULL;

    servant->exon_ids = g_new(char*,e+2);
    for(l=0;l<e;l++)
      servant->exon_ids[l] = temp_exon[l];
    servant->exon_ids[l] = NULL;
    servant->has_done_transcript = 1;
  }
  
  /** return buffer. All in temp buffer **/
    
  
  no = i;
  retval = CORBA_sequence_Ensembl_artemis_Feature__alloc();
  if( no != 0 ) {
    retval->_buffer = (Ensembl_artemis_Feature *) CORBA_sequence_Ensembl_artemis_Feature_allocbuf(no);
  }
  
  retval->_maximum = no;
  retval->_length = no;
  for(i=0;i<no;i++) {
    retval->_buffer[i] = temp_buffer[i];
  }
  
  
  CORBA_sequence_set_release(retval,1);
  return retval;
}

static Ensembl_artemis_Sequence
impl_Ensembl_artemis_Entry_getSequence(impl_POA_Ensembl_artemis_Entry *
				       servant, CORBA_Environment * ev)
{
   Ensembl_artemis_Sequence retval;

   char sqlbuffer[1024];
   MYSQL_RES * result;
   MYSQL_ROW row;
   int state;
   gpointer data;

   sprintf(sqlbuffer,"SELECT length from contig where id = '%s'",servant->contig_id);

   state = mysql_query(servant->connection,sqlbuffer);
   result = mysql_store_result(servant->connection);
   mysql_free_result(result);
   retval = new_Ensembl_artemis_Sequence(servant->poa,servant->connection,servant->contig_id,servant->soma,ev);
   return retval;
}












