#include "artemis-mysql-impl.h"

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
  int verbose;        /* should we chat to stderr? */
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
} impl_POA_Ensembl_artemis_Feature;

typedef struct {
  POA_Ensembl_artemis_Entry servant;
  PortableServer_POA poa;
  MYSQL * connection;
  char * contig_id;
  long transcript_number;
  int verbose; /* should we chat to stderr? */
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

impl_POA_Ensembl_artemis_Sequence * new_impl_EA_Sequence(PortableServer_POA poa,MYSQL * c,char * c_id,int verbose)
{
   impl_POA_Ensembl_artemis_Sequence *newservant;


   newservant = g_new0(impl_POA_Ensembl_artemis_Sequence, 1);
   newservant->servant.vepv = &impl_Ensembl_artemis_Sequence_vepv;
   newservant->poa = poa;
   newservant->connection =c ;
   newservant->contig_id = c_id;
   newservant->verbose  = verbose;

   return newservant;
}

Ensembl_artemis_Sequence new_Ensembl_artemis_Sequence(PortableServer_POA poa,MYSQL *c,char * c_id,int verbose,CORBA_Environment * ev)
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
  
  
  if( c == NULL ) {
    fprintf(stderr,"Passed in NULL connection. Cannot build sequence object on null connection");
    /* yikes- should do something!*/
  }

  if( c_id == NULL ) {
    fprintf(stderr,"Passed in NULL dna database id - horrifying!");
    /* do something */
  }

  /*
   * Check in exists in database *now* not later!
   */

  fprintf(stderr,"About to make select call...\n");
  
  sprintf(sqlbuffer,"SELECT length from contig where id = '%s'",c_id);

  fprintf(stderr,"Before state\n");
  state = mysql_query(c,sqlbuffer);
  fprintf(stderr,"Before store\n");
  result = mysql_store_result(c);
  no = mysql_num_rows(result);
  if( no == 0 ) {
    fprintf(stderr,"No sequence of this name %s",c_id);
    return;
  }


  fprintf(stderr,"Before fetch\n");
  row = mysql_fetch_row(result);

  fprintf(stderr,"Making implementation\n");
  newservant = new_impl_EA_Sequence(poa,c,c_id,verbose);
  
  fprintf(stderr,"Made newservant...\n");

  newservant->length = atol(row[0]);

  fprintf(stderr,"Made and stored! %d\n",newservant->length);

  fprintf(stderr,"About to free!\n");
  mysql_free_result(result);
  fprintf(stderr,"Freed!\n");


  
  /*
   * ok. We can rock and roll now 
   */
   POA_Ensembl_artemis_Sequence__init((PortableServer_Servant) newservant, ev);
   fprintf(stderr,"Freed!\n");
   objid = PortableServer_POA_activate_object(poa, newservant, ev);
   CORBA_free(objid);
   fprintf(stderr,"Freed!\n");
   retval = PortableServer_POA_servant_to_reference(poa, newservant, ev);
   fprintf(stderr,"About to return!\n");
   return retval;
}


static void
impl_Ensembl_artemis_Sequence__destroy(impl_POA_Ensembl_artemis_Sequence *
				       servant, CORBA_Environment * ev)
{
   PortableServer_ObjectId *objid;

   objid = PortableServer_POA_servant_to_id(servant->poa, servant, ev);
   PortableServer_POA_deactivate_object(servant->poa, objid, ev);
   CORBA_free(objid);

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
   /* fprintf(stderr,"Returned %s",str); */

   /** check for DNA sequences in here? **/

   temp = calloc(end-start+2,sizeof(char));
   strncpy(temp,str+start-1,end-start+1);
   temp[end-start+1] = '\0';

   retval = CORBA_string_dup(temp);
   return retval;
}

static CORBA_long
impl_Ensembl_artemis_Sequence_length(impl_POA_Ensembl_artemis_Sequence *
				     servant, CORBA_Environment * ev)
{
   CORBA_long retval;
   fprintf(stderr,"Going to return %d\n",servant->length);
   retval = servant->length;
   return retval;
}

impl_POA_Ensembl_artemis_Feature * new_EA_Feature(PortableServer_POA poa,char * location,char * transcript_id)
{
   impl_POA_Ensembl_artemis_Feature *newservant;

   newservant = g_new0(impl_POA_Ensembl_artemis_Feature, 1);
   newservant->servant.vepv = &impl_Ensembl_artemis_Feature_vepv;
   newservant->poa = poa;
   newservant->location = location;
   newservant->transcript_id = transcript_id;

   return newservant;
}

Ensembl_artemis_Feature new_Ensembl_artemis_Feature(PortableServer_POA poa,char * location,char * transcript_id,CORBA_Environment * ev) 
{
   Ensembl_artemis_Feature retval;
   impl_POA_Ensembl_artemis_Feature *newservant;
   PortableServer_ObjectId *objid;

   newservant = new_EA_Feature(poa,location,transcript_id);
   POA_Ensembl_artemis_Feature__init((PortableServer_Servant) newservant, ev);
   objid = PortableServer_POA_activate_object(poa, newservant, ev);
   CORBA_free(objid);
   retval = PortableServer_POA_servant_to_reference(poa, newservant, ev);

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

   POA_Ensembl_artemis_Feature__fini((PortableServer_Servant) servant, ev);
   g_free(servant);
}

static CORBA_char *
impl_Ensembl_artemis_Feature_getKey(impl_POA_Ensembl_artemis_Feature *
				    servant, CORBA_Environment * ev)
{
  
   CORBA_char *retval;
   fprintf(stderr,"Getting into getKey\n");
   retval = CORBA_string_dup("CDS");
   return retval;
}

static CORBA_char *
impl_Ensembl_artemis_Feature_getLocation(impl_POA_Ensembl_artemis_Feature *
					 servant, CORBA_Environment * ev)
{
   CORBA_char *retval;
   fprintf(stderr,"Getting into getLocation\n");
   retval = CORBA_string_dup(servant->location);
   return retval;
}

static Ensembl_artemis_QualifierList *
impl_Ensembl_artemis_Feature_getQualifiers(impl_POA_Ensembl_artemis_Feature *
					   servant, CORBA_Environment * ev)
{
   Ensembl_artemis_QualifierList *retval;

   retval = CORBA_sequence_Ensembl_artemis_Qualifier__alloc();
   retval->_buffer = (Ensembl_artemis_Qualifier * ) calloc (1,sizeof(Ensembl_artemis_Qualifier));
   retval->_length = 1;
   retval->_maximum = 1;

   retval->_buffer[0] = new_EA_Qualifier("transcript_id",servant->transcript_id);

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

impl_POA_Ensembl_artemis_Entry * new_EA_Entry(PortableServer_POA poa, MYSQL * c,char * c_id,long trans_number)
{
  impl_POA_Ensembl_artemis_Entry * newservant;

  newservant = g_new0(impl_POA_Ensembl_artemis_Entry, 1);
  newservant->servant.vepv = &impl_Ensembl_artemis_Entry_vepv;
  newservant->poa = poa;

  g_assert(c);
  g_assert(c_id);

  newservant->connection = c;
  newservant->contig_id  = c_id;
  newservant->transcript_number = trans_number;

  return newservant;
}

Ensembl_artemis_Entry new_Ensembl_artemis_Entry(PortableServer_POA poa, MYSQL * c,char * c_id,CORBA_Environment * ev)
{
   Ensembl_artemis_Entry retval;
   impl_POA_Ensembl_artemis_Entry *newservant;
   PortableServer_ObjectId *objid;

   newservant = new_EA_Entry(poa,c,c_id,0);

   POA_Ensembl_artemis_Entry__init((PortableServer_Servant) newservant, ev);
   objid = PortableServer_POA_activate_object(poa, newservant, ev);
   CORBA_free(objid);
   retval = PortableServer_POA_servant_to_reference(poa, newservant, ev);

   return retval;
}
  


static void
impl_Ensembl_artemis_Entry__destroy(impl_POA_Ensembl_artemis_Entry * servant,
				    CORBA_Environment * ev)
{
   PortableServer_ObjectId *objid;

   objid = PortableServer_POA_servant_to_id(servant->poa, servant, ev);
   PortableServer_POA_deactivate_object(servant->poa, objid, ev);
   CORBA_free(objid);

   POA_Ensembl_artemis_Entry__fini((PortableServer_Servant) servant, ev);
   g_free(servant);
}

static CORBA_char *
impl_Ensembl_artemis_Entry_getName(impl_POA_Ensembl_artemis_Entry * servant,
				   CORBA_Environment * ev)
{
   CORBA_char *retval;
   fprintf(stderr,"Going to return name %s\n",servant->contig_id);
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


   c = servant->connection;
   fprintf(stderr,"Got in...\n",no);


   sprintf(sqlbuffer,"SELECT p1.transcript,p1.rank,p2.start,p2.end,p2.strand,p2.id,p2.created,p2.modified from exon_transcript as p1,exon as p2 where p2.contig = '%s' and p1.exon = p2.id order by p1.transcript,p1.rank",servant->contig_id);
   fprintf(stderr,"Going to issue %s\n",sqlbuffer);

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
   fprintf(stderr,"Going to call loop...%d\n",state);

   if( mysql_num_rows(result) == 0 ) {
     fprintf(stderr,"Got none...!\n");
     /* return empty list */
     retval = CORBA_sequence_Ensembl_artemis_Feature__alloc();
     retval->_buffer = NULL;
     retval->_length = 0;
     retval->_maximum = 0;
     
     return retval;
   }

   while ( row = mysql_fetch_row(result) ) {
     fprintf(stderr,"Got a row for transcript %s\n",trans_id == NULL ? "NoTranscript" : trans_id);

     /* transcript id */
     if( trans_id == NULL || strcmp(trans_id,row[0]) != 0 ) {
       /* put away old feature object */
       if( seen != 0 ) {
	 strcat(loc,")");
	 fprintf(stderr,"Going to be making %s %s\n",trans_id,loc);
	 temp_buffer[i++] = new_Ensembl_artemis_Feature(servant->poa,g_strdup(loc),trans_id,ev);
       } 

       /* get ready for the next transcript */
       trans_id = g_strdup(row[0]);
       loc[0] = '\0';
       strcat(loc,"join(");
       touched = 0;
       seen =1;
     }

     fprintf(stderr,"Now processing %s %s  %s %s %s %s\n",trans_id,row[4],row[2],row[3],row[6],row[7]);
     temp_buffer[i++] = new_EA_Exon_Feature(servant->poa,row[5],row[6],row[7],atol(row[2]),atol(row[3]),1,0,ev);

     if( touched == 1 ) {
       strcat(loc,",");
     } else {
       touched = 1;
     }


     if( strcmp(row[4],"1") ) {
       fprintf(stderr,"Going to use forward\n");
       strcat(loc,row[2]);
       strcat(loc,"..");
       strcat(loc,row[3]);
     } else {
       fprintf(stderr,"Going to use backward\n");
       strcat(loc,"complement(");
       strcat(loc,row[2]);
       strcat(loc,"..");
       strcat(loc,row[3]);
       strcat(loc,")");
     }
     
   }

   /* last feature */
   if( seen != 0 ) {
     strcat(loc,")");
     fprintf(stderr,"Going to be making %d %s %s\n",i,trans_id,loc);
     temp_buffer[i++] = new_Ensembl_artemis_Feature(servant->poa,g_strdup(loc),trans_id,ev);
   } 

   fprintf(stderr,"Out of loop...\n");

   no = i;
   retval = CORBA_sequence_Ensembl_artemis_Feature__alloc();
   if( no != 0 ) {
     retval->_buffer = (Ensembl_artemis_Feature *) calloc (no,sizeof(Ensembl_artemis_Feature));
   }

   retval->_maximum = no;
   retval->_length = no;
   for(i=0;i<no;i++) {
     retval->_buffer[i] = temp_buffer[i];
   }


   /*CORBA_sequence_set_release(retval,1);*/
   fprintf(stderr,"Passing in buffer with %d\n",CORBA_sequence_get_release(retval));
   return retval;
}

static Ensembl_artemis_Sequence
impl_Ensembl_artemis_Entry_getSequence(impl_POA_Ensembl_artemis_Entry *
				       servant, CORBA_Environment * ev)
{
   Ensembl_artemis_Sequence retval;

   fprintf(stderr,"Getting into give sequence\n");

   retval = new_Ensembl_artemis_Sequence(servant->poa,servant->connection,servant->contig_id,1,ev);
   fprintf(stderr,"Made new sequence!\n");

   return retval;
}












