#include "ensembl-seq-server.h"

#include <stdio.h>

/*** App-specific servant structures ***/

typedef struct
{
  POA_BioSource_Seq servant;
  PortableServer_POA poa;
  char * database_dna_id;
  MYSQL * connection;
}
impl_POA_BioSource_Seq;

typedef struct
{
  POA_BioSource_SeqDB servant;
  PortableServer_POA poa;
  MYSQL * connection;
}
impl_POA_BioSource_SeqDB;


/*** Implementation stub prototypes ***/

static void impl_BioSource_Seq__destroy(impl_POA_BioSource_Seq * servant,
					CORBA_Environment * ev);
static CORBA_char *impl_BioSource_Seq_seq(impl_POA_BioSource_Seq * servant,
					  CORBA_Environment * ev);

static CORBA_char *impl_BioSource_Seq_subseq(impl_POA_BioSource_Seq * servant,
					     CORBA_long start,
					     CORBA_long end,
					     CORBA_Environment * ev);

static CORBA_char *impl_BioSource_Seq_id(impl_POA_BioSource_Seq * servant,
					 CORBA_Environment * ev);

static CORBA_char *impl_BioSource_Seq_accession(impl_POA_BioSource_Seq *
						servant,

						CORBA_Environment * ev);

static CORBA_char *impl_BioSource_Seq_description(impl_POA_BioSource_Seq *
						  servant,

						  CORBA_Environment * ev);

static BioSource_seqtype
impl_BioSource_Seq_type(impl_POA_BioSource_Seq * servant,

			CORBA_Environment * ev);

static void
impl_BioSource_Seq_release(impl_POA_BioSource_Seq * servant,
			   CORBA_Environment * ev);

static void impl_BioSource_SeqDB__destroy(impl_POA_BioSource_SeqDB * servant,
					  CORBA_Environment * ev);
static BioSource_Seq
impl_BioSource_SeqDB_get_Seq_by_id(impl_POA_BioSource_SeqDB * servant,
				   CORBA_char * id, CORBA_Environment * ev);

static BioSource_Seq
impl_BioSource_SeqDB_get_Seq_by_acc(impl_POA_BioSource_SeqDB * servant,
				    CORBA_char * acc, CORBA_Environment * ev);

static void
impl_BioSource_SeqDB_release(impl_POA_BioSource_SeqDB * servant,
			     CORBA_Environment * ev);

/*** epv structures ***/

static PortableServer_ServantBase__epv impl_BioSource_Seq_base_epv = {
   NULL,			/* _private data */
   NULL,			/* finalize routine */
   NULL,			/* default_POA routine */
};
static POA_BioSource_Seq__epv impl_BioSource_Seq_epv = {
   NULL,			/* _private */
   (gpointer) & impl_BioSource_Seq_seq,

   (gpointer) & impl_BioSource_Seq_subseq,

   (gpointer) & impl_BioSource_Seq_id,

   (gpointer) & impl_BioSource_Seq_accession,

   (gpointer) & impl_BioSource_Seq_description,

   (gpointer) & impl_BioSource_Seq_type,

};
static POA_BioSource_ReleaseableObject__epv
   impl_BioSource_Seq_BioSource_ReleaseableObject_epv = {
   NULL,			/* _private */
   (gpointer) & impl_BioSource_Seq_release,
};
static PortableServer_ServantBase__epv impl_BioSource_SeqDB_base_epv = {
   NULL,			/* _private data */
   NULL,			/* finalize routine */
   NULL,			/* default_POA routine */
};
static POA_BioSource_SeqDB__epv impl_BioSource_SeqDB_epv = {
   NULL,			/* _private */
   (gpointer) & impl_BioSource_SeqDB_get_Seq_by_id,

   (gpointer) & impl_BioSource_SeqDB_get_Seq_by_acc,

};
static POA_BioSource_ReleaseableObject__epv
   impl_BioSource_SeqDB_BioSource_ReleaseableObject_epv = {
   NULL,			/* _private */
   (gpointer) & impl_BioSource_SeqDB_release,
};

/*** vepv structures ***/

static POA_BioSource_Seq__vepv impl_BioSource_Seq_vepv = {
   &impl_BioSource_Seq_base_epv,
   &impl_BioSource_Seq_BioSource_ReleaseableObject_epv,
   &impl_BioSource_Seq_epv,
};
static POA_BioSource_SeqDB__vepv impl_BioSource_SeqDB_vepv = {
   &impl_BioSource_SeqDB_base_epv,
   &impl_BioSource_SeqDB_BioSource_ReleaseableObject_epv,
   &impl_BioSource_SeqDB_epv,
};

/*** Stub implementations ***/

/*
 * new implementation structure for a bioseq. Does not activate it
 *
 */

static impl_POA_BioSource_Seq * new_impl_BioSource_Seq(PortableServer_POA poa,
						       MYSQL * c,
						       char * id)
{
   impl_POA_BioSource_Seq *newservant;


   newservant = g_new0(impl_POA_BioSource_Seq, 1);
   newservant->servant.vepv = &impl_BioSource_Seq_vepv;
   newservant->poa = poa;
   newservant->connection = c;
   newservant->database_dna_id = id;

   return newservant;
}

BioSource_Seq   new_EnsEMBL_BioSource_Seq(PortableServer_POA poa,MYSQL * connection,char * dna_database_id,CORBA_Environment * ev)
{
  BioSource_Seq retval;
  impl_POA_BioSource_Seq *newservant;
  PortableServer_ObjectId *objid;

  /* MySQL stuff */
  char sqlbuffer[1024];
  MYSQL_RES * result;
  int state;


  if( connection == NULL ) {
    /*g_warn("Passed in NULL connection. Cannot build sequence object on null connection");*/
    /* yikes- should do something!*/
  }

  if( dna_database_id == NULL ) {
    /*g_warn("Passed in NULL dna database id - horrifying!");*/
    /* do something */
  }

  /*
   * Check in exists in database *now* not later!
   */
  
  fprintf(stderr,"About to make connection!\n");
  sprintf(sqlbuffer,"SELECT contig from dna where contig = '%s'",dna_database_id);
  state = mysql_query(connection,sqlbuffer);
  result = mysql_store_result(connection);
  fprintf(stderr,"Made and stored!\n");
  if( mysql_num_rows(result) != 1 ) {
    fprintf(stderr,"No sequences!\n");
    /*g_warn("Bad news - in asking for dna_id %s - no rows %d!",dna_database_id,mysql_num_rows(result));*/
  }
  fprintf(stderr,"About to free!\n");
  mysql_free_result(result);
  fprintf(stderr,"Freed!\n");

  newservant = new_impl_BioSource_Seq(poa,connection,dna_database_id);
  
  /*
   * ok. We can rock and roll now 
   */
   POA_BioSource_Seq__init((PortableServer_Servant) newservant, ev);
   fprintf(stderr,"Freed!\n");
   objid = PortableServer_POA_activate_object(poa, newservant, ev);
   CORBA_free(objid);
   fprintf(stderr,"Freed!\n");
   retval = PortableServer_POA_servant_to_reference(poa, newservant, ev);
   fprintf(stderr,"About to return!\n");
   return retval;
}


static BioSource_Seq
impl_BioSource_Seq__create(PortableServer_POA poa, CORBA_Environment * ev)
{
   BioSource_Seq retval;
   impl_POA_BioSource_Seq *newservant;
   PortableServer_ObjectId *objid;

   newservant = g_new0(impl_POA_BioSource_Seq, 1);
   newservant->servant.vepv = &impl_BioSource_Seq_vepv;
   newservant->poa = poa;
   POA_BioSource_Seq__init((PortableServer_Servant) newservant, ev);
   objid = PortableServer_POA_activate_object(poa, newservant, ev);
   CORBA_free(objid);
   retval = PortableServer_POA_servant_to_reference(poa, newservant, ev);

   return retval;
}

static void
impl_BioSource_Seq__destroy(impl_POA_BioSource_Seq * servant,
			    CORBA_Environment * ev)
{
   PortableServer_ObjectId *objid;

   objid = PortableServer_POA_servant_to_id(servant->poa, servant, ev);
   PortableServer_POA_deactivate_object(servant->poa, objid, ev);
   CORBA_free(objid);

   POA_BioSource_Seq__fini((PortableServer_Servant) servant, ev);
   g_free(servant);
}

static CORBA_char *
impl_BioSource_Seq_seq(impl_POA_BioSource_Seq * servant,
		       CORBA_Environment * ev)
{
   CORBA_char *retval;
   char sqlbuffer[1024];
   MYSQL_RES * result;
   MYSQL_ROW row;
   int state;
   char * str;

   sprintf(sqlbuffer,"SELECT sequence from dna where contig = '%s'",servant->database_dna_id);
   state = mysql_query(servant->connection,sqlbuffer);
   result = mysql_store_result(servant->connection);
   /* FIXME: error trapping here */

   row = mysql_fetch_row(result);

   str = row[0];
   /* fprintf(stderr,"Returned %s",str); */

   /** check for DNA sequences in here? **/

   retval = CORBA_string_dup(str);
   return retval;
}

static CORBA_char *
impl_BioSource_Seq_subseq(impl_POA_BioSource_Seq * servant,
			  CORBA_long start,
			  CORBA_long end, CORBA_Environment * ev)
{
   CORBA_char *retval;

   return retval;
}

static CORBA_char *
impl_BioSource_Seq_id(impl_POA_BioSource_Seq * servant,
		      CORBA_Environment * ev)
{
   CORBA_char *retval;

   return CORBA_string_dup(servant->database_dna_id);
}

static CORBA_char *
impl_BioSource_Seq_accession(impl_POA_BioSource_Seq * servant,
			     CORBA_Environment * ev)
{
   CORBA_char *retval;

   return CORBA_string_dup(servant->database_dna_id);
}

static CORBA_char *
impl_BioSource_Seq_description(impl_POA_BioSource_Seq * servant,
			       CORBA_Environment * ev)
{
   CORBA_char *retval;

   return CORBA_string_dup("");
}

static BioSource_seqtype
impl_BioSource_Seq_type(impl_POA_BioSource_Seq * servant,
			CORBA_Environment * ev)
{
   BioSource_seqtype retval;
   return BioSource_DNA;
}

static void
impl_BioSource_Seq_release(impl_POA_BioSource_Seq * servant,
			   CORBA_Environment * ev)
{
}

static BioSource_SeqDB
impl_BioSource_SeqDB__create(PortableServer_POA poa, CORBA_Environment * ev)
{
   BioSource_SeqDB retval;
   impl_POA_BioSource_SeqDB *newservant;
   PortableServer_ObjectId *objid;

   newservant = g_new0(impl_POA_BioSource_SeqDB, 1);
   newservant->servant.vepv = &impl_BioSource_SeqDB_vepv;
   newservant->poa = poa;
   POA_BioSource_SeqDB__init((PortableServer_Servant) newservant, ev);
   objid = PortableServer_POA_activate_object(poa, newservant, ev);
   CORBA_free(objid);
   retval = PortableServer_POA_servant_to_reference(poa, newservant, ev);

   return retval;
}

static void
impl_BioSource_SeqDB__destroy(impl_POA_BioSource_SeqDB * servant,
			      CORBA_Environment * ev)
{
   PortableServer_ObjectId *objid;

   objid = PortableServer_POA_servant_to_id(servant->poa, servant, ev);
   PortableServer_POA_deactivate_object(servant->poa, objid, ev);
   CORBA_free(objid);

   POA_BioSource_SeqDB__fini((PortableServer_Servant) servant, ev);
   g_free(servant);
}

static BioSource_Seq
impl_BioSource_SeqDB_get_Seq_by_id(impl_POA_BioSource_SeqDB * servant,
				   CORBA_char * id, CORBA_Environment * ev)
{
   BioSource_Seq retval;

   return retval;
}

static BioSource_Seq
impl_BioSource_SeqDB_get_Seq_by_acc(impl_POA_BioSource_SeqDB * servant,
				    CORBA_char * acc, CORBA_Environment * ev)
{
   BioSource_Seq retval;

   return retval;
}

static void
impl_BioSource_SeqDB_release(impl_POA_BioSource_SeqDB * servant,
			     CORBA_Environment * ev)
{
}










