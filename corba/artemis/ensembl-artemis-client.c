#include "artemis.h"
#include <stdio.h>
#include <popt.h>


char * dbior = "db.ior";
int getall   = 0;
int enstart = -1;
int enend = -1;

static const
struct poptOption options[] = {
  {"dbior", 'i', POPT_ARG_STRING, &dbior, 0, "File is the database ior", "filename"},
  {"all", 'a', POPT_ARG_NONE, &getall, 0, "Get all entries", NULL},
  {"start",'s',POPT_ARG_INT,&enstart,0,"Start point in all entries","entry-start-point"},
  {"end",'e',POPT_ARG_INT,&enend,0,"End point in all entries","entry-end-point"},
  POPT_AUTOHELP
  {NULL, '\0', 0, NULL, 0, NULL, NULL}
};


int main (int argc,char ** argv)
{
  CORBA_Environment ev;
  CORBA_ORB orb;
  char * ior;
  int len,i,j;
  char * seqchar;
  char * seqid;
  char filebuffer[1024];
  char * end;
  FILE * ifp;
  char * arg;

  Ensembl_artemis_DB db;
  Ensembl_artemis_Entry entry;
  Ensembl_artemis_Sequence seq;
  Ensembl_artemis_FeatureList * ftl;
  Ensembl_artemis_QualifierList * ql;

  Ensembl_artemis_EntryNameList * enlist;

  char ** list;
  
  poptContext pcon;
  int rc;
  
  pcon=poptGetContext("ensembl-artemis-server", argc, argv, options, 0);
  
  if(argc < 1) {
    poptPrintUsage(pcon, stdout, 0);
    return 0;
  }
  
  if((rc=poptGetNextOpt(pcon)) < -1) {
    g_print("ensembl-artemis-server: bad argument %s: %s\n", 
	    poptBadOption(pcon, POPT_BADOPTION_NOALIAS),
	    poptStrerror(rc));
    exit(0);
  }

  CORBA_exception_init(&ev);
  orb = CORBA_ORB_init(&argc, argv, "orbit-local-orb", &ev);

  if( ev._major != CORBA_NO_EXCEPTION ) {
    fprintf(stderr,"Unable to bring up orb. exiting");
    exit(1);
  }

  ifp = fopen(dbior,"r");
  if( ifp == NULL ) {
    g_error("No %s file! you must have a file with dbior...",dbior);
    exit(1);
  }
  
  fgets(filebuffer,1024,ifp);
  ior = g_strdup(filebuffer);
  
  /** strip off any trailing newlines etc **/
  end = ior + strlen(ior)-1;
  for(;!isalnum(*end);end--) 
    ;
  *(end+1)='\0';
  
  /*
   * Actually get the object. So easy!
   */
  
  db = CORBA_ORB_string_to_object(orb, ior, &ev);
  
  if( ev._major != CORBA_NO_EXCEPTION ) {
    fprintf(stderr,"Unable to connect to db object, despite getting an ORB. exiting");
    exit(1);
  }

  /*
   * Narrow to db object
   */

  /*
   * process along argv
   */

  if( getall == FALSE ) {
    list = poptGetArgs(pcon);
  } else {
    enlist = Ensembl_artemis_DB_getallEntryNames(db,&ev);
    list  = g_new(char*,enlist->_length);
    j = 0;
    for(i=0;i<enlist->_length;i++) {
      if( enstart != -1 && enend != -1 ) {
	if( enstart <= i && i < enend ) {
	  list[j++] = g_strdup(enlist->_buffer[i]);
	} 
      } else {
	list[j++] = g_strdup(enlist->_buffer[i]);
      }
    }
    list[j] = NULL;
  }

  while( (arg=(*list++)) !=NULL) {
    entry = Ensembl_artemis_DB_getEntry(db,arg,&ev);
    /* check exception */
    if( ev._major != CORBA_NO_EXCEPTION ) {
      fprintf(stderr,"No entry provided for %s. Exception \n",arg);
      if( ev._major == CORBA_SYSTEM_EXCEPTION) {
	fprintf(stderr,"System exception %s\n",CORBA_exception_id(&ev));
      }
      continue;
    }
    fprintf(stderr,"About to call getAllFeatures....\n");
    ftl = Ensembl_artemis_Entry_getAllFeatures(entry,&ev);
    fprintf(stderr,"Got length of %d\n",ftl->_length);
    
      for(i=0;i<ftl->_length;i++) {
	fprintf(stdout,"Attempting %d\n",i);
	fprintf(stdout,"Got feature at %s with key %s\n",
		Ensembl_artemis_Feature_getLocation(ftl->_buffer[i],&ev),
		Ensembl_artemis_Feature_getKey(ftl->_buffer[i],&ev)
		);
	ql = Ensembl_artemis_Feature_getQualifiers(ftl->_buffer[i],&ev);
	for(j=0;j<ql->_length;j++) {
	  fprintf(stdout,"   %s %s\n",ql->_buffer[j].name,ql->_buffer[j].values._buffer[0]);
	}
	
      }
      
      fprintf(stderr,"Out...\n",i);
      
      seq = Ensembl_artemis_Entry_getSequence(entry,&ev);
      fprintf(stderr,"Got seq\n");
      len = Ensembl_artemis_Sequence_length(seq,&ev);
      fprintf(stderr,"Got length %d\n",len);
      seqchar = Ensembl_artemis_Sequence_getSubSequence(seq,1,1000,&ev);
      
      fprintf(stderr,">%s Transcript:%d\n",Ensembl_artemis_Entry_getName(entry,&ev),
	      Ensembl_artemis_Entry_getFeatureCount(entry,&ev));
      
      /* print 60 char long lines */
      len = strlen(seqchar);
      for(i=0;i<len;i++) {
	if( i%60 == 0 && i != 0 ) {
	  fputc('\n',stdout);
	} 
	fputc(toupper(seqchar[i]),stdout);
      }
  }

}




