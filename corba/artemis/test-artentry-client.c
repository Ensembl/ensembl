#include "artemis.h"
#include <stdio.h>

int main (int argc,char ** argv)
{
  CORBA_Environment ev;
  CORBA_ORB orb;
  char * ior;
  int len,i;
  Ensembl_artemis_Entry entry;
  Ensembl_artemis_Sequence seq;
  Ensembl_artemis_FeatureList * ftl;
  char * seqchar;
  char * seqid;
  char filebuffer[1024];
  char * end;
  FILE * ifp;

  
  CORBA_exception_init(&ev);
  orb = CORBA_ORB_init(&argc, argv, "orbit-local-orb", &ev);
  
        
  ifp = fopen("entry.ior","r");
  if( ifp == NULL ) {
    g_error("No seq.ior file! you must have an seq.ior");
    exit(-1);
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
  
  entry = CORBA_ORB_string_to_object(orb, ior, &ev);
  
  ftl = Ensembl_artemis_Entry_getAllFeatures(entry,&ev);
  fprintf(stdout,"Got length of %d\n",ftl->_length);

  for(i=0;i<ftl->_length;i++) {
    fprintf(stdout,"Attempting %d\n",i);
    fprintf(stdout,"Got feature at %s with key %s\n",
	    Ensembl_artemis_Feature_getLocation(ftl->_buffer[i],&ev),
	    Ensembl_artemis_Feature_getKey(ftl->_buffer[i],&ev)
	    );
  }

  fprintf(stdout,"Out...\n",i);


  seq = Ensembl_artemis_Entry_getSequence(entry,&ev);
  fprintf(stderr,"Got seq\n");
  len = Ensembl_artemis_Sequence_length(seq,&ev);
  fprintf(stderr,"Got length\n");
  seqchar = Ensembl_artemis_Sequence_getSubSequence(seq,1,len,&ev);

  fprintf(stdout,">%s Transcript:%d\n",Ensembl_artemis_Entry_getName(entry,&ev),
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


