#include "artemis.h"
#include <stdio.h>

int main (int argc,char ** argv)
{
  CORBA_Environment ev;
  CORBA_ORB orb;
  char * ior;
  int len,i;
  Ensembl_artemis_Sequence seq;
  char * seqchar;
  char * seqid;
  char filebuffer[1024];
  char * end;
  FILE * ifp;

  
  CORBA_exception_init(&ev);
  orb = CORBA_ORB_init(&argc, argv, "orbit-local-orb", &ev);
  
        
  ifp = fopen("seq.ior","r");
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
  
  seq = CORBA_ORB_string_to_object(orb, ior, &ev);
  len = Ensembl_artemis_Sequence_length(seq);

  seqchar = Ensembl_artemis_Sequence_getSubSequence(seq,1,len,&ev);

  fprintf(stdout,">SomeSequence\n",);

  /* print 60 char long lines */
  len = strlen(seqchar);
  for(i=0;i<len;i++) {
    if( i%60 == 0 && i != 0 ) {
      fputc('\n',stdout);
    } 
    fputc(toupper(seqchar[i]),stdout);
  }
  
}


