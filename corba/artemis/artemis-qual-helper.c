
#include "artemis-qual-helper.h"



Ensembl_artemis_Qualifier new_EA_Qualifier(char * key,char * value)
{
  Ensembl_artemis_Qualifier q;

  q.name = CORBA_string_dup(key);
  q.values._buffer = (CORBA_char **) calloc(1,sizeof(CORBA_char *));
  q.values._length =1 ;
  q.values._buffer[0] = CORBA_string_dup(value);
  
  return q;
}
