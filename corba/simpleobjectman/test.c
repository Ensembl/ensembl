#include "simpleobjectmanager.h"
#include <stdio.h>
#include <stdlib.h>

int remove_temp_data(gpointer data)
{
  fprintf(stderr,"killing %d\n",*(int*)data);
  return 0;
}

int main(int argc,char ** argv)
{
  int i;
  int * data;
  SimpleObjectManager * som;
  SimpleObjectManagerAdaptor soma;

  som = new_SimpleObjectManager(stderr,0,0,15,"local",NULL);
  soma = SimpleObjectManager_get_Adaptor(som);

  for(i=0;i<10;i++) {
    data = malloc(sizeof(int));
    *data = i;

    fprintf(stderr,"adding %d\n",i);

    SimpleObjectManagerAdaptor_activate(&soma,data,remove_temp_data);
  }


  for(i=0;i<40;i++) {
    fprintf(stderr,"Before sleep %d\n",i);
    system("sleep 10");
    SimpleObjectManager_verify(som);

  }
}



