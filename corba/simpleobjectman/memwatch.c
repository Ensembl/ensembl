#include "glib.h"
#include <stdio.h>


struct mw_block_str {
  char * file;
  int line;
  int number;
  void * alloc;
  int isfree;
};

static struct mw_block_str mw_block[4096];
int next_block = 0;

void show_alloc_blocks(FILE * ofp)
{
  int i;

  for(i=0;i<next_block;i++) {
    if( mw_block[i].isfree == 0 ) {
      fprintf(stderr,"%4d [%20s:%5d] %d\n",i,mw_block[i].file,mw_block[i].line,mw_block[i].number);
    }
  }
}

	      
void * mem_watch_new(int size,int number,char * file,int line)
{
  mw_block[next_block].file = file;
  mw_block[next_block].line = line;
  mw_block[next_block].number = number*size;
  mw_block[next_block].alloc = (void*)g_new0(char,number*size);
  mw_block[next_block].isfree = 0;
  next_block++;
  
  return mw_block[next_block-1].alloc;
}

void mem_watch_free(void * alloc)
{
  int i;

  if( alloc == NULL ) {
    fprintf(stderr,"Passed a NULL pointer to free");
  }

  for(i=0;i<next_block;i++) {
    if( mw_block[i].alloc == alloc ) {
      mw_block[i].isfree = 1;
      fprintf(stderr,"Freeing block %d\n",i);
      break;
    }
  }
  if( i >= next_block ) {
    fprintf(stderr,"No block allocated for %d\n",(int)alloc);
  }

  g_free(alloc);

  return;
}
      

