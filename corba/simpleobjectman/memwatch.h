
#ifndef MEM_WATCH_HEADER
#define MEM_WATCH_HEADER

#include <stdio.h>

#ifdef  g_new
#undef  g_new
#endif

#define g_new(class,number) (class *)mem_watch_new(sizeof(class),number,__FILE__,__LINE__)

#ifdef  g_new0
#undef  g_new0
#endif

#define g_new0(class,number) (class *)mem_watch_new(sizeof(class),number,__FILE__,__LINE__)

#ifdef g_free
#undef g_free
#endif

#define g_free(pos) mem_watch_free(pos);
void * mem_watch_new(int bytes,int number,char * file,int line);
void show_alloc_blocks(FILE * ofp);
#endif
