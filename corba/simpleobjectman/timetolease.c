#include "timetolease-private.h"
#include <stdio.h>


/*
 * Builds a new TOL manager.
 *
 * All we have to do is allocate the structure and the list and pass it back
 *
 */
TOLManager * new_TOLManager(void)
{
  TOLManager * out;

  out = g_new(TOLManager,1);
  out->list= NULL;
  out->list = g_list_alloc();
  
  return out;
}

/*
 * Removes the allocated memory for the TOLManager.
 *
 * Does not (at the moment) call release for all the
 * things kept in the list
 *
 */ 
void         delete_TOLManager(TOLManager * tolm)
{
  g_assert(tolm);
  g_assert(tolm->list);

  g_list_free(tolm->list);
  g_free(tolm);
}

/*
 * Wrapper function, to both verify the list and add a new data
 *
 * Could be the only one you use... ;)
 */
int TOL_add(TOLManager * tolm,gpointer data,int (*remove_data)(gpointer),int life_in_seconds,char * unique_id)
{
  g_assert(tolm);
  g_assert(data);
  g_assert(remove_data);
  TOL_verify(tolm);
  return TOL_add_only(tolm,data,remove_data,life_in_seconds,unique_id);
}

/*
 * Just adds a new point to the list 
 */

int TOL_add_only(TOLManager * tolm,gpointer data,int (*remove_data)(gpointer),int life_in_seconds,char * unique_id)
{
  TOLWrapper * wrap;

  g_assert(tolm);
  g_assert(tolm->list);

  wrap = g_new(TOLWrapper,1);
  wrap->data = data;
  wrap->remove_data = remove_data;
  wrap->created = time(NULL);
  wrap->life_time = life_in_seconds;
  if( unique_id != NULL ) {
    wrap->unique_id = g_strdup(unique_id);
  }

  g_list_append(tolm->list,(gpointer)wrap);
  return 1;
}

int TOL_verify(TOLManager * tolm)
{
  GList * u;
  GList * temp;
  TOLWrapper * wrap;
  double diff;

  g_assert(tolm);
  g_assert(tolm->list);

  for(u = g_list_next(tolm->list);u;) {
    wrap = (TOLWrapper*)(u->data);
    diff =  difftime(time(NULL),wrap->created);
    if( diff > wrap->life_time ) {
      temp = u;
      u = g_list_next(u);
      if( (*wrap->remove_data)(wrap->data) != 0 ) {
	g_error("Removing data error in TimeToLease");
	return 1;
      }
      /* move the pointer on before we do anything else */
      if( wrap->unique_id != NULL ) {
	g_free(wrap->unique_id);
      }

      /* remove the actual pointer */
      g_free(temp->data);
      temp = g_list_remove(tolm->list,temp->data);
      /* free it */
      /*g_free(temp);*/
      
    } else {
      u = g_list_next(u);
    }
  }

  return 0;
}

void free_TOLWrapper(TOLWrapper * wrap)
{
  g_assert(wrap);
  if( wrap->unique_id != NULL ) {
    g_free(wrap->unique_id);
  }
  if( (*wrap->remove_data)(wrap->data) != 0 ) {
    g_error("Removing data error in TimeToLease");
    return;
  }
  
  /* remove the actual pointer */
  g_free(wrap);

}
  

int TOL_remove_unique(TOLManager * tolm,char * unique_id)
{
  GList * u;
  GList * temp;
  TOLWrapper * wrap;

  for(u = g_list_next(tolm->list);u;) {
    wrap = (TOLWrapper*)(u->data);
    temp = u;
    u = g_list_next(temp);
    if( wrap->unique_id != NULL && strcmp(wrap->unique_id,unique_id) == 0 ) {
      free_TOLWrapper(wrap);
      temp = g_list_remove(tolm->list,temp->data);
      return 1;
    }
  }

  return 0;
}


int TOL_remove_firstin(TOLManager * tolm,int number)
{
  GList * u;
  GList * temp;
  TOLWrapper * wrap;
  int i;

  g_assert(tolm);
  g_assert(tolm->list);
  for(i=0,u = g_list_next(tolm->list);u && i < number;i++) {

    wrap = (TOLWrapper*)(u->data);
    /* move the pointer on before we do anything else */
    temp = u;
    u = g_list_next(u);

    free_TOLWrapper(wrap);
    temp = g_list_remove(tolm->list,temp->data);
  }

  return i;
}

int TOL_list_length(TOLManager * tolm)
{
  g_assert(tolm);

  return g_list_length(tolm->list);
}











