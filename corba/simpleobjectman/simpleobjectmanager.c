#include "simpleobjectmanager.h"

typedef struct {
  void *_private;
  void * vepv;
  PortableServer_POA poa;
} POA_Dummy;

struct SimpleObjectManager_struct {
  FILE * log_file; /* could be NULL */
  int    use_syslog; 
  TOLManager * tolm;
  int log_flags;
  int object_lifetime;
  char * context;
  int max_objects;
  int block_size;
  GHashTable * uid_hash;
  CORBA_Environment * ev;
  int anonymous_id;
  int allow_cache;
};

/*
 * This makes a new SimpleObjectManager
 *   FILE * log_file - log file to send things to
 *   use_syslog      - use system log not the log file
 *                   (log_file has to be non null *or* use_syslog == 1)
 *   log_flags       - default flags for logging messages. We reuse the Glib error types.
 *   object_lifetime - default life time for objects
 *   context         - optional additional context for the error messages. Allows different somas
 *                     in one orb
 *   max_objects     - maximum number of objects in one orb
 *   blocksize       - number of objects to remove when max is hit
 *   ev              - The activated CORBA_Environment of the ORB
 *
 */

SimpleObjectManager * new_SimpleObjectManager(FILE * log_file,int use_syslog,int log_flags,int object_lifetime,char * context,int max_objects,int blocksize,int allow_cache,CORBA_Environment * ev)
{
  SimpleObjectManager * out;

  if( use_syslog == 1 ) {
    g_error("Have not implemented sys logging yet. Apologies!");
  }

  fprintf(stderr,"About to allocated som\n");
  out = g_new(SimpleObjectManager,1);

  fprintf(stderr,"Allocated...%d\n",(int)out);

  out->log_file = log_file;
  fprintf(stderr,"Assigning out...\n");
  out->use_syslog = use_syslog;
  out->log_flags = log_flags;
  out->object_lifetime = object_lifetime;
  if( context != NULL ) {
    out->context = g_strdup(context);
  }
  fprintf(stderr,"About to allocate tol m\n");

  out->tolm = new_TOLManager();
  out->ev = ev;
  out->max_objects = max_objects;
  out->block_size = blocksize;
  out->uid_hash   = g_hash_table_new(g_str_hash,g_str_equal);
  out->anonymous_id = 1;
  out->allow_cache  = allow_cache;

  fprintf(stderr,"About to return\n");
  return out;
}


static char * merge_class_uid(const char * class,const char * uid) 
{
  char * out;
  g_assert(uid);

  if ( class == NULL ){
    class = "NoClass";
  }

  out = g_new(char,strlen(class) + strlen(uid) + 2);
  out[0] = '\0';
  strcat(out,class);
  strcat(out,":");
  strcat(out,uid);
  
  return out;
}

static char * get_anonymous_uid(SimpleObjectManager * som)
{
  char buffer[256];
  /* FIXME: get lock on obj id now */
  sprintf(buffer,"anonymous-som-id:%d",som->anonymous_id++);
  /* release lock */

  return g_strdup(buffer);
}


/*
 * Removes the soma
 */
void                  delete_SimpleObjectManager(SimpleObjectManager * som)
{
  g_assert(som);
  if(som->context != NULL ) {
    g_free(som->context);
  }

  g_free(som);
}


/*
 * Provides an adaptor from the SimpleObjectManager. Once you have one adaptor,
 * you can assign by value many new cases
 */

SimpleObjectManagerAdaptor SimpleObjectManager_get_Adaptor(SimpleObjectManager * som)
{
  SimpleObjectManagerAdaptor out;

  out.som = som;
  out.log_flags = 0;
  out.ev = som->ev;

  return out; /* duh */
}

/*
 * log message. We reuse glib types.
 */

int SimpleObjectManagerAdaptor_log_message(SimpleObjectManagerAdaptor * soma,int message_type,char * fmt,...)
{
  char buffer[1024]; /* FIXME - Can we do this without a fixed length buffer?*/
  va_list args;

  g_assert(soma);
  g_assert(soma->som);
  g_assert(soma->som->log_file);
  g_assert(fmt);
  va_start (args, fmt);
  vsprintf(buffer,fmt,args);	
  if( soma->som->context != NULL ) {
    fprintf(soma->som->log_file,"SOMA[%s]: %s\n",soma->som->context,buffer);
  } else {
    fprintf(soma->som->log_file,"SOMA:%s\n",buffer);
  }

  va_end (args);
}

static int SimpleObjectManagerAdaptor_remove(gpointer data)
{
  return SimpleObjectManagerAdaptor_dismiss((SimpleObjectManagerAdaptor*)data,1);
}


int SimpleObjectManagerAdaptor_dismiss(SimpleObjectManagerAdaptor * soma,int should_free_object)
{
  char * uid;
  uid = soma->unique_id;

  if( uid != NULL ) {
    g_hash_table_remove(soma->som->uid_hash,uid);
  }

  if( should_free_object == TRUE ) 
    if ( soma->unique_id == NULL ) {
      SimpleObjectManagerAdaptor_log_message(soma,G_LOG_LEVEL_ERROR,"Problem - attempting to remove an object without an object id. LEAKING MEMORY");
    } else {
      TOL_remove_unique(soma->som->tolm,soma->unique_id);
      if( (*soma->remove_data)(soma->data) != 0) {
	SimpleObjectManagerAdaptor_log_message(soma,G_LOG_LEVEL_ERROR,"Unable to call remove data correctly!");
      }
    }

  /* weird but here it goes: The object itself has free'd the soma memory
   * (or should have at least). The only thing allocated by the soma was the unique
   * id. We can free this here. Revolting eh!
   */

  if( uid != NULL ) {
    g_free(uid);
  }
	
  return 0;
}


CORBA_Object SimpleObjectManagerAdaptor_reactivate(SimpleObjectManagerAdaptor * soma)
{
  CORBA_Object retval;
  POA_Dummy * d;

  g_assert(soma);
  g_assert(soma->som);
  g_assert(soma->som->tolm);
  /* remove from tolm then add it back */
  if( soma->unique_id == NULL ) {
    SimpleObjectManagerAdaptor_log_message(soma,G_LOG_LEVEL_ERROR,"Unable to reactivate object with no unique id");
    return NULL;
  }

  TOL_remove_unique(soma->som->tolm,soma->unique_id);
  TOL_add(soma->som->tolm,soma,SimpleObjectManagerAdaptor_remove,soma->som->object_lifetime,soma->unique_id);

  d = (POA_Dummy*)soma->data;
  retval = PortableServer_POA_servant_to_reference(d->poa, d, soma->ev);

  if( soma->ev->_major != CORBA_NO_EXCEPTION ) {
    fprintf(stderr,"Bugger. Got an exception!\b");
  }

  return retval;
}
  
  

/*
 * Provides the activation of the object with the adaptor
 *
 */

int SimpleObjectManagerAdaptor_activate(SimpleObjectManagerAdaptor * soma,const char * class,const char * uid,CORBA_Object obj,gpointer data,int (*remove_object)(gpointer data))
{

  g_assert(soma);
  g_assert(soma->som);
  g_assert(soma->som->tolm);
  if( TOL_list_length(soma->som->tolm) >= soma->som->max_objects ) {
    SimpleObjectManagerAdaptor_log_message(soma,G_LOG_LEVEL_MESSAGE,"Hit max object level - reaping some objects");
    TOL_remove_firstin(soma->som->tolm,soma->som->block_size);
  }
  if( uid != NULL && soma->som->allow_cache == 1) {
    soma->unique_id = merge_class_uid(class,uid);
    if( g_hash_table_lookup(soma->som->uid_hash,soma->unique_id) != NULL ) {
      SimpleObjectManagerAdaptor_log_message(soma,G_LOG_LEVEL_MESSAGE,"Identical object %s but already in cache, and not reactivated. Have to assign anonymous id to object",soma->unique_id);
      g_free(soma->unique_id);
      soma->unique_id = get_anonymous_uid(soma->som);
    }       
  } else {
    soma->unique_id = get_anonymous_uid(soma->som);
  }
  soma->data        = data;
  soma->remove_data = remove_object;
  soma->corba_object = obj;
  g_hash_table_insert(soma->som->uid_hash,soma->unique_id,soma);
  TOL_add(soma->som->tolm,data,remove_object,soma->som->object_lifetime,soma->unique_id);

  return 0;
  
}

/*
 * Runs the object manager 
 *
 */
int SimpleObjectManager_verify(SimpleObjectManager * som)
{
  g_assert(som);
  g_assert(som->tolm);

  TOL_verify(som->tolm);
}

/*
 * Finds a unique object
 */
SimpleObjectManagerAdaptor * SimpleObjectManagerAdaptor_find_object(SimpleObjectManagerAdaptor * soma,const char * class,const char * unique_id)
{
  gpointer out;
  char * key;

  g_assert(soma);
  g_assert(soma->som);
  g_assert(soma->som->uid_hash);

  if( soma->som->allow_cache == 0 ) {
    return NULL;
  }

  key = merge_class_uid(class,unique_id);
  out = g_hash_table_lookup(soma->som->uid_hash,key);

  g_free(key);
  return (SimpleObjectManagerAdaptor *) out;
}
  
  




