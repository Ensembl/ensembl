
#ifndef SIMPLEOBJECTMANAGER_HEADER
#define SIMPLEOBJECTMANAGER_HEADER

/*
 * SimpleObjectManager. Written by Ewan Birney <birney@sanger.ac.uk>
 *
 * Lots of problems in this.
 *
 * Copyright Ewan Birney. GPL'd. No warranty. You know the rest.
 *
 */

/*#define MEMW_DEBUG 1 */

#include "timetolease.h"
#include <orb/orbit.h>
#include <stdio.h>

#ifdef MEMW_DEBUG
#include "memwatch.h"
#endif

/*
 * This assertion is for server side fatal errors. Although they are fatal errors,
 * and you want to get out of them, you shouldn't bring down the server because of it
 * (it might have some perfectly happy objects in it) and you certainly shouldn't bring
 * down the client...
 *
 * We use a system exception here. Which is sort of cheating.
 */
#define S_ASSERT(soma,assertion) if(assertion) {\
                     SimpleObjectManagerAdaptor_log_message(&soma,\
                     G_LOG_LEVEL_ERROR,\
                     "Assertion failed at %s:%s",__FILE__,__LINE__);\
                     CORBA_exception_set_system(ev,ex_CORBA_UNKNOWN,CORBA_COMPLETED_NO);\
                     return;}

/*
 * You need to call this bascially after *every* function call ;)
 */
#define S_CHECK_AND_RETHROW(ev) if(ev->_major != CORBA_NO_EXCEPTION){return;}


/*
 * SimpleObjectManager is where we keep the information about
 * the log files this is writing to, any context we want, etc
 *
 * It is an opaque type so that we can change in underneath at will
 *
 */
typedef struct SimpleObjectManager_struct SimpleObjectManager;

/*
 * SimpleObjectManagerAdapter is the per-instance copy of a "smart
 * pointer" to the resources required for Object Management. This
 * includes the time-to-lease manager, the simple object manager
 * and flags about logging messages on a per object basis.
 *
 * This is not opaque, as objects will need access to these structures
 * themselves, to register with the time-to-lease manager and to switch
 * error logging on and off as desired.
 *
 * The structure is delibrately a lightweight structure, so you should be
 * able to copy/assign by value. By doing this, it allows further attributes
 * to be added to the Adaptor without breaking stuff.
 */

typedef struct SimpleObjectManagerAdapter_struct {
  SimpleObjectManager * som;
  int log_flags;
  CORBA_Environment * ev;
  char * unique_id;
  int (*remove_data)(gpointer data);
  gpointer data;
  CORBA_Object corba_object;
} SimpleObjectManagerAdaptor;

/*
 * This makes a new SimpleObjectManager
 *   FILE * log_file - log file to send things to
 *   use_syslog      - use system log not the log file
 *                   (log_file has to be non null *or* use_syslog == 1)
 *   log_flags       - default flags for logging messages. We reuse the Glib error types.
 *   object_lifetime - default life time for objects
 *   context         - optional additional context for the error messages. Allows different somas
 *                     in one orb
 *   max_objects     - max objects
 *   block_size      - number of objects to reap when over max object size
 *   ev              - The activated CORBA_Environment of the ORB
 *
 */

SimpleObjectManager * new_SimpleObjectManager(FILE * log_file,int use_syslog,int log_flags,int object_lifetime,char * context,int max_objects,int block_size,int allow_cache,CORBA_Environment * ev);

/*
 * Removes the soma
 */
void                  delete_SimpleObjectManager(SimpleObjectManager * som);

/*
 * Provides an adaptor from the SimpleObjectManager. Once you have one adaptor,
 * you can assign by value many new cases
 */

SimpleObjectManagerAdaptor SimpleObjectManager_get_Adaptor(SimpleObjectManager * som);

/*
 * log message. We reuse glib types.
 */

int SimpleObjectManagerAdaptor_log_message(SimpleObjectManagerAdaptor * soma,int message_type,char * fmt,...);

/*
 * Provides the activation of the object with the adaptor
 *
 * Class and unique_id *together* have to be unique.
 *
 * The class can be NULL, in which case it is set as the class "NoClass"
 *
 * Unique id can be NULL, in which case it sets soma->unique_id to be NULL.
 * Otherwise, sets soma->unique_id = g_strdup(unique_id), and registers it
 * with the SimpleObjectManager so you can retrieve it again.
 */

int SimpleObjectManagerAdaptor_activate(SimpleObjectManagerAdaptor * soma,const char * class,const char * unique_id,CORBA_Object obj,gpointer data,int (*remove_object)(gpointer data));

/*
 * This should be called in the main loop, if you want to trigger
 * SOM management, which calls the timetolease manager
 */

int SimpleObjectManager_verify(SimpleObjectManager * som);

SimpleObjectManagerAdaptor * SimpleObjectManagerAdaptor_find_object(SimpleObjectManagerAdaptor * soma,const char * class,const char * unique_id);

CORBA_Object SimpleObjectManagerAdaptor_reactivate(SimpleObjectManagerAdaptor * soma);

#endif








