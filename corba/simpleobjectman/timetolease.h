
#ifndef TIMETOLEASE_HEADER
#define TIMETOLEASE_HEADER

/*
 * This file is GPL'd. No warranty. And it ain't brilliant either
 *
 * Bugs and credit to Ewan Birney <birney@sanger.ac.uk>
 *
 */


/*
 * The overview of this package.
 *
 * The idea is that you want to store some objects with a finite "life-span"
 * on them. When you store each object, you provide
 *   a) the data for the object
 *   b) the routine to run on that data to kill it 
 *   c) the life time in seconds for this object
 *   d) an (optional) unique id so you can zap specific objects if you wish
 *
 * You give these things to the TOL manager, and then it does its
 * stuff, reaping the object when it is past its life time date.
 *
 * The TOL manager is *not* multithreaded. Instead you have to call
 * TOL_verify which is the routine which actually checks each active
 * object as to whether it needs reaping or not. Therefore you need to
 * figure out where you will call TOL_verify during whatever you do. A
 * common time to call TOL_verify is when you add new objects (while
 * you are not adding objects, presumably the memory/resources load is
 * ok). To help you out, there is a function which verifies and adds a
 * new object.
 * 
 */

/*
 * BUGS/KNOWN problems etc
 *
 * a) I am not sure if I am using the glist stuff correctly
 *
 * b) glist should be kept ordered by time of death, so not all objects need to
 *    be looked at.
 *
 * c) Not thread safe.
 *
 */


#include "glib.h"
#include <time.h>
#ifdef MEMW_DEBUG
#include "memwatch.h"
#endif

/*
 * TOLManager is an opaque object that represents
 * the holder for a set of objects
 *
 */

typedef struct _TOLManager TOLManager;

/*
 * Makes a new TOLManager
 */

TOLManager * new_TOLManager(void);

/*
 * Removes a TOLManager. Does not remove the data inside
 * probably should do
 */
void         delete_TOLManager(TOLManager * tolm);

/*
 * Both verifies and adds a single data. The remove_data function will
 * be called past its lifetime, and should return 0 on success and non
 * zero on failure. Failure will cause a g_error to be called 
 */

int TOL_add(TOLManager * tolm,gpointer data,int (*remove_data)(gpointer),int life_in_seconds,char * unique_id);

/*
 * Only adds a new data point
 */

int TOL_add_only(TOLManager * tolm,gpointer data,int (*remove_data)(gpointer),int life_in_seconds,char * unique_id);

/*
 * Looks at each object and objects over their lifetime get zapped
 */
int TOL_verify(TOLManager * tolm);


/*
 * Treats the TOL as a first in, first out queue,
 * removing the objects first added to the queue.
 *
 */
int TOL_remove_firstin(TOLManager * tolm,int number);

/*
 * Number of objects in the list
 */
int TOL_list_length(TOLManager * tolm);

/*
 * Removes a unique id. Returns 1 on sucess, 0 on failure
 */
int TOL_remove_unique(TOLManager * tolm,char * unique_id);


#endif







