
#ifndef TIMETOLEASE_PRIVATE
#define TIMETOLEASE_PRIVATE

#include "timetolease.h"

struct _TOLManager {
  GList * list;
};

typedef struct {
  time_t created;
  int    life_time;
  gpointer data;
  int (*remove_data)(gpointer);
  char * unique_id;
} TOLWrapper;


#endif
