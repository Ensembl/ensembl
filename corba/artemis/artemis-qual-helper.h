
#ifndef ARTEMIS_QUAL_HELPER_HEADER
#define ARTEMIS_QUAL_HELPER_HEADER

#include "artemis.h"

/*
 * This makes a new single key/value qualifier pair. Don't forget
 * this are structures which get passed directly into the ORB
 */
Ensembl_artemis_Qualifier new_EA_Qualifier(char * key,char * value);

#endif
