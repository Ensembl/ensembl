#include "artemis-db-mysql.h"
#include <stdio.h>
#include <signal.h>
#include <stdlib.h>

#include <popt.h>


char  * host       = "localhost";
char  * user       = "ensemblro";
char  * pass       = "ensemblropass";
char  * db         = "ensdev";

MYSQL mysql;

int verbose        = 0; 

static const
struct poptOption options[] = {
  {"host", 'h', POPT_ARG_STRING, &host, 0, "Host for the MySQL connection", "hostname"},
  {"user", 'u', POPT_ARG_STRING, &user, 0, "User for the MySQL connection", "dbusername"},
  {"pass", 'p', POPT_ARG_STRING, &pass, 0, "Password for the MySQL connection (- means NULL)", "password"},
  {"db",'d',POPT_ARG_STRING, &pass, 0, "Database name for the MySQL connection", "db"},
  {"verbose",'v',POPT_ARG_NONE, &verbose, 0, "Verbose reporting on STDERR", NULL},
  POPT_AUTOHELP
  {NULL, '\0', 0, NULL, 0, NULL, NULL}
};


MYSQL * ea_connect ( void ) 
{
    MYSQL * connection;
    
    if( strcmp(pass,"-") == 0 ) {
	pass = NULL;
    }
    
    mysql_init(&mysql);
    if ( verbose ) {
	fprintf(stderr,"Connecting with %s %s %s\n",host,user,db);
    }
    
    connection = mysql_real_connect(&mysql,host,user,pass,db,0,0,0); 

    if( connection == NULL ) {
	g_error("Unable to make connection to mysql with host %s, user %s, password %s and database %s",host,user,pass == NULL ? "NoPassword" : pass,db);
	exit(0);
    }
    return connection;
}


int main (int argc, char *argv[])
{
    MYSQL * connection;
    char sqlbuffer[1024];
   MYSQL_RES * result;
   MYSQL_ROW row;
   int state;
   
   char * c_id = "Z69666.00001";

  PortableServer_ObjectId objid = {0, sizeof("EnsemblArtemisServer"), "EnsemblArtemisServer"};
  PortableServer_POA poa;
  FILE * ifp;
  CORBA_Environment ev;
  char *retval;
  CORBA_ORB orb;
  Ensembl_artemis_DB eadb;
  
  poptContext pcon;
  int rc;
  
  pcon=poptGetContext("ensembl-artemis-server", argc, argv, options, 0);
  /*  poptSetOtherOptionHelp(pcon, "<IDL files>");*/
  
  if(argc < 1) {
    poptPrintUsage(pcon, stdout, 0);
    return 0;
  }
  
  if((rc=poptGetNextOpt(pcon)) < -1) {
    g_print("ensembl-artemis-server: bad argument %s: %s\n", 
	    poptBadOption(pcon, POPT_BADOPTION_NOALIAS),
	    poptStrerror(rc));
    exit(0);
  }
  
  connection = ea_connect();

  if( verbose ) {
    fprintf(stderr,"Made connection to MySQL successfully...\n");
  }

  CORBA_exception_init(&ev);
  orb = CORBA_ORB_init(&argc, argv, "orbit-local-orb", &ev);
  poa = (PortableServer_POA)CORBA_ORB_resolve_initial_references(orb, "RootPOA", &ev);
  PortableServer_POAManager_activate(PortableServer_POA__get_the_POAManager(poa, &ev), &ev);

  if( verbose ) {
    fprintf(stderr,"Built ORB successfully...\n");
  }


  eadb = new_EA_Database(poa,connection,verbose,&ev);

  if( verbose ) {
    fprintf(stderr,"Built Ensembl Database object...\n");
  }
  
  retval = CORBA_ORB_object_to_string(orb, eadb, &ev);
  ifp = fopen("db.ior","w");
  fprintf(ifp,"%s\n", retval); 
  fclose(ifp);

  if( verbose ) {
    fprintf(stderr,"Written out db.ior file...\n");
  }

  CORBA_free(retval);
  

  if( verbose ) {
    fprintf(stderr,"...Waiting for DB requests\n");
  }
  CORBA_ORB_run(orb, &ev);
  return 0;



}

