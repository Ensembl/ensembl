A master copy of the misc_set table is now kept in the production
database.

The production database is kept on the ens-staging1 server and is called
"ensembl_production".

The master table for the external_db table is called "master_misc_set"
and should be copied as-is to any new Core or Core-like database.

For new entries, please insert them into this master table (after
notifying the current release coordinator).


For the release coordinator:

In misc-scripts/production_database/scripts lives a script called
"push_master_tables.pl".  This script is used to synchronise all master
tables over all Core and Core-like databases on all staging servers.

Please run that script with its "--about" and "--help" command line
switches for more information.


$Id$
