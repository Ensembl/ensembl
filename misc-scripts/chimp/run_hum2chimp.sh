#!/bin/sh

/usr/local/ensembl/bin/perl human2chimp.pl  \
-hdbname homo_sapiens_core_26_35 -hchromosome 1 -hhost ecs2 -huser ensro -hport 3364 \
-cdbname pan_troglodytes_core_26_1 -chost ecs2 -cport 3364 -cuser ensro \
-hassembly NCBI35 -cassembly CHIMP1 -logfile chimp.test.log \
-store -duser ensadmin -dpass ensembl -dhost ecs2 -dport 3364 -ddbname pan_troglodytes_core_26_1

