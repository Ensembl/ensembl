#!/bin/sh

/usr/local/ensembl/bin/perl human2chimp.pl  \
-hdbname mcvicker_human_genes -hhost ecs4 -huser ensro -hport 3350 \
-cdbname mcvicker_chimp_agp -chost ecs4 -cport 3350 -cuser ensro \
-hassembly NCBI34 -cassembly BROAD1 -logfile chimp.test.log \
-store -duser ensadmin -dpass ensembl -dhost ecs4 -dport 3350 -ddbname mcvicker_chimp_test

