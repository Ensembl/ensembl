#!/bin/sh

/usr/local/ensembl/bin/perl human2chimp.pl -hdbname homo_sapiens_core_20_34b -hhost ecs4 -huser ensro -hport 3351 -cdbname mcvicker_chimp_human_merge -chost ecs4 -cport 3350 -cuser ensro -hassembly NCBI34 -cassembly BROAD1 -hchromosome  9 -logfile chimp.log
