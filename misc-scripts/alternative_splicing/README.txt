*********************************************************
*	Alternative Splicing Events Computations              *
*							                                          *
*	Authors:					                                    *
*	  Ian Longden                                         *
*   Gautier Koscielny				                            * 
*   Please contact Ensembl helpdesk for more            *
*	  information.                                        *
*                                                       *
*       Last update: 10/10/2011                         *
*                                                       *
*********************************************************

Table of contents
	0. Preamble
	1. Dump gene models in GFF
	2. Calculate the AS events
	3. Populate database with AS events
	4. All in one: pipeline script used at Sanger Institute

0) Preamble
-----------

Alternative splicing (AS) events computation is a three step process.
First, the transcript structures are fetch from the core database and
stored in a text file in GFF.
Then, a binary program called altSpliceFinder will compute all the events
and generate a GFF output of the events.
Finally, the events are uploaded in the core database. 

There is now a script that does all the job for you for every release of Ensembl. 
If you're interested in computing Alternative splicing events for other species,
you can follow section 1) and 2) of this document.
If you want to run the pipeline from start to end in one go, there is a script 
described in section 4) of this document.
  

1) Dump gene models in GFF 
--------------------------

use ensembl/misc-scripts/alternative_splicing/Fetch_gff.pl,  i.e.

perl Fetch_gff.pl -dbhost host1 -dbuser ro -dbname ianl_homo_sapiens_core_55_37 > ianl_homo_sapiens_core_55_37.gff


2) Calculate the AS events
--------------------------

2.1 prerequisite

the AltSplicingToolkit must be installed to compute the splicing events.
Please refer to the documentation in 

ensembl/misc-scripts/alternative_splicing/AltSplicingToolkit/INSTALL

to install this toolkit.

2.2 altSpliceFinder

Run the altSpliceFinder binary program from the toolkit

altSpliceFinder -i ianl_homo_sapiens_core_55_37.gff -o ianl_homo_sapiens_core_55_37_AS_events.gff


3) Populate database with AS events
-----------------------------------

cat ianl_homo_sapiens_core_55_37_AS_events.gff | \
perl load_alt_splice_gff.pl -user admin -pass XXX -dbname ianl_homo_sapiens_core_55_37 -host host1

4) All in one: pipeline script at Sanger
----------------------------------------

Run the script as_event_computations.sh with the following parameters:

Usage:
    as_event_computations.sh -h <dbhost> [-P <dbport>] -u <dbuser> [-p <dbpass>] 
                             -s <species> [-d <dbname>] [-o <output_dir>]

If the species name is passed, the script will find the corresponding core database on <dbhost>.
If the database name <dbname> is passed, the script will use this database as the core database.
By default, all intermediate results will be written in the /tmp directory.
Please use the -o parameter to pass a different existing writable directory.

--
