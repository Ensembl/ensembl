*********************************************************
*				Alternative Splicing Events Computations				*
*																												*
*				Authors:																				*				
*				Ian Longden																			*
*       Gautier Koscielny																* 
*       Please contact Ensembl helpdesk for more				*
*				information.																		*
*																												*
*       Last update: 23/02/2010													*
*																												*
*********************************************************

Alternative splicing (AS) events computation is a three step process.
First, the transcript structures are fetch from the core database and
stored in a text file in GFF.
Then, a binary program called altSpliceFinder will compute all the events
and generate a GFF output of the events.
Finally, the events are uploaded in the core database. 

1) create the GFF file (takes a while)
--------------------------------------

use ensembl/misc-scripts/alternative_splicing/Fetch_gff.pl,  i.e.

perl Fetch_gff.pl -dbhost host1 -dbuser ro -dbname ianl_homo_sapiens_core_55_37 > ianl_homo_sapiens_core_55_37.gff


2) calculate the AS events
--------------------------

2.1 prerequisite

the AltSplicingToolkit must be installed to compute the splicing events.
Please refer to the documentation in 

ensembl/misc-scripts/alternative_splicing/AltSplicingToolkit/INSTALL

to install this toolkit.

2.2 altSpliceFinder

Run the altSpliceFinder binary program from the toolkit

altSpliceFinder -i ianl_homo_sapiens_core_55_37.gff -o ianl_homo_sapiens_core_55_37_AS_events.gff


3) load them into core database
-------------------------------

cat ianl_homo_sapiens_core_55_37_AS_events.gff | \
perl load_alt_splice_gff.pl -user admin -pass XXX -dbname ianl_homo_sapiens_core_55_37 -host host1

--
