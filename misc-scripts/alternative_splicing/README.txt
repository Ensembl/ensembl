1) create the gff file (takes a while)


use ensembl-personal/gk4/altSplicing/scripts/Fetch_gff.pl,  i.e.

perl Fetch_gff.pl -dbhost host1 -dbuser ro -dbname ianl_homo_sapiens_core_55_37 > ianl_homo_sapiens_core_55_37.gff


2) calculate the alternatives and load them into core database.


cat ianl_homo_sapiens_core_55_37.gff |
 ~/ensembl-live/ensembl-personal/gk4/altSplicing/scripts/Find_events.sh | 
perl load_alt_slice_gff.pl -user admin -pass XXX -dbname ianl_homo_sapiens_core_55_37 -host host1





