#this is a SQL file intended to be run after mapping the assemblies
#once you are happy with the results, run it like this
#
#mysql -h your_host -u the_user -p database_updated < cleanup_tmp_table.sql
#will remove the *bak and *tmp* tables created during the mapping process

DROP TABLE assembly_bak;
DROP TABLE coord_system_bak;
DROP TABLE meta_bak;
DROP TABLE seq_region_bak;
DROP TABLE tmp_align;
