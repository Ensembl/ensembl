# configuration file for the gene_comparison_script.pl
# based on the config files for the pipeline
#
# Written by Eduardo Eyras
# eae@sanger.ac.uk

# Databases in 120 schema (schema as for the 120 release)
#
# ensembl110_new_schema@ensrv3, path_type='UCSC', genetype = '1'
#
# homo_sapiens_core_120@ensrv3, path_type='UCSC', genetype='ensembl'
#
# chr20_120@ecs1b, path_type='Sanger_02', genetypes='HUMACE-Novel_CDS','HUMACE-Known',etc...
#


BEGIN {
package main;

%db1_conf = (
	     'host'      => "ecs1b",

	     'dbname'    => "chr20_120",

	     'path'      => "Sanger_02",

	     'user'      => "ensro",
	     
	     # genetypes should be array_ref, to allow for multiple types
	     'genetypes'  => ["HUMACE-Novel_CDS","HUMACE-Known"],
	     
);

	      
%db2_conf = (
	     'host'      => "ensrv3",

	     'dbname'    => "homo_sapiens_core_120",
	     	     
	     'path'      => "UCSC",
	    	    
	     'user'      => "ensro",

	     'genetypes'  => ["ensembl"], 
	     
);

}

1;
