# configuration file for run_GeneComparison
# based on the config files for the pipeline
#
# Written by Eduardo Eyras
# eae@sanger.ac.uk


# Databases in 120 schema (schema as for the 120 release)
#
# prediction release 110: ensembl110_new_schema_2@ecs1c, path_type='UCSC', genetype = '1'

# prediction release 120: homo_sapiens_core_120@ensrv3, path_type='UCSC', genetype='ensembl'
#
# prediction release 121: ens_NCBI_26_temp@ecs1f, path='NCBI_26', genetype = 'ensembl' 
#
# chr20_120@ecs1b, path_type='Sanger_02', genetypes='HUMACE-Novel_CDS','HUMACE-Known',etc...
#
# genetypes is an array_ref, to allow for multiple types

BEGIN {
package main;

# configuration for the annotation (or benchmark) genes

%db1_conf = (
	     ######################################################
	     # chr20 annotation database (schema as release 120 ) #
	     ######################################################
	     #'host'      => "ecs1b",
	     #'dbname'    => "chr20_120",
	     #'path'      => "Sanger_02",   # to be used for predictions on UCSC
	     #'path'      => "NCBI_26",     # for release 121
	     #'user'      => "ensro",
	     #'genetypes'  => ["HUMACE-Novel_CDS","HUMACE-Known"],
	     #'genetypes' => ["HUMACE-Pseudogene"],
	     

	     ######################################################
	     # mouse_4.2.1 gene-build
	     ######################################################
	     'host'      => "ecs1f",
	     'dbname'    => "mouse_Sanger_Nov01_denormalised",
	     'path'      => "CHR",    
	     'user'      => "ensro",
	     'genetypes' => ["ensembl"],
	     
	     #'host'      =>  ,
	     #'dbname'    =>  ,
	     #'path'      =>  ,
	     #'user'      =>  ,
	     #'genetypes' =>  ,
	     
);

# configuration for the prediction genes
	      
%db2_conf = (
	     ##############################
	     # parameters for release 110 #
	     ##############################
	     #'host'      => "ecs1c",                    
	     #'dbname'    => "ensembl110_new_schema_2",
	     #'path'      => "UCSC",                     
	     #'user'      => "ensro",
	     #'genetypes' => ["1"],                  

	     ##############################
	     # parameters for release 120 #
	     ##############################
	     #'host'      => "ensrv3",                   
	     #'dbname'    => "homo_sapiens_core_120",    
	     #'path'      => "UCSC",                     
	     #'user'      => "ensro",
	     #'genetypes' => ["ensembl"],               
	     
	     ##############################
	     # parameters for release 121 #
	     ##############################
	     #'host'      => "ecs1e",                   
	     #'dbname'    => "ens_NCBI_26",              
	     #'path'      => "NCBI_26",                  
	     #'user'      => "ensro",
	     #'genetypes' => ["ensembl"],              
	     
	     ##################################
	     # parameters for genes from ESTs #
	     ##################################
	     #'host'      => "ecs1e",                    
	     #'dbname'    => "ens_UCSC_0801_est90",      
	     #'path'      => "UCSC",                     
	     #'user'      => "ensro",
	     #'genetypes' => ["genomewise_clustering_twice"], 

	     ##################################
	     # Mouse_4.2.1 EST-genes
	     ##################################
	     'host'      => 'ecs1b',                    
	     'dbname'    => 'mouse_Sanger_Nov01_est',
	     'path'      => 'CHR',
	     'user'      => 'ensro',
	     'genetypes' => ["genomewise"], 
	     
	     #'host'      =>  ,
	     #'dbname'    =>  ,
	     #'path'      =>  ,
	     #'user'      =>  ,
	     #'genetypes' =>  ,
);

}

1;
