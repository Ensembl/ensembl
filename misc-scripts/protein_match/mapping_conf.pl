=head1 mappping_conf.pl

=head2 Description

This script gives the basic configuration needed by the mapping. This configuration script will be used by each script of the mapping. 

For some documentatio, see below.

=head2 Contact

Emmanuel Mongin (mongin@ebi.ac.uk)

=cut





BEGIN {
package main;

%mapping_conf = ( 

             ################################ 
	     # Files location (Input/Output)#
             ################################


             #Location of the query peptide file (eg: Ensembl predicted protein) 
             #'query'        => '/work1/mongin/mapping/primary/ensembl110.pep',
             'query'       => '',   
             
             #Location of the sptr file in fasta format
	     #'sptr_fa'      => '/work1/mongin/mapping/primary/HS.f',
	     'sptr_fa'      => '',
	     
             #Location of the sptr file in Swiss-Prot format
	     #'sptr_swiss'      => '/work1/mongin/mapping/primary/HS.SPTR',
	     'sptr_swiss'      => '',
	     
             #Location of the Refseq (proteins) file in fasta format
	     #'refseq_fa'    => '/work1/mongin/mapping/primary/refseq.fa',
	     'refseq'    => '',
	     
             #Location of the Refseq (proteins) file in Genbank format
	     #'refseq_gnp'    => '/work1/mongin/mapping/primary/refseq.gnp',
	     'refseq_gnp'    => '',
	     
             #Location of the file containing all refseq and all SP in fasta format (This file will be produced by runni             ng prepare_proteome.pl)
             #'human_fa'    => '/work1/mongin/mapping/kate/refseq_p.fa',
	     'human_fa'    => '',

             #ens1 and ens4, location of files used for Hugo mapping (http://www.gene.ucl.ac.uk/public-files/nomen/), th             is files will be used only for human
	     #'ens1'      => '/work1/mongin/mapping/primary/ens1.txt',
	     'ens1'      => '',

	     #'ens4'      => '/work1/mongin/mapping/primary/ens4.txt',
	     'ens4'      => '',

                          
             #Output file containing the mapping of SP and refseq sequences to external databases
             #'x_map'  => '/work1/mongin/mapping/outputs/xmap_out1.txt',
             'x_map_out'  => '',

             #Output file from pmatch.pl and input file for maps2db.pl
             #'human_map'  => '/work1/mongin/mapping/outputs/pmatch_human1.txt',
             'x_map_out'  => '',

             ###################
             #Database handling#
             ###################

             #DB name
             #'db' => 'proteintest',
             'db' => '',

             #Host name
             #'host' => 'ecs1d',
             'host' => '',

             #User
             'dbuser' => '',

             #Password
             'password' => '',
             
             #####################
             #Executable location#
             #####################

             #Location for pmatch binaries
             #'pmatch' => '/nfs/disk65/ms2/bin/pmatch'
             'pmatch' => ' ',

             ##############################
             #Organism related information#
             ##############################

             #Name of the organism studied. Current keywords used(or planned to be used): human, drosophila, mouse
             #You can adapt the other scripts given the organisms (eg: do some specific x_mapping for a given organism)
             'organism' => ''
             

 );


}

1;
