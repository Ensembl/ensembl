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
             
             #Location of the sptr file, this file will be used as an input to grep the specific sp entries to the organism using grep_sp_entries.pl. This file is supposed to be in SP format
             #'total_sptr'  => '/work1/mongin/mapping/primary/sptr.txl',
             'total_sptr'  => '/work1/mongin/mouse/data/old/spMouse.sp',

             #Location of the sptr file in fasta format containing the entries specific to the organism
	     #'sptr_fa'      => '/work1/mongin/mapping/primary/HS.f',
	     'sptr_fa'      => '/work1/mongin/mouse/data/old/tmp.fa',
	     
             #Location of the sptr file in Swiss-Prot format containing the entries specific to the organism
	     #'sptr_swiss'      => '/work1/mongin/mapping/primary/HS.SPTR',
	     'sptr_swiss'      => '/work1/mongin/mouse/data/old/tmp.swiss',
	     
             #Location of the Refseq (proteins) file in fasta format
	     #'refseq_fa'    => '/work1/mongin/mapping/primary/refseq.fa',
	     'refseq'    => '',
	     
             #Location of the Refseq (proteins) file in Genbank format
	     #'refseq_gnp'    => '/work1/mongin/mapping/primary/refseq.gnp',
	     'refseq_gnp'    => '',
	     
             #Location of the file containing all refseq and all SP in fasta format (This file will be produced by runni             ng prepare_proteome.pl)
             #'human_fa'    => '/work1/mongin/mapping/kate/refseq_p.fa',
	     'pmatch_input_fa'    => '',

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
             'pmatch_out'  => '',

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
             'pmatch' => '',

             ##############################
             #Organism related information#
             ##############################

             #Name of the organism studied. Current keywords used(or planned to be used): human, drosophila, mouse
             #You can adapt the other scripts given the organisms (eg: do some specific x_mapping for a given organism)
             #'organism' => 'human'
             'organism' => '',
             

             #OX (Organism taxonomy cross-reference) number
             'ox' => '9606'
             #'ox' => ''

 );


}

1;
