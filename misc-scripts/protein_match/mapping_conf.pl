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

             #################
             #General options#
             #################

             #'check'      => 'yes',
             'check'      => 'no',


             #The mapping to known genes is assymetrical. This is due to the fact that's our gene prediction is quite fragmented compared to the manually curated genes             

             #'query_idt'  => 40,
             'query_idt'    => 40,

             #'target_idt  => 10,
             'target_idt'  => 40,

             #For the mapping to predicted gene, the mapping is better to be more or less symetrical (the values given are indication on what should be used)

             #'query_idt'  => 40,
             'pred_query_idt'    => ,

             #'target_idt  => 50,
             'pred_target_idt'  => ,

             #Location of the statistic file (only neede if you run get_stats.pl)
             #'statistic'  => '/work1/mongin/mapping/stats.txt',
             'statistic_file'  => '',        


             ################################ 
	     # Files location (Input/Output)#
             ################################


             #Location of the query peptide file (eg: Ensembl predicted protein) 
             #'query'        => '/work1/mongin/mapping/primary/ensembl110.pep',
             'query'       => '/acari/work4/mongin/final_build/release_mapping/Primary/final.fa',   
             
             #Location of the sptr file, this file will be used as an input to grep the specific sp entries to the organism using grep_sp_entries.pl. This file is supposed to be in SP format
             #'total_sptr'  => '/work1/mongin/mapping/primary/sptr.txl',
             'total_sptr'  => '',

             #Location of the sptr file in fasta format containing the entries specific to the organism
	     #'sptr_fa'      => '/work1/mongin/mapping/primary/HS.f',
	     'sptr_fa'      => '/acari/work4/mongin/final_build/release_mapping/Primary/sptr_ano_gambiae_19_11_02_formated.fa',
	     
             #Location of the sptr file in Swiss-Prot format containing the entries specific to the organism
	     #'sptr_swiss'      => '/work1/mongin/mapping/primary/HS.SPTR',
	     'sptr_swiss'      => '/acari/work4/mongin/final_build/release_mapping/Primary/sptr_ano_gambiae_19_11_02.swiss',
	     	     
             #Location of the file containing all refseq and all SP in fasta format (This file will be produced by runni             ng prepare_proteome.pl)
             #'human_fa'    => '/work1/mongin/mapping/kate/refseq_p.fa',
	     'pmatch_input_fa'    => '/acari/work4/mongin/final_build/release_mapping/Primary/sptr_ano_gambiae_19_11_02_formated.fa',

             #Output file containing the mapping of SP and refseq sequences to external databases
             #'x_map'  => '/work1/mongin/mapping/outputs/xmap_out1.txt',
             'x_map_out'  => '/acari/work4/mongin/final_build/release_mapping/Output/x_map.out',

             #Output file from pmatch.pl and input file for maps2db.pl
             #'human_map'  => '/work1/mongin/mapping/outputs/pmatch_human1.txt',
             'pmatch_out'  => '/acari/work4/mongin/final_build/release_mapping/Output/scanwise.out',


             #Location of the Refseq (proteins) file in fasta format
	     #'refseq_fa'    => '/work1/mongin/mapping/primary/refseq.fa',
	     'refseq_fa'    => '',
	     
             #Location of the Refseq (proteins) file in Genbank format
	     #'refseq_gnp'    => '/work1/mongin/mouse/mapping/primary/mouse.gnp',
	     'refseq_gnp'  => '',

             ############################################
             #Organism specific files for the X_mapping #
             ############################################
                  
                  #######
                  #Human#
                  #######

                  #ens1 and ens4, location of files used for Hugo mapping (http://www.gene.ucl.ac.uk/public-files/nomen/),                   th is files will be used only for human
	          #'ens1'      => '/work1/mongin/mapping/primary/ens1.txt',
	          'ens1'      => '',

	          #'ens4'      => '/work1/mongin/mapping/primary/ens4.txt',
	          'ens4'      => '',

                  #Location of the file in .gnp format for the NCBI prediction
                  #'refseq_pred' => '',
                  'refseq_pred' => '',

                  #Location of the file for GO mapping (gene_association.goa)
                  #'go' => '',
                  'go' => '',
                  
                  #######
                  #Mouse#
                  #######

                  #The files needed for the mouse X_mapping can be obatained there: ftp://ftp.informatics.jax.org/pub/informatics/reports/   
                  #2 files are needed MRK_SwissProt.rpt and MRK_LocusLink.rpt
                  
                   #File containing MGI/SP mapping (MRK_SwissProt.rpt)
                   #'mgi_sp'  => '/work1/mongin/mouse/mapping/primary/MRK_SwissProt.rpt',                 
                   'mgi_sp'  => '',
                  
                   #File containing MGI/LocusLink mapping (MRK_LocusLink.rpt)
                   #'mgi_locus'  => '/work1/mongin/mouse/mapping/primary/MRK_LocusLink.rpt',                   
                   'mgi_locus'  => '',
                                      


             ###################
             #Database handling#
             ###################

             #DB name
             #'db' => 'proteintest',
             'db' => 'anopheles_arne_core_9_2',

             #Host name
             #'host' => 'ecs1d',
             'host' => 'ecs1b',

             #User
             'dbuser' => 'ensadmin',

             #Password
             'password' => 'ensembl',
             
             #####################
             #Executable location#
             #####################

             #Location for pmatch binaries
             #'pmatch' => '/nfs/disk65/ms2/bin/pmatch'
             'pmatch' => '/nfs/disk65/ms2/bin/pmatch',

             

             ##############################
             #Organism related information#
             ##############################

             #Name of the organism studied. Current keywords used(or planned to be used): human, drosophila, mouse
             #You can adapt the other scripts given the organisms (eg: do some specific x_mapping for a given organism)
             #'organism' => 'human'
             'organism' => 'anopheles',
             

             #OX (Organism taxonomy cross-reference) number
             #'ox' => '9606'
             #'ox' => '10090'
             #'ox' => '7227'
             'ox'  => ''

 );


}

1;





