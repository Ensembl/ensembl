=head1 mappping_conf.pl

=head2 Description

This script gives the basic configuration needed by the mapping. This configuration script will be used by each script of the mapping. Each field is described bellow.

    + organism: The name of the organism which has to be mapped. eg: human, mouse, worm

    + sptr:     Location of the SPTR file

    + refseq:   Location of the Refseq file (only valide for the human organism)

    + (ens1, ens2, ens5) : Various files produced by Hugo (http://www.gene.ucl.ac.uk/public-files/nomen/) which would be used to map SP to hugo and refseq to Hugo (only valide for the human organism)

    

=cut





BEGIN {
package main;

%mapping_conf = ( 
	     # Files location (Input/Output)

             'query'        => '/work4/mongin/data/ensembl.pep',
             #'query'       => '/work1/mongin/mapping/targetted/total_targetted.pep',   
        
	     'sptr_fa'      => '/work1/mongin/mapping/primary/HS.f',
	     #'sptr_fa'      => '',
	     
	     'sptr_swiss'      => '/work1/mongin/mapping/primary/HS.SPTR',
	     #'sptr_swiss'      => '',
	     
	     'refseq_fa'    => '/work1/mongin/mapping/primary/refseq.fa',
	     #'refseq'    => '',
	     
	     'refseq_gnp'    => '/work1/mongin/mapping/primary/refseq.gnp',
	     #'refseq_gnp'    => '',
	    
             #File containing all refseq and all SP in fasta format
              'human_fa'    => '/work1/mongin/mapping/primary/h_prot.fa',
	     #'human_fa'    => '',

             
	     'ens1'      => '/work1/mongin/mapping/primary/ens1.txt',
	     #'ens1'      => '',

	     'ens4'      => '/work1/mongin/mapping/primary/ens4.txt',
	     # 'ens4'      => '',

             'refseq_map' => '/work1/mongin/mapping/tests/refseq_map.txt',
             #'refseq_map' => '',

             'sp_map' => '/work1/mongin/mapping/tests/refseq_map3.txt',
             #'sp_map' => '',
             
             #Output file containing the mapping of SP and refseq sequences to external databases
             'x_map'  => '/work1/mongin/mapping/tests/xmap_out1.txt',
             #'x_map_out'  => '',

             #Output file from pmatch.pl and input file for maps2db.pl
             'human_map'  => '/work1/mongin/mapping/tests/human_out1.txt',
             #'x_map_out'  => '',

             #Database handling

             #DB name
             'db' => 'prot_pipeline_test',
             #'db' => '',

             #Host name
             'host' => 'ecs1b',
             #'host' => '',


             #Location for pmatch
             'pmatch' => '/nfs/disk65/ms2/bin/pmatch'

 );


}

1;
