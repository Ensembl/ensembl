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
             #'query'       => '',   
        
	     'sptr_fa'      => '/work1/mongin/mapping/kate/refseq_p.fa',
	     #'sptr_fa'      => '',
	     
	     'sptr_swiss'      => '/work1/mongin/mapping/primary/sptr.swiss',
	     #'sptr_swiss'      => '',
	     
	     'refseq_fa'    => '/work1/mongin/mapping/primary/refseq.fa',
	     #'refseq'    => '',
	     
	     'refseq_gnp'    => '/work1/mongin/mapping/primary/refseq.gnp',
	     #'refseq_gnp'    => '',
	    
	     'ens1'      => '/work1/mongin/mapping/primary/ens1.txt',
	     #'ens1'      => '',

	     'ens4'      => '/work1/mongin/mapping/primary/ens4.txt',
	     # 'ens4'      => '',

             'refseq_map' => '/work1/mongin/mapping/tests/refseq_map.txt',
             #'refseq_map' => '',

             'sp_map' => '/work1/mongin/mapping/tests/refseq_map3.txt',
             #'sp_map' => '',

             'x_map'  => '/work1/mongin/mapping/tests/xmap_out1.txt',
             #'x_map_out'  => '',

             #Database handling
             'db' => 'prot_pipeline_test',
             #'db' => '',

             'host' => 'ecs1b',
             #'host' => '',

             

 );


}

1;
