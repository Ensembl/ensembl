#
# BioPerl module for Bio::EnsEMBL::KnownMapping::MappingConf;
#
# Cared for by EnsEMBL (ensembl-dev@ebi.ac.uk)
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::KnownMapping::MappingConf - imports global variables used by EnsEMBL gene building

=head1 SYNOPSIS

    use Bio::EnsEMBL::KnownMapping::MappingConf;
    use Bio::EnsEMBL::KnownMapping::MappingConf qw(  );

=head1 DESCRIPTION

MappingConf is a pure ripoff of humConf written by James Gilbert.

humConf is based upon ideas from the standard perl Env environment
module.

It imports and sets a number of standard global variables into the
calling package, which are used in many scripts in the human sequence
analysis system.  The variables are first decalared using "use vars",
so that it can be used when "use strict" is in use in the calling
script.  Without arguments all the standard variables are set, and
with a list, only those variables whose names are provided are set.
The module will die if a variable which doesn\'t appear in its
C<%MappingConf> hash is asked to be set.

The variables can also be references to arrays or hashes.

Edit C<%MappingConf> to add or alter variables.

All the variables are in capitals, so that they resemble environment
variables.

=head1 CONTACT

=cut


package Bio::EnsEMBL::KnownMapping::MappingConf;

use strict;
use vars qw( %MappingConf );

# Hash containing config info
%MappingConf = (

             #################
             #General options#
             #################

             #CHECK      => 'yes',
             CHECK        => 'no',

             #QUERY_IDT   => 50,
             QUERY_IDT    => 40,

             #TARGET_IDT  => 50,
             TARGET_IDT    => 20,

             ################################ 
	     # Files location (Input/Output)#
             ################################


             #Location of the query peptide file (eg: Ensembl predicted protein) 
             #QUERY        => '/work1/mongin/mapping/primary/ensembl110.pep',
             QUERY       => '/acari/work4/mongin/anopheles_gambiae/mapping/Primary/common_prot_dataset_dumped.fa',   
             

             #Location of the sptr file in fasta format containing the entries specific to the organism
	     #'sptr_fa'      => '/work1/mongin/mapping/primary/HS.f',
	     SPTR_FASTA      => '/acari/work4/mongin/anopheles_gambiae/mapping/Primary/ano_knw_formated.fa',
	     
             #Location of the sptr file in Swiss-Prot format containing the entries specific to the organism
	     #'sptr_swiss'      => '/work1/mongin/mapping/primary/HS.SPTR',
	     SPTR_SWISS      => '/acari/work4/mongin/anopheles_gambiae/mapping/Primary/ano_knw.swiss',
	     	     
             #Location of the file containing all refseq and all SP in fasta format (This file will be produced by runni             ng prepare_proteome.pl)
             #'pmatch_input_fa'    => '/work1/mongin/mapping/kate/refseq_p.fa',
	     KNOWN_PROTEINS_FASTA    => '/acari/work4/mongin/anopheles_gambiae/mapping/Primary/total.fa',

             #Output file containing the mapping of SP and refseq sequences to external databases
             #'x_map'  => '/work1/mongin/mapping/outputs/xmap_out1.txt',
             'x_map_out'  => '/acari/work4/mongin/anopheles_gambiae/mapping/Output/x_map.out',

#             'pmatch_out'  => '/acari/work4/mongin/NCBI_26/mapping/Output/pmatch_out_test.txt',

             #Output file from pmatch.pl and input file for maps2db.pl
             #'human_map'  => '/work1/mongin/mapping/outputs/pmatch_human1.txt',
             PMATCH_OUT  => '/acari/work4/mongin/anopheles_gambiae/mapping/Output/pmatch.out',


             #Location of the Refseq (proteins) file in fasta format
	     #'refseq_fa'    => '/work1/mongin/mapping/primary/refseq.fa',
	     REFSEQ_FASTA    => '',
	     
             #Location of the Refseq (proteins) file in Genbank format
	     #'refseq_gnp'    => '/work1/mongin/mouse/mapping/primary/mouse.gnp',
	     REFSEQ_GNP  => '',


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
                  #'refseq_pred_gnp' => '',
                  'refseq_pred_gnp' => '',

                  
                  #Location of the file in .gnp format for the NCBI prediction
                  #'refseq_pred_fa' => '',
                  'refseq_pred_fa' => '',

                  #Location of the file for GO mapping (gene_association.goa)
                  #'go' => '',
                  'go' => '',

                  #######
                  #Mouse#
                  #######

                  #The files needed for the mouse X_mapping can be obatained there: ftp://ftp.informatics.jax.org/pub/informatics/reports/   
                  #2 files are needed MRK_SwissProt.rpt and MRK_LocusLink.rpt
                  
                   #File containing MGI/SP mapping (MRK_SwissProt.rpt)
                   #'mgi_sp'  => '/primary/MRK_SwissProt.rpt',                 
                   'mgi_sp'  => '',
                  
                   #File containing MGI/LocusLink mapping (MRK_LocusLink.rpt)
                   #'mgi_locus'  => '/primary/MRK_LocusLink.rpt',                   
                   'mgi_locus'  => '',
                                      

                   ###########
                   #Anopheles#
                   ###########

                   'celera_map' => '/nfs/acari/rust/anopheles_gambiae/data/celera_transcripts_mapping.dat',


             ###################
             #Database handling#
             ###################

             #DB name
             #'db' => 'proteintest',
             DBNAME => 'anopheles_gambiae_core_4_1',

             #Host name
             #'host' => 'ecs1d',
             DBHOST => 'ecs1d',

             #User
             DBUSER => 'ecs1dadmin',

             #Password
             PASSWORD => 'TyhRv',
             
             #####################
             #Executable location#
             #####################

             #Location for pmatch binaries
             #'pmatch' => '/nfs/disk65/ms2/bin/pmatch'
             PMATCH => '/nfs/disk65/ms2/bin/pmatch',

             ##############################
             #Organism related information#
             ##############################

             #Name of the organism studied. Current keywords used(or planned to be used): human, drosophila, mouse
             #You can adapt the other scripts given the organisms (eg: do some specific x_mapping for a given organism)
             #'organism' => 'human'
             ORGANISM => 'anopheles',
             
);


sub import {
    my ($callpack) = caller(0); # Name of the calling package
    my $pack = shift; # Need to move package off @_

    # Get list of variables supplied, or else
    # all of GeneConf:
    my @vars = @_ ? @_ : keys( %MappingConf );
    return unless @vars;

    # Predeclare global variables in calling package
    eval "package $callpack; use vars qw("
         . join(' ', map { '$'.$_ } @vars) . ")";
    die $@ if $@;


    foreach (@vars) {
	if ( defined $MappingConf{ $_ } ) {
            no strict 'refs';
	    # Exporter does a similar job to the following
	    # statement, but for function names, not
	    # scalar variables:
	    *{"${callpack}::$_"} = \$MappingConf{ $_ };
	} else {
	    die "Error: MappingConf: $_ not known\n";
	}
    }
}

1;

