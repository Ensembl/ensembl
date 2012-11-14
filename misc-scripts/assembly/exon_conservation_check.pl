#!/usr/bin/env perl

=head1 SUMMARY

This script is intended to highlight issues with an assembly mapping, by inspecting
the equivalent sequence for each exon. The resulting log is grep-suitable and keyed
for severity.

=head1 SYNOPSIS

perl exon_conservation_check.pl <many arguments>

    --dbname, db_name=NAME              database name NAME
    --host, --dbhost, --db_host=HOST    database host HOST
    --port, --dbport, --db_port=PORT    database port PORT
    --user, --dbuser, --db_user=USER    database username USER
    --pass, --dbpass, --db_pass=PASS    database passwort PASS
    --assembly=ASSEMBLY                 assembly version ASSEMBLY
    --altdbname=NAME                    alternative database NAME
    --altassembly=ASSEMBLY              alternative assembly version ASSEMBLY

Optional options
    --logfile, --log=FILE               log to FILE (default: *STDOUT)
    --logpath=PATH                      write logfile to PATH (default: .)
=cut

use strict;
use warnings;

use FindBin qw($Bin);
use lib "$Bin";

use AssemblyMapper::Support;
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Pod::Usage;
use Bio::EnsEMBL::DBSQL::OntologyDBAdaptor;
use Bio::EnsEMBL::Utils::BiotypeMapper;

my $support = AssemblyMapper::Support->new(
   
);
unless ($support->parse_arguments(@_)) {
    warn $support->error if $support->error;
    pod2usage(1);
}
$support->connect_dbs;


my $onto_db_adaptor = Bio::EnsEMBL::DBSQL::OntologyDBAdaptor->new( 
    -DBNAME => $support->ref_dba->dbc->dbname,
    -DBCONN => $support->ref_dba->dbc,
);
my $biotype_mapper = new Bio::EnsEMBL::Utils::BiotypeMapper($onto_db_adaptor);

$support->log_stamped("Beginning analysis.\n");
$support->log("!! = Very bad, %% = Somewhat bad, ?? = No mapping, might be bad\n");
my $ok = $support->iterate_chromosomes(
    prev_stage => '40-fix_overlaps',
    this_stage =>     '41-exon-conservation',
    worker =>         \&compare_exons,
);

$support->log_stamped("Finished.\n");

sub compare_exons {
    my ($asp) = @_;
    
    my $R_chr   = $asp->ref_chr;
    my $A_chr   = $asp->alt_chr;

    my $R_slice = $asp->ref_slice;
    my $A_slice = $asp->alt_slice;
    
    my $old_exons = $A_slice->get_all_Exons;
    while (my $old_exon = shift @$old_exons) {
        # Establish equivalent locations on old and new DBs
        my $coord_system = $old_exon->slice->coord_system;

        my $new_slice = $R_slice->adaptor->fetch_by_region(
            $coord_system->name,
            $old_exon->seq_region_name,
            $old_exon->start,
            $old_exon->end,
            $old_exon->strand,
            $coord_system->version
        );
        
        # make a shadow exon for the new database
        my $shadow_exon = new Bio::EnsEMBL::Feature (
            -start => $old_exon->seq_region_start,
            -end => $old_exon->seq_region_end,
            -strand => $old_exon->strand,
            -slice => $new_slice,
        );
        # project new exon to new assembly
        my $projected_exon = $shadow_exon->transform($A_slice->coord_system->name,$A_slice->coord_system->version,$R_slice);

        # Note that Anacode database naming patterns interfere with normal Registry adaptor fetching,
        # hence we must go around the houses somewhat when looking for the appropriate source gene.
        #warn "!! fetching gene adaptor ".$old_exon->adaptor->species.":".$old_exon->adaptor->dbc->dbname."Gene";
        my $old_gene_adaptor = $old_exon->adaptor->db->get_GeneAdaptor();
        my $gene_list = $old_gene_adaptor->fetch_nearest_Gene_by_Feature($old_exon, undef, undef);
        if (scalar(@$gene_list) >1) {warn "Encountered many possible genes for the exon."}
        my $parent_gene = $gene_list->[0];
        # compare sequences if a projection exists
        if ($projected_exon && $projected_exon->seq ne $old_exon->seq) {
            # Now we have a problem - the feature's sequence was not conserved between assemblies.
            # Determine severity of the problem
            
            my $group_list = $biotype_mapper->belongs_to_groups($parent_gene->biotype);
            my $warned = 0;
            foreach my $group (@$group_list) {
                if ($group eq 'protein_coding') {
                    # Maximum badness.
                    $support->log("!! ".$parent_gene->stable_id." - ".$old_exon->stable_id.
                        " projected - ".$projected_exon->start.":".$projected_exon->end."\n"
                    );
                    $warned = 1;
                }
            }
            unless ($warned) {
                # Middle badness. 
                $support->log("%% ".$parent_gene->stable_id." - ".$old_exon->stable_id."\n");
            }
        }
        
        if (! $projected_exon) {
            # No projection possibly bad news
            $support->log("?? ".$parent_gene->stable_id." - ".$old_exon->stable_id."\n");
        }
        
    }
       
}

