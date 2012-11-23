#!/usr/bin/env perl

use strict;
use warnings;

use FindBin qw($Bin);
use lib "$Bin";

#Bring in pipeline as well for their helpers
use lib "$Bin/../../../ensembl-pipeline/scripts/Finished/assembly";
use lib "$Bin/../../../ensembl-analysis/modules";

#runtime include normally
require AssemblyMapper::Support;
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Pod::Usage;
use Bio::EnsEMBL::DBSQL::OntologyDBAdaptor;
use Bio::EnsEMBL::Utils::BiotypeMapper;

#Genebuilder utils
require Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils;
Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils->import(qw/clone_Transcript/);
require Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonUtils;
Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonUtils->import(qw/clone_Exon/);

my $support = AssemblyMapper::Support->new(
  extra_options => [
    qw/
      check_transcripts!
      check_exons!
      /
  ]
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
$support->log("EXON KEY       : !! = Very bad (pc mismatch), %% = Somewhat bad (mismatch), ?? = No mapping, might be bad, ¤¤ = eval error\n");
$support->log("TRANSCRIPT KEY : @@ = Very bad (pc translation mismatch), ££ = Very bad (pc transcript mismatch), ** = Somewhat bad (mismatch), XX = No mapping, might be bad, ±± = eval error\n");

my $total_exons = 0;
my $total_transcripts = 0;

$support->iterate_chromosomes(
  prev_stage => '40-fix_overlaps',
  this_stage => '41-conservation',
  worker     => \&compare,
);

$support->log(sprintf("Total exons : %d\n", $total_exons));
$support->log(sprintf("Total transcripts : %d\n", $total_transcripts));

$support->log_stamped("Finished.\n");

sub compare {
  my ($asp) = @_;
  if ($support->param('check_exons')) {
    compare_exons($asp);
  }
  if ($support->param('check_transcripts')) {
    compare_transcripts($asp);
  }
  return 1;
}

sub compare_exons {
  my ($asp) = @_;

  my $R_chr = $asp->ref_chr;
  my $A_chr = $asp->alt_chr;

  my $R_slice = $asp->ref_slice;
  my $A_slice = $asp->alt_slice;

  my $new_slice_adaptor = $R_slice->adaptor();

  my $old_exons = $A_slice->get_all_Exons;

  $support->log(sprintf("Total exons %d\n", scalar(@{$old_exons})));

  while (my $old_exon = shift @$old_exons) {
    $total_exons++;
    eval {exon($old_exon, $R_slice, $new_slice_adaptor)};
    if($@) {
      $support->log(sprintf("¤¤ | ID : %s | EVAL ERROR (%s)\n", $old_exon->stable_id(), $@));
    }
  }
  return;
}

sub exon {
  my ($old_exon, $R_slice, $new_slice_adaptor) = @_;

  # Establish equivalent locations on old and new DBs
  my $new_A_slice = new_Slice($old_exon, $new_slice_adaptor);

  # make a shadow exon for the new database
  my $shadow_exon = clone_Exon($old_exon);
  $shadow_exon->slice($new_A_slice);

  # project new exon to new assembly
  my $projected_exon = $shadow_exon->transform($R_slice->coord_system->name, $R_slice->coord_system->version, $R_slice);

  # Note that Anacode database naming patterns interfere with normal Registry adaptor fetching,
  # hence we must go around the houses somewhat when looking for the appropriate source gene.
  #warn "!! fetching gene adaptor ".$old_exon->adaptor->species.":".$old_exon->adaptor->dbc->dbname."Gene";
  my $old_gene_adaptor = $old_exon->adaptor->db->get_GeneAdaptor();
  my $gene_list = $old_gene_adaptor->fetch_nearest_Gene_by_Feature($old_exon, undef, undef);
  if (scalar(@$gene_list) > 1) { warn "Encountered many possible genes for the exon." }
  my $parent_gene = $gene_list->[0];

  my $state;
  my $location;
  my $difference;
  my $total_length;

  # compare sequences if a projection exists
  if ($projected_exon) {
    my $old_seq       = $old_exon->seq()->seq();
    my $projected_seq = $projected_exon->seq()->seq();
    $total_length = length($old_seq);

    if ($projected_seq ne $old_seq) {

      # Now we have a problem - the feature's sequence was not conserved between assemblies.
      # Determine severity of the problem
      $difference = diff(\$old_seq, \$projected_seq);

      my $group_list = $biotype_mapper->belongs_to_groups($parent_gene->biotype);
      foreach my $group (@$group_list) {
        if ($group eq 'protein_coding') {
          $state = '!!';
          last;
        }
      }
      if (!$state) {

        # Middle badness.
        $state = '%%';
      }
    }

    $location = sprintf('%d : %d', $projected_exon->start(), $projected_exon->end());
  }

  if (!$projected_exon) {
    $state        = '??';
    $location     = 'None';
    $total_length = 0;
    $difference   = 0;
  }

  if ($state) {
    $support->log(sprintf('%s | ID : %s | Gene ID: %s | Location : %s | Differences : %d | Length : %d', $state, $old_exon->stable_id(), $parent_gene->stable_id(), $location, $difference, $total_length) . "\n");
  }

  return;
}

sub compare_transcripts {
  my ($asp)             = @_;
  my $R_slice           = $asp->ref_slice;
  my $A_slice           = $asp->alt_slice;
  my $new_slice_adaptor = $R_slice->adaptor();

  my $old_transcripts = $A_slice->get_all_Transcripts();
  $support->log(sprintf("Total transcripts %d\n", scalar(@{$old_transcripts})));
  while (my $old_transcript = shift @$old_transcripts) {
    $total_transcripts++;
    eval {transcript($old_transcript, $R_slice, $new_slice_adaptor)};
    if($@) {
      $support->log(sprintf("±± | ID : %s | EVAL ERROR (%s)\n", $old_transcript->stable_id(), $@));
    }
  }

  return;
}

sub transcript {
  my ($old_transcript, $R_slice, $new_slice_adaptor) = @_;

  my $new_A_slice = new_Slice($old_transcript, $new_slice_adaptor);
  my $shadow_transcript = clone_Transcript($old_transcript);
  foreach my $e (@{$shadow_transcript->get_all_Exons()}) {
    $e->slice($new_A_slice);
  }
  foreach my $se (@{$shadow_transcript->get_all_supporting_features()}) {
    $se->slice($new_A_slice);
  }
  $shadow_transcript->slice($new_A_slice);
  my $projected_transcript = $shadow_transcript->transform($R_slice->coord_system->name, $R_slice->coord_system->version, $R_slice);

  my $state;
  my $location;
  my $total_length;
  my $difference;
  if ($projected_transcript) {

    #Check if it was protein coding
    my $group_list = $biotype_mapper->belongs_to_groups($projected_transcript->biotype);
    my $is_pc      = 0;
    foreach my $group (@$group_list) {
      if ($group eq 'protein_coding') {
        $is_pc = 1;
        last;
      }
    }

    #Now check for protein sequence mis-match
    if ($is_pc) {
      my $old_seq = $old_transcript->translate()->seq();
      my $new_seq = $projected_transcript->translate()->seq();
      $total_length = length($old_seq);
      if ($old_seq ne $new_seq) {
        $state = '@@';
        $difference = diff(\$old_seq, \$new_seq);
      }
    }

    if (!$state) {
      my $old_seq = $old_transcript->spliced_seq();
      my $new_seq = $projected_transcript->spliced_seq();
      $total_length = length($old_seq);
      if ($old_seq ne $new_seq) {
        $state = ($is_pc) ? '££' : '**';
        $difference = diff(\$old_seq, \$new_seq);
      }
    }

    $location = sprintf('%d : %d', $projected_transcript->start(), $projected_transcript->end());
  }
  else {
    $state        = '**';
    $location     = 'None';
    $total_length = 0;
    $difference   = 0;
  }

  if ($state) {
    $support->log(sprintf('%s | ID : %s | Location : %s | Differences : %d | Length : %d', $state, $old_transcript->stable_id(), $location, $difference, $total_length) . "\n");
  }

  return;
}

sub new_Slice {
  my ($feature, $new_slice_adaptor) = @_;
  my $old_slice    = $feature->slice();
  my $coord_system = $old_slice->coord_system;
  my $new_slice    = $new_slice_adaptor->fetch_by_region($coord_system->name, $old_slice->seq_region_name, $old_slice->start, $old_slice->end, $old_slice->strand, $coord_system->version);
  return $new_slice;
}

#Count number of diffs between strings
sub diff {
  my ($s1, $s2) = @_;
  my $diff  = ${$s1} ^ ${$s2};
  my $total = 0;
  $total += ($+[1] - $-[1]) while $diff =~ m{ ([^\x00]+) }xmsg;
  return $total;
}

__END__

=pod 

=head1 NAME

exon_conservation_check.pl

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
    
    --check_exons                       Test exons for their robustness 
                                        between assembly mappings
                                        
    --check_transcripts                 Expand the test to check if we still 
                                        get full length transcript encoding
=cut
