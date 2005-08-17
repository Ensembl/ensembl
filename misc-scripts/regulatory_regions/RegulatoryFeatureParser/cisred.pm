package RegulatoryFeatureParser::cisred;

use strict;

# Parse data from cisred database dump; download zip file from
# http://www.cisred.org/files/Databases/cisREDdb-Hsap-1.1b/SQL/data.zip
# Only features.txt is parsed
#
# Format:

# <id> <batch_id> <seqname> <source> <feature> <start> <end> <score> <strand> <frame> <gene_id> <consensus>
#
# 1	1	7	cisRed	con.omops.w10_1	75069868	75069877	0.0370026894185291	-	.	ENSG00000006606	rGGTTkGGGG
# 2	1	7	cisRed	con.omops.w6_1	75069913	75069918	0.0124691222715891	-	.	ENSG00000006606	GCCTGG

use RegulatoryFeatureParser::BaseParser;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

use vars qw(@ISA);
@ISA = qw(RegulatoryFeatureParser::BaseParser);

# Parse file and return hashref containing:
#
# - arrayref of features
# - arrayref of factors

sub parse {

  my ($self, $db_adaptor, $file) = @_;

  my %result;

  print "Parsing $file with cisred parser\n";

  my $feature_internal_id = ($self->find_max_id($db_adaptor, "regulatory_feature")) + 1;
  my $highest_factor_id = ($self->find_max_id($db_adaptor, "regulatory_factor")) + 1;

  my $analysis_adaptor = $db_adaptor->get_AnalysisAdaptor();
  my $slice_adaptor = $db_adaptor->get_SliceAdaptor();

  my @features;
  my @factors;
  my %factor_ids_by_name; # name -> factor_id
  my %feature_objects;
  my %anal;

  # TODO - regulatory_factor_transcripts

  my $stable_id_to_internal_id = $self->build_stable_id_cache($db_adaptor);

  open (FILE, "<$file") || die "Can't open $file";

  while (<FILE>) {

    next if ($_ =~ /^\s*\#/ || $_ =~ /^\s*$/);

    my %feature;

    my ($id, $batch_id, $seqname, $source, $feature, $start, $end, $score, $strand, $frame, $gene_id, $consensus) = split;

    my $strand = ($strand =~ /\+/ ? 1 : -1);

    # ----------------------------------------
    # Factor

    # $seq is the name of a factor - if it's already there, find its ID, otherwise add it
    my $factor_id = $factor_ids_by_name{$feature};
    if (!$factor_id) {
      my %factor;
      $factor_id = $highest_factor_id + 1;
      $factor{INTERNAL_ID} = $factor_id;
      $factor{NAME} = $feature;
      $factor{TYPE} = $feature; # TODO - error checking that $feature is one of the enums?
      push @factors, \%factor;
      $factor_ids_by_name{$feature} = $factor_id;
      $highest_factor_id = $factor_id;
      #print join("  ", ("Factor: ", $factor{ID}, $factor{NAME}, $factor{TYPE})) . "\n";
    }

    $feature{FACTOR_ID} = $factor_id;

    # ----------------------------------------
    # Analysis

    my $analysis = $anal{$source};
    if (!$analysis) {
      $analysis = $analysis_adaptor->fetch_by_logic_name($source);
      $anal{$source} = $analysis;
    }

    if (!$analysis) {
      print STDERR "Can't get analysis for $source, skipping\n";
      next;
    }

    $feature{ANALYSIS_ID} = $analysis->dbID();

    # ----------------------------------------
    # Seq_region ID and co-ordinates

    my $seqname_slice = $slice_adaptor->fetch_by_region(undef, $seqname, $start, $end, $strand);

    if (!$seqname_slice) {
      print STDERR "Can't get slice for $seqname:$start:$end:$strand\n";
      next;
    }

    my $seq_region_id = $slice_adaptor->get_seq_region_id($seqname_slice);

    if (!$seq_region_id) {
      print STDERR "Can't get seq_region_id for chromosome $seqname\n";
      next;
    }

    $feature{SEQ_REGION_ID} = $seq_region_id;
    $feature{START} = $seqname_slice->start();
    $feature{END} = $seqname_slice->end();
    $feature{STRAND} = $seqname_slice->strand();

    # ----------------------------------------
    # Feature name

    # For cisRed, individual features don't have a unique name, so create
    # a composite one. Also set influence.

    $feature{NAME} = $gene_id .":" . $feature;
    $feature{INFLUENCE} = "unknown"; # TODO - what does cisRed store?

    # ----------------------------------------
    # Ensembl object - always a gene in cisRed

    my $ensembl_id = $stable_id_to_internal_id->{gene}->{$gene_id};

    if (!$ensembl_id) {
      print STDERR "Can't get ensembl internal ID for $gene_id, skipping\n";
      next;
    }

    $feature{ENSEMBL_TYPE} = "Gene";
    $feature{ENSEMBL_ID} = $ensembl_id;

    # ----------------------------------------
    # Feature internal ID
    # note this is not the "id" referred to above

    $feature{INTERNAL_ID} = $feature_internal_id++;

    # ----------------------------------------
    # Evidence

    $feature{EVIDENCE} = "";

    # ----------------------------------------
    # Add to object to be returned

    push @features, \%feature;

    #print "Feature: ";
    #foreach my $field (keys %feature) {
    #  print $field . ": " . $feature{$field} . " ";
    #}
    #print "\n";

  }

  close FILE;

  $result{FEATURES} = \@features;
  $result{FACTORS} = \@factors;

  print "Parsed " . scalar(@{$result{FEATURES}}) . " features and " . scalar(@{$result{FACTORS}}) . " factors\n";

  return \%result;

}



sub new {

  my $self = {};
  bless $self, "RegulatoryFeatureParser::cisred";
  return $self;

}

1;
