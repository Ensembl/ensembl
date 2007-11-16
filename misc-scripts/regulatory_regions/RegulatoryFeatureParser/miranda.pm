package RegulatoryFeatureParser::miranda;

use strict;

# Parse data from miRanda analyses; format:
#  #<GROUP>	<SEQ>	<METHOD>	<FEATURE>	<CHR>	<START>	<END>	<STRAND>	<PHASE>	<SCORE>	
#  Similarity	hsa-miR-23b	miRanda	miRNA_target	1	919788	919807	+	.	69	transcript id "ENST00000310998"
#  Similarity	hsa-miR-23a	miRanda	miRNA_target	1	919787	919807	+	.	71	transcript id "ENST00000310998"

use RegulatoryFeatureParser::BaseParser;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

use vars qw(@ISA);
@ISA = qw(RegulatoryFeatureParser::BaseParser);

# Parse file and return hashref containing:
#
# - arrayref of features
# - arrayref of factors

sub parse {

  my ($self, $db_adaptor, $file, $old_assembly, $new_assembly) = @_;

  my %result;

  print "Parsing $file with miranda parser\n";

  my $feature_internal_id = ($self->find_max_id($db_adaptor, "regulatory_feature")) + 1;
  my $highest_factor_id = ($self->find_max_id($db_adaptor, "regulatory_factor")) + 1;

  my $analysis_adaptor = $db_adaptor->get_AnalysisAdaptor();
  my $slice_adaptor = $db_adaptor->get_SliceAdaptor();

  my @features;
  my @factors;
  my %factor_ids_by_name; # name -> factor_id
  my %feature_objects;
  my %anal;

  # TODO - regulatory_factor_coding

  my $stable_id_to_internal_id = $self->build_stable_id_cache($db_adaptor);

  # this object is only used for projection
  my $dummy_analysis = new Bio::EnsEMBL::Analysis(-logic_name => 'CisRedProjection');

  open (FILE, "<$file") || die "Can't open $file";

  while (<FILE>) {

    next if ($_ =~ /^\s*\#/ || $_ =~ /^\s*$/);

    my %feature;

    my ($group, $seq, $method, $feature, $chr, $start, $end, $str, $phase, $score, $pvalue, $type, $id_ignore, $id) = split;
    my $strand = ($str =~ /\+/ ? 1 : -1);
    $id =~ s/[\"\']//g;  # strip quotes

    # ----------------------------------------
    # Feature name

    # For miRNA_target, individual features don't have a unique name, so create
    # a composite one. Also set influence.

    $feature{NAME} = $id .":" . $seq;
    $feature{INFLUENCE} = "negative";

    # ----------------------------------------
    # Factor

    # $seq is the name of a factor - if it's already there, find its ID, otherwise add it
    my $factor_id = $factor_ids_by_name{$seq};
    if (!$factor_id) {
      my %factor;
      $factor_id = $highest_factor_id + 1;
      $factor{INTERNAL_ID} = $factor_id;
      $factor{NAME} = $seq;
      $factor{TYPE} = $feature; # TODO - error checking that $feature is one of the enums?
      push @factors, \%factor;
      $factor_ids_by_name{$seq} = $factor_id;
      $highest_factor_id = $factor_id;
      #print join("  ", ("Factor: ", $factor{ID}, $factor{NAME}, $factor{TYPE})) . "\n";
    }

    $feature{FACTOR_ID} = $factor_id;

    # ----------------------------------------
    # Analysis

    my $analysis = $anal{$method};
    if (!$analysis) {
      $analysis = $analysis_adaptor->fetch_by_logic_name($method);
      $anal{$method} = $analysis;
    }

    if (!$analysis) {
      print STDERR "Can't get analysis for $method, skipping\n";
      next;
    }

    $feature{ANALYSIS_ID} = $analysis->dbID();

    # ----------------------------------------
    # Seq_region ID and co-ordinates

    my $chr_slice = $slice_adaptor->fetch_by_region('chromosome', $chr, undef, undef, undef, $old_assembly);

    if (!$chr_slice) {
      print STDERR "Can't get slice for $chr:$start:$end:$strand\n";
      next;
    }

    my $seq_region_id = $slice_adaptor->get_seq_region_id($chr_slice);

    if (!$seq_region_id) {
      print STDERR "Can't get seq_region_id for chromosome $chr\n";
      next;
    }

    $feature{SEQ_REGION_ID} = $seq_region_id;

    # project if necessary
    if ($new_assembly) {

      #print join("\t", "OLD: ", $start, $end, $strand, $chr, $feature{NAME}) . "\n";

      my $projected_feature = $self->project_feature($start, $end, $strand, $chr, $chr_slice, $dummy_analysis, $new_assembly, $slice_adaptor, $feature{NAME});

      $start = $projected_feature->start();
      $end = $projected_feature->end();
      $strand = $projected_feature->strand();

      $feature{SEQ_REGION_ID} = $slice_adaptor->get_seq_region_id($projected_feature->feature_Slice);

      #print join("\t", "NEW: ", $start, $end, $strand, $chr, $feature{NAME}) . "\n";

    }

    $feature{START} = $start;
    $feature{END} = $end;
    $feature{STRAND} = $strand;

    # ----------------------------------------
    # Ensembl object

    my ($ensembl_type) = $type =~ /(gene|transcript|translation)/;
    if (!$ensembl_type) {
      print STDERR "Can't get ensembl type from $type, skipping\n";
      next;
    }

    my $ensembl_id = $stable_id_to_internal_id->{$ensembl_type}->{$id};
    $ensembl_type = ucfirst(lc($ensembl_type));

    if (!$ensembl_id) {
      print STDERR "Can't get ensembl internal ID for $id, skipping\n";
      next;
    }

    $feature{ENSEMBL_TYPE} = $ensembl_type;
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
  bless $self, "RegulatoryFeatureParser::miranda";
  return $self;

}

1;
