

package RegulatoryFeatureParser::enhancer;

use strict;

# Parse data from LBL enhancers, see http://enhancer.lbl.gov/cgi-bin/imagedb.pl?show=1;search.result=yes;form=search;search.form=no;action=search;search.sequence=1
# e.g.
#
# >chr16:84987588-84988227 | element 1 | positive  | neural tube[12/12] | hindbrain (rhombencephalon)[12/12] | limb[3/12] | cranial nerve[8/12]
# AACTGAAGGGACCCCGTTAGCATAtaaacaaaaggtggggggtagccccgagcctcttct
# ctgacagccagtggcggcagtgatgaatttgtgaagttatctaattttccactgttttaa
# ttagagacttgggctctgaggcctcgcagctggcttctttgtgctgtattctgttgcctg
# acagag

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

  print "Parsing $file with enhancer parser\n";

  my $feature_internal_id = ($self->find_max_id($db_adaptor, "regulatory_feature")) + 1;
  my $highest_factor_id = ($self->find_max_id($db_adaptor, "regulatory_factor")) + 1;

  my $analysis_adaptor = $db_adaptor->get_AnalysisAdaptor();
  my $slice_adaptor = $db_adaptor->get_SliceAdaptor();

  my @features;
  my %feature_objects;
  my %anal;

  # this object is only used for projection
  my $dummy_analysis = new Bio::EnsEMBL::Analysis(-logic_name => 'EnhancerProjection');


  # use separate analyses to distinguish positive and negative influence
  my $analysis_positive = $analysis_adaptor->fetch_by_logic_name("enhancer_positive")->dbID();
  my $analysis_negative = $analysis_adaptor->fetch_by_logic_name("enhancer_negative")->dbID();

  if (!$analysis_positive) {
    print STDERR "Can't get analysis for enhancer_positive, skipping\n";
    exit(1);
  }

  if (!$analysis_negative) {
    print STDERR "Can't get analysis for enhancer_negative, skipping\n";
    exit(1);
  }


  # link to a blank factor so the SQL still works
  my $blank_factor_id = $self->get_blank_factor_id($db_adaptor);

  # read file

  open (FILE, "<$file") || die "Can't open $file";

  while (<FILE>) {

    next if ($_ !~ /^>/); # only read headers

    my %feature;

    # >chr16:84987588-84988227 | element 1 | positive  | neural tube[12/12] | hindbrain (rhombencephalon)[12/12] | limb[3/12] | cranial nerve[8/12]
    my ($coords, $element, $posneg, @stuff) = split /\s+\|\s+/;;

    # parse co-ordinates
    my ($chr, $start, $end) = $coords =~ /chr([^:]+):(\d+)-(\d+)/;

    # parse element name
    my ($element_number) = $element =~ /\s*element\s*(\d+)/;

    # ----------------------------------------
    # Feature name

    $feature{NAME} = "LBNL-$element_number";

    # ----------------------------------------
    # Analysis

    $feature{ANALYSIS_ID} = $posneg eq 'positive' ? $analysis_positive : $analysis_negative;

    # ----------------------------------------
    # Seq_region ID and co-ordinates

    my $chr_slice;

    if ($old_assembly) {
      $chr_slice = $slice_adaptor->fetch_by_region('chromosome', $chr, undef, undef, undef, $old_assembly);
    } else {
      $chr_slice = $slice_adaptor->fetch_by_region('chromosome', $chr);
    }

    if (!$chr_slice) {
      print STDERR "Can't get slice for chromosme $chr\n";
      next;
    }

    my $seq_region_id = $slice_adaptor->get_seq_region_id($chr_slice);

    if (!$seq_region_id) {
      print STDERR "Can't get seq_region_id for chromosome $chr\n";
      next;
    }

    $feature{SEQ_REGION_ID} = $seq_region_id;

    # Assume these are all on the positive strand
    my $strand = 1;

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
    # Feature internal ID
    # note this is not the "id" referred to above

    $feature{INTERNAL_ID} = $feature_internal_id++;

    # ----------------------------------------
    # Evidence

    $feature{EVIDENCE} = "";

    # ----------------------------------------
    # Factor - blank so that the SQL for retrieving still works

    $feature{FACTOR_ID} = $blank_factor_id;

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

  print "Parsed " . scalar(@{$result{FEATURES}}) . " features\n";

  return \%result;

}


sub new {

  my $self = {};
  bless $self, "RegulatoryFeatureParser::enhancer";
  return $self;

}

1;
