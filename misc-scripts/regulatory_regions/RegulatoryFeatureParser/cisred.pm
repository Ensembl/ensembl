package RegulatoryFeatureParser::cisred;

use strict;

use File::Basename;

# To get files for CisRed data, download the following 2 files (e.g. via wget):
#
# http://www.cisred.org/content/databases_methods/human_2/data_files/motifs.txt
#
# http://www.cisred.org/content/databases_methods/human_2/data_files/search_regions.txt

# Format of motifs.txt (note group_name often blank)

# name	chromosome	start	end	strand	group_name
# craHsap1	X	138337029	138337034	1
# craHsap2	X	138338145	138338150	1
# craHsap3	X	138338363	138338368	1
# craHsap4	X	138338388	138338395	1

# Format of search_regions.txt
# name	chromosome	start	end	strand	ensembl_gene_id
# 1	17	39822200	39824467	-1	ENSG00000005961
# 8	17	23151483	23153621	-1	ENSG00000007171
# 14	1	166434638	166437230	-1	ENSG00000007908
# 19	1	23602820	23605631	-1	ENSG00000007968


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

  my $analysis_adaptor = $db_adaptor->get_AnalysisAdaptor();
  my $slice_adaptor = $db_adaptor->get_SliceAdaptor();

  my $stable_id_to_internal_id = $self->build_stable_id_cache($db_adaptor);

  print "Parsing $file with cisred parser\n";

  # ----------------------------------------
  # We need a "blank" factor for those features which aren't assigned factors
  # Done this way to maintain referential integrity

  my $blank_factor_id = $self->get_blank_factor_id($db_adaptor);

  # ----------------------------------------

  my $feature_internal_id = ($self->find_max_id($db_adaptor, "regulatory_feature")) + 1;
  my $highest_factor_id = ($self->find_max_id($db_adaptor, "regulatory_factor")) + 1;

  my @features;
  my @factors;
  my %factor_ids_by_name; # name -> factor_id
  my %feature_objects;

  # ----------------------------------------
  # Analysis - need one for each type of feature
  my %analysis;

  foreach my $anal ("cisRed", "cisred_search") {  # TODO - add other types as necessary

    my $analysis_obj = $analysis_adaptor->fetch_by_logic_name($anal);

    die "Can't get analysis for $anal, skipping" if (!$analysis_obj);

    $analysis{$anal} = $analysis_obj->dbID();
    print "Analysis ID for $anal is " . $analysis{$anal} . "\n";

  }

  # this object is only used for projection
  my $dummy_analysis = new Bio::EnsEMBL::Analysis(-logic_name => 'CisRedProjection');

  # ----------------------------------------
  # Parse motifs.txt file

  print "Parsing features from $file\n";

  my $skipped = 0;

  my $coords_changed = 0;

  open (FILE, "<$file") || die "Can't open $file";
  <FILE>; # skip header
  while (<FILE>) {

    next if ($_ =~ /^\s*\#/ || $_ =~ /^\s*$/);

    my %feature;
    # name	chromosome	start	end	strand	group_name   ensembl_gene_id
    my ($motif_name, $chromosome, $start, $end, $strand, $group_name, $gene_id) = split (/\t/);
    ($gene_id) = $gene_id =~ /(ENS.*G\d{11})/;

    # ----------------------------------------
    # Feature name & analysis

    $feature{NAME} = $motif_name;
    $feature{INFLUENCE} = "unknown"; # TODO - what does cisRed store?
    $feature{ANALYSIS_ID} = $analysis{cisRed};

    # ----------------------------------------
    # Factor

    # If $group_id is present in %group_sizes we want to create or reuse a factor.
    # If not, this feature is not associated with any factor.
    # If the factor is to be created, its name is crtHsapXX where XX is group_id.
    # TODO - other species prefixes

    $feature{FACTOR_ID} = $blank_factor_id;
    if ($group_name && $group_name ne '' && $group_name !~ /\s/) {
      my $factor_id = $factor_ids_by_name{$group_name};

      if (!$factor_id) { # create one
	my %factor;
	$factor_id = $highest_factor_id + 1;
	$factor{INTERNAL_ID} = $factor_id;
	$factor{NAME} = $group_name;
	$factor{TYPE} = 'NULL';
	push @factors, \%factor;
	$factor_ids_by_name{$factor{NAME}} = $factor_id;
	$highest_factor_id = $factor_id;
	#print join("  ", ("Factor: ", $factor{ID}, $factor{NAME}, $factor{TYPE})) . "\n";
      }
      $feature{FACTOR_ID} = $factor_id;
    }

    # ----------------------------------------
    # Seq_region ID and co-ordinates, projected if necessary

    my $chr_slice = $slice_adaptor->fetch_by_region('chromosome', $chromosome, undef, undef, undef, $old_assembly);

    if (!$chr_slice) {
      print STDERR "Can't get slice for $chromosome:$start:$end\n";
      next;
    }

    my $seq_region_id = $slice_adaptor->get_seq_region_id($chr_slice);

    if (!$seq_region_id) {
      print STDERR "Can't get seq_region_id for chromosome $chromosome\n";
      next;
    }

    $feature{SEQ_REGION_ID} = $seq_region_id;

    # project if necessary
    if ($new_assembly) {

      #print join("\t", "OLD: ", $start, $end, $strand, $chromosome, $motif_name) . "\n";

      my $projected_feature = $self->project_feature($start, $end, $strand, $chromosome, $chr_slice, $dummy_analysis, $new_assembly, $slice_adaptor, $motif_name);

      $coords_changed++ if ($projected_feature->start() != $start || $projected_feature->end() != $end);

      $start = $projected_feature->start();
      $end = $projected_feature->end();
      $strand = $projected_feature->strand();

      $feature{SEQ_REGION_ID} = $slice_adaptor->get_seq_region_id($projected_feature->feature_Slice);

      #print join("\t", "NEW: ", $start, $end, $strand, $chromosome, $motif_name) . "\n";

    }

    $feature{START} = $start;
    $feature{END} = $end;
    $feature{STRAND} = $strand;


    # ----------------------------------------
    # Ensembl object - always a gene in cisRed

    my $ensembl_id = $stable_id_to_internal_id->{gene}->{$gene_id};

    if (!$ensembl_id) {
      print STDERR "Can't get ensembl internal ID for $gene_id, skipping\n";
      print join("-", $motif_name, $chromosome, $start, $end, $strand, $group_name, $gene_id, "\n");
      $skipped++;
      next;
    }

    $feature{ENSEMBL_TYPE} = "Gene";
    $feature{ENSEMBL_ID} = $ensembl_id;

    # ----------------------------------------
    # Feature internal ID

    $feature{INTERNAL_ID} = $feature_internal_id++;

    # ----------------------------------------
    # Evidence

    $feature{EVIDENCE} = "";

    # ----------------------------------------
    # Add to object to be returned

    push @features, \%feature;

    #print "Feature: ";
    #foreach my $field (keys %feature) {
    #	print $field . ": " . $feature{$field} . " ";
    #}
    #print "\n";

  }

  close FILE;

  # ----------------------------------------
  # Search regions 
  # read search_regions.txt from same location as $file
  my $search_regions_file = dirname($file) . "/search_regions.txt";

  my $skipped_sr = 0;

  my @search_regions;
  print "Parsing search regions from $search_regions_file\n";
  open (SEARCH_REGIONS, "<$search_regions_file") || die "Can't open $search_regions_file";
  <SEARCH_REGIONS>; # skip header
  while (<SEARCH_REGIONS>) {
    chomp;
    my ($id, $chromosome, $start, $end, $strand, $ensembl_gene_id) = split;
    my $gene_id = $stable_id_to_internal_id->{gene}->{$ensembl_gene_id};
    if (!$gene_id) {
      warn("Can't get internal ID for $ensembl_gene_id\n");
      $skipped_sr++;
      next;
    }
    my $sr_chr_slice = $slice_adaptor->fetch_by_region('chromosome', $chromosome, undef, undef, undef, $old_assembly);
    if (!$sr_chr_slice) {
      print STDERR "Can't get slice for $chromosome:$start:$end\n";
      next;
    }
    my $sr_seq_region_id = $slice_adaptor->get_seq_region_id($sr_chr_slice);
    if (!$sr_seq_region_id) {
      print STDERR "Can't get seq_region_id for chromosome $chromosome\n";
      next;
    }

    my $name = "CisRed_Search_$id";

    # project if necessary
    if ($new_assembly) {

      #print join("\t", "OLD: ", $start, $end, $strand, $chromosome, $name) . "\n";

      my $projected_region = $self->project_feature($start, $end, $strand, $chromosome, $sr_chr_slice, $dummy_analysis, $new_assembly, $slice_adaptor, "CisRed_Search_$id");

      $start = $projected_region->start();
      $end = $projected_region->end();
      $strand = $projected_region->strand();

      #print join("\t", "NEW: ", $start, $end, $strand, $chromosome, $name) . "\n";

    }

    my %search_region;
    $search_region{NAME} = $name;
    $search_region{SEQ_REGION_ID} = $sr_seq_region_id;
    $search_region{START} = $start;
    $search_region{END} = $end;
    $search_region{STRAND} = ($strand =~ /\+/ ? 1 : -1);
    $search_region{ENSEMBL_OBJECT_TYPE} = 'Gene';
    $search_region{ENSEMBL_OBJECT_ID} = $gene_id;
    $search_region{ANALYSIS_ID} = $analysis{cisred_search};
    push @search_regions, \%search_region;

  }
  close(SEARCH_REGIONS);


  # ----------------------------------------

  $result{FEATURES} = \@features;
  $result{FACTORS} = \@factors;
  $result{SEARCH_REGIONS} = \@search_regions;

  print "Parsed " . scalar(@{$result{FEATURES}}) . " features, " . scalar(@{$result{FACTORS}}) . " factors and " . scalar(@{$result{SEARCH_REGIONS}}) . " search regions\n";

  print "Skipped $skipped features and $skipped_sr search regions due to missing stable ID - internal ID mappings\n";

  print "$coords_changed features had their co-ordinates changed as a result of assembly mapping.\n" if ($new_assembly);

  return \%result;

}

sub new {

  my $self = {};
  bless $self, "RegulatoryFeatureParser::cisred";
  return $self;

}

1;
