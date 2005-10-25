package RegulatoryFeatureParser::cisred;

use strict;

use File::Basename;

# To get files for CisRed data, connect to db.cisred.org as anonymous
# and use the queries below.
#
# Note that you can't connect directly to it from the Sanger machines
# (firewall issues) so I do it from the EBI and then scp the data across
# to the Sanger machine.										
#
#   mysql -u anonymous -h db.cisred.org -e 'select f.id as motif_id, f.seqname as chromosome, f.start, f.end, f.ensembl_gene_id,g.group_id from features f, group_content g where f.id=g.feature_id' cisred_Hsap_1_2e > motifs.txt
#
#   mysql -u anonymous -h db.cisred.org -e 'select distinct(group_id), count(*) from group_content where group_id !=0 and group_id != -1 group by group_id having count(*) > 1' cisred_Hsap_1_2e > group_sizes.txt
#
# The queries should take ~ 1 minute and ~ 2 seconds respectively.
#
# For the second query, if a group_id has an entry in this file then the regulatory_factor should be assigned as "crtHsapXX" where XX is the group_id. If there's no entry, don't assign a regulatory_factor.
#
# Note all features are unstranded

# Format of motifs.txt
# motif_id	chromosome	start	        end	        ensembl_gene_id	  group_id
# 6	        10	        43087959	43087968	ENSG00000198915	  -1
# 15	        10	        43087044	43087053	ENSG00000198915	  -1
# 20	        10	        43087045	43087054	ENSG00000198915	  -1
# 37	        10	        43086873	43086884	ENSG00000198915	  -1
# 51	        5	        42855740	42855752	ENSG00000198865	  -1
# 54	        10	        43082459	43082470	ENSG00000198915	  -1

# Format of group_sizes.txt
# group_id	count(*)
# 1	        8
# 2	        7
# 3	        4


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

  # ----------------------------------------
  # We need a "blank" factor for those features which aren't assigned factors
  # Done this way to maintain referential integrity

  my $blank_factor_id = $self->get_blank_factor_id($db_adaptor);

  # ----------------------------------------

  my $feature_internal_id = ($self->find_max_id($db_adaptor, "regulatory_feature")) + 1;
  my $highest_factor_id = ($self->find_max_id($db_adaptor, "regulatory_factor")) + 1;

  my $analysis_adaptor = $db_adaptor->get_AnalysisAdaptor();
  my $slice_adaptor = $db_adaptor->get_SliceAdaptor();

  my @features;
  my @factors;
  my %factor_ids_by_name; # name -> factor_id
  my %feature_objects;

  # TODO - regulatory_factor_coding

  my $stable_id_to_internal_id = $self->build_stable_id_cache($db_adaptor);

  # read group_sizes.txt from same location as $file
  my $group_sizes_file = dirname($file) . "/group_sizes.txt";

  my %group_sizes;
  print "Parsing group sizes from $group_sizes_file\n";
  open (GROUP_SIZES, "<$group_sizes_file") || die "Can't open $group_sizes_file";
  <GROUP_SIZES>; # skip header
  while (<GROUP_SIZES>) {
    my ($file_group_id, $size) = split;
    $group_sizes{$file_group_id} = $size;
  }
  close(GROUP_SIZES);

  # ----------------------------------------
  # Analysis

  my $analysis = $analysis_adaptor->fetch_by_logic_name("cisRed");

  if (!$analysis) {
    print STDERR "Can't get analysis for cisRed, skipping\n";
    next;
  }

  my $analysis_id = $analysis->dbID();

  # ----------------------------------------
  # Parse motifs.txt file

  print "Parsing features from $file\n";

  open (FILE, "<$file") || die "Can't open $file";
  <FILE>; # skip header
  while (<FILE>) {

    next if ($_ =~ /^\s*\#/ || $_ =~ /^\s*$/);

    my %feature;

    my ($motif_id, $chromosome, $start, $end, $gene_id, $group_id) = split;

    # ----------------------------------------
    # Feature name
    # Name is craHsap + motif_id 
    # TODO - other species

    $feature{NAME} = "craHsap" . $motif_id;
    $feature{INFLUENCE} = "unknown"; # TODO - what does cisRed store?
    $feature{ANALYSIS_ID} = $analysis_id;

    # ----------------------------------------
    # Factor

    # If $group_id is present in %group_sizes we want to create or reuse a factor.
    # If not, this feature is not associated with any factor.
    # If the factor is to be created, its name is crtHsapXX where XX is group_id.
    # TODO - other species prefixes

    $feature{FACTOR_ID} = $blank_factor_id;
    if ($group_sizes{$group_id}) {
      my $factor_id = $factor_ids_by_name{$feature{NAME}};

      if (!$factor_id) { # create one
	my %factor;
	$factor_id = $highest_factor_id + 1;
	$factor{INTERNAL_ID} = $factor_id;
	$factor{NAME} = "crtHsap" . $group_id;
	$factor{TYPE} = $factor{NAME}; # TODO - error checking that type is one of the enums?
	push @factors, \%factor;
	$factor_ids_by_name{$factor{NAME}} = $factor_id;
	$highest_factor_id = $factor_id;
	#print join("  ", ("Factor: ", $factor{ID}, $factor{NAME}, $factor{TYPE})) . "\n";
      }
      $feature{FACTOR_ID} = $factor_id;
    }

    # ----------------------------------------
    # Seq_region ID and co-ordinates

    my $chr_slice = $slice_adaptor->fetch_by_region(undef, $chromosome, $start, $end);

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
    $feature{START} = $chr_slice->start();
    $feature{END} = $chr_slice->end();

    # Note CisRed should be unstranded
    $feature{STRAND} = $chr_slice->strand();

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

  $result{FEATURES} = \@features;
  $result{FACTORS} = \@factors;

  print "Parsed " . scalar(@{$result{FEATURES}}) . " features and " . scalar(@{$result{FACTORS}}) . " factors\n";

  return \%result;

}


sub get_blank_factor_id () {

  my ($self, $db_adaptor) = @_;

  my $sth = $db_adaptor->dbc->prepare("SELECT regulatory_factor_id FROM regulatory_factor WHERE name=''");
  $sth->execute();

  my ($factor_id) = $sth->fetchrow_array();

  if ($factor_id) {
    print "Found existing blank factor, id=$factor_id\n";
  } else {
     $db_adaptor->dbc->do("INSERT INTO regulatory_factor (name) VALUES ('')");
     $sth->execute();
     ($factor_id) = $sth->fetchrow_array();
     print "Created new blank factor, id=$factor_id\n";
  }

  return $factor_id;

}


sub new {

  my $self = {};
  bless $self, "RegulatoryFeatureParser::cisred";
  return $self;

}

1;
