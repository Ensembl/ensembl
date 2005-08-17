package RegulatoryFeatureParser::BaseParser;

use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;


# Base functionality for regualatory feature parsers


# --------------------------------------------------------------------------------
# Delete existing regulatory features etc

sub delete_existing {

  my ($db_adaptor, $type) = @_;

  my $t = lc($type);

  print "Deleting existing features & related data for type $type\n";

  # Delete any regulatory_factor_transcript entries first
  my $sth = $db_adaptor->dbc->prepare("DELETE rft FROM regulatory_feature rfeat, regulatory_factor_transcript rft, analysis a WHERE rfeat.regulatory_factor_id=rft.regulatory_factor_id AND a.analysis_id=rfeat.analysis_id AND LOWER(a.analysis_id)=?");
  $sth->execute($t);

  # now delete interlinked regulatory_feature, regulatory_factor and regulatory_feature_object entries
  $sth = $db_adaptor->dbc->prepare("DELETE rfeat, rfact, rfo FROM regulatory_feature rfeat, regulatory_factor rfact, regulatory_feature_object rfo, analysis a WHERE rfeat.regulatory_feature_id=rfo.regulatory_feature_id AND rfeat.regulatory_factor_id=rfact.regulatory_factor_id AND rfeat.analysis_id=a.analysis_id AND LOWER(a.logic_name)=?");
  $sth->execute($t);

  # delete dangling regulatory_factors
  $sth = $db_adaptor->dbc->prepare("DELETE rfact FROM regulatory_feature rfeat, regulatory_factor rfact, analysis a WHERE rfeat.regulatory_factor_id=rfact.regulatory_factor_id AND a.analysis_id=rfeat.analysis_id AND LOWER(a.analysis_id)=?");
  $sth->execute($t);

  # and finally any dangling regulatory features
  $sth = $db_adaptor->dbc->prepare("DELETE rfeat FROM regulatory_feature rfeat, analysis a WHERE a.analysis_id=rfeat.analysis_id AND LOWER(a.analysis_id)=?");
  $sth->execute($t);

}

# --------------------------------------------------------------------------------
# Check that the type specified corresponds to an analysis

sub validate_type {

  my ($db_adaptor, $type) = @_;

  my $sth = $db_adaptor->dbc->prepare("SELECT analysis_id FROM analysis WHERE LOWER(logic_name)=?");

  $sth->execute(lc($type));
  if ($sth->fetchrow_array()) {
    print "Type $type is valid\n";
  } else {
    print "Type $type is not valid\n";
    return 0;
  }

  return 1;

}


# --------------------------------------------------------------------------------
# Find the maximum existing ID in a table.
sub find_max_id {

  my ($self, $db_adaptor, $table) = @_;

  my $row = @{$db_adaptor->dbc->db_handle->selectall_arrayref("SELECT MAX(${table}_id) FROM $table")}[0];
  my $max_id = @{$row}[0];
  if (!defined $max_id) {
    print "No existing ${table}_ids, will start from 1\n";
    $max_id = 0;
  } else {
    print "Maximum existing ${table}_id = $max_id\n";
  }

  return $max_id;

}

# --------------------------------------------------------------------------------
# Build a cache of ensembl stable ID -> internal ID
# Return hashref keyed on {$type}{$stable_id}
# Type is always all lower case

sub build_stable_id_cache {

  my ($self, $db_adaptor) = @_;

  my %stable_id_to_internal_id;

  foreach my $type ('gene', 'transcript', 'translation') { # Add exon here if required

    print "Caching stable ID -> internal ID links for ${type}s\n";

    my $sth = $db_adaptor->dbc->prepare("SELECT ${type}_id, stable_id FROM ${type}_stable_id");
    $sth->execute();
    my ($internal_id, $stable_id);
    $sth->bind_columns(\$internal_id, \$stable_id);

    while ($sth->fetch) {

      $stable_id_to_internal_id{$type}{$stable_id} = $internal_id;

    }

  }

  return \%stable_id_to_internal_id;

}

# --------------------------------------------------------------------------------
# Upload features and factors etc to database

sub upload_features_and_factors {

  my ($self, $db_adaptor, $objects) = @_;

  my $dbc = $db_adaptor->dbc;

  my $feature_sth = $dbc->prepare("INSERT INTO regulatory_feature (regulatory_feature_id, name, seq_region_id, seq_region_start, seq_region_end, seq_region_strand, analysis_id, regulatory_factor_id) VALUES(?,?,?,?,?,?,?,?)");
  my $factor_sth = $dbc->prepare("INSERT INTO regulatory_factor (regulatory_factor_id, name, type) VALUES(?,?,?)");
  my $feature_object_sth = $dbc->prepare("INSERT INTO regulatory_feature_object (regulatory_feature_id, ensembl_object_type, ensembl_object_id, influence, evidence) VALUES(?,?,?,?,?)");

  print "Uploading " . scalar(@{$objects->{FEATURES}}) . " features ...\n";

  foreach my $feature (@{$objects->{FEATURES}}) {

    $feature_sth->execute($feature->{INTERNAL_ID},
			  $feature->{NAME},
			  $feature->{SEQ_REGION_ID},
			  $feature->{START},
			  $feature->{END},
			  $feature->{STRAND},
			  $feature->{ANALYSIS_ID},
			  $feature->{FACTOR_ID});

    $feature_object_sth->execute($feature->{INTERNAL_ID},
				  $feature->{ENSEMBL_TYPE},
				  $feature->{ENSEMBL_ID},
				  $feature->{INFLUENCE},
				  $feature->{EVIDENCE});

  }

  print "Uploading " . scalar(@{$objects->{FACTORS}}) . " factors ...\n";

  foreach my $factor (@{$objects->{FACTORS}}) {

    $factor_sth->execute($factor->{INTERNAL_ID},
			 $factor->{NAME},
			 $factor->{TYPE});

  }

  print "Done\n";

}

# --------------------------------------------------------------------------------

1;
