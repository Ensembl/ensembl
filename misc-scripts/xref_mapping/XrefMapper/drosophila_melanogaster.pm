package XrefMapper::drosophila_melanogaster;

use  XrefMapper::BasicMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };

sub get_set_lists {

  return [["ExonerateGappedBest1", ["drosophila_melanogaster","*"]]];

}

sub gene_description_filter_regexps {

  return ();

}

# Special logic for drosophila display_xrefs:
#
# gene: flybase_name if present, else gadfly_gene_cgid
#
# transcript: flybase_name if present, else gadfly_transcript_cgid

sub build_transcript_display_xrefs {

  my ($self) = @_;

  $self->build_xref_to_source_mappings();

  $self->build_display_xrefs("Transcript", "flybase_name", "gadfly_transcript_cgid");

}

sub build_gene_display_xrefs {

  my ($self) = @_;

  $self->build_display_xrefs("Gene", "flybase_name", "gadfly_gene_cgid");

}

sub build_display_xrefs {

  my ($self, $type, $first_source, $second_source) = @_;

  print "Building " . lc($type) . " display_xrefs for drosophila\n";

  my ($self, $xref_id_offset) = @_;

  my $dir = $self->core->dir();

  $self->read_direct_xref_mappings();

  print "Building " . lc($type) . " display_xrefs\n";
  my $n = 0;

  # go through each object/xref mapping and store the best ones as we go along
  my %obj_to_best_xref;

  OBJECT: foreach my $key (keys %object_xref_mappings) {

    my ($obj_type, $object_id) = split /\|/, $key;

    next if ($obj_type !~ /$type/i);

    # if the ensembl_object has more than one associated xref,
    # use the $first_source if present, otherwise the $second_source
    my @xrefs = @{$object_xref_mappings{$key}};
    my ($best_xref);
    foreach my $xref (@xrefs) {

      my $source = $xref_to_source{$xref};
      if ($source) {

	if ($source =~ /$first_source/) {
	  $obj_to_best_xref{$key} = $xref;
	  next OBJECT;
	}

	if ($source =~ /$second_source/) {
	  $obj_to_best_xref{$key} = $best_xref;
          print "Found $second_source for $type $object_id\n";
	}

      } else {
	warn("Couldn't find a source for xref id $xref " . $xref_accessions{$xref_id} . "\n");
      }
    }

  }

  open (DX, ">$dir/" . lc($type) . "_display_xref.sql");
  open (DX_TXT, ">$dir/" . lc($type) . "_display_xref.txt");

  foreach my $key (keys %obj_to_best_xref) {

    my ($obj_type, $object_id) = split /\|/, $key;
    $best_xref = $obj_to_best_xref{$key};
    # Write record with xref_id_offset
    print DX "UPDATE " . lc($type) . " SET display_xref_id=" . ($best_xref+$xref_id_offset) . " WHERE " . lc($type) . "_id=" . $object_id . ";\n";
    print DX_TXT ($best_xref+$xref_id_offset) . "\t" . $object_id . "\n";
    $n++;

  }

  close(DX);
  close(DX_TXT);

  print "Wrote $n " . lc($type) . " display_xref entries to " . lc($type) . "_display_xref.sql\n";

}

# Cache the source of each xref

sub build_xref_to_source_mappings {

  my $self = shift;

  # get a list of xref sources; format:
  # key: xref_id value: source_name
  # lots of these; if memory is a problem, just get the source ID (not the name)
  # and look it up elsewhere
  # note %xref_to_source is global

  print "Building xref->source mapping table\n";
  my $sql = "SELECT x.xref_id, s.name FROM source s, xref x WHERE x.source_id=s.source_id";
  my $sth = $self->xref->dbc->prepare($sql);
  $sth->execute();

  my ($xref_id, $source_name);
  $sth->bind_columns(\$xref_id, \$source_name);

  while ($sth->fetch()) {
    $xref_to_source{$xref_id} = $source_name;
  }

  print "Got " . scalar(keys %xref_to_source) . " xref-source mappings\n";

}


# build mappings from direct xrefs

sub read_direct_xref_mappings {

  my ($self) = @_;

  # will need stable_id->internal_id mapping
  $stable_id_to_internal_id = $self->build_stable_id_to_internal_id_hash() if (!$stable_id_to_internal_id);

  print "Reading direct xref mappings\n";

  # SQL / statement handle for getting all direct xrefs
  my $xref_sql = "SELECT dx.general_xref_id, dx.ensembl_stable_id, dx.type, dx.linkage_xref, x.accession, x.version, x.label, x.description, x.source_id, x.species_id FROM direct_xref dx, xref x WHERE dx.general_xref_id=x.xref_id";
  my $xref_sth = $self->xref->dbc->prepare($xref_sql);

  $xref_sth->execute();

  my ($xref_id, $ensembl_stable_id, $type, $linkage_xref, $accession, $version, $label, $description, $source_id, $species_id);
  $xref_sth->bind_columns(\$xref_id, \$ensembl_stable_id, \$type, \$linkage_xref,\ $accession, \$version, \$label, \$description, \$source_id, \$species_id);

  my $n = 0;

  # Get source-external_db mappings
  %source_to_external_db = $self->map_source_to_external_db();

  while ($xref_sth->fetch()) {

    my $external_db_id = $source_to_external_db{$source_id};
    if ($external_db_id) {

      my $ensembl_internal_id = $stable_id_to_internal_id->{$type}->{$ensembl_stable_id};

      if ($ensembl_internal_id) {

	my $key = $type . "|" . $ensembl_internal_id;
	push @{$object_xref_mappings{$key}}, $xref_id;
	$n++;

      }
    }
  }

  print "Read $n direct_xref mappings\n";

}

1;
