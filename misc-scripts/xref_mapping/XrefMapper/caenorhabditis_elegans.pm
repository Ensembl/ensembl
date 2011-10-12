package XrefMapper::caenorhabditis_elegans;

use  XrefMapper::BasicMapper;

use vars qw(@ISA);

@ISA = qw(XrefMapper::BasicMapper);


sub get_set_lists {

  return [["ExonerateGappedBest1", ["caenorhabditis_elegans","*"]]];

}



# Elegans is imported from WormBase. The gene and transcript stable IDs
# are the WormBase identifiers. The display_xref_ids for genes and
# transcripts are calculated directly rather than via the more complex
# priority-based method in BasicMapper.pm

sub build_display_xrefs {

  my ($self, $type, $external_db) = @_;

  print "Setting $type display_xrefs from $type stable IDs\n";
  my $dir = $self->core()->dir();

  my $sql = "UPDATE $type t, xref x, external_db e SET t.display_xref_id=x.xref_id WHERE t.stable_id=x.dbprimary_acc AND e.external_db_id=x.external_db_id AND e.db_name=\'${external_db}\'\n";

  open (SQL, ">$dir/${type}_display_xref.sql");

  print SQL $sql;

  close(SQL);

}


sub build_transcript_display_xrefs {

  my ($self) = @_;

  $self->build_display_xrefs("transcript", "wormbase_transcript");

}

sub build_gene_display_xrefs {

  my ($self) = @_;

  $self->build_display_xrefs("gene", "wormbase_gene");

}


sub gene_description_filter_regexps {

  return ('[0-9A-Z]+\.\d*[A-Z]* PROTEIN[ \.]',
	  '\(\d[A-Z]\d+\)\.',
	  '\([0-9A-Z]+\.\d*[A-Z]* PROTEIN\)[ \.]',
	  '^\(*HYPOTHETICAL\s+.*',
	  '^\s*\(FRAGMENT\)\.?\s*$' );

}

1;
