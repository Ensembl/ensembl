=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

package XrefParser::WormbaseDatabaseStableIDParser;

# Read gene and transcript stable IDs from a WormBase-imported database and create
# xrefs and direct xrefs for them.
# Note the only things that are actually WormBase specific here are the source names
# for the wormbase_gene and wormbase_transcript sources.

use strict;
use warnings;
use Carp;
use base qw( XrefParser::DatabaseParser );

sub run {
  my ($self, $ref_arg) = @_;

  my $source_id    = $ref_arg->{source_id} || confess "Need a source_id";
  my $species_id   = $ref_arg->{species_id};
  my $dsn          = $ref_arg->{dsn};
  my $verbose      = $ref_arg->{verbose};

  if((!defined $source_id) or (!defined $species_id) or (!defined $dsn) ){
    croak "Need to pass source_id, species_id and dsn as pairs";
  }
  $verbose |=0;

  my $db = $self->connect($dsn); # source db (probably core)
  my $xref_db = $self->dbi();

  my $xref_sth = $xref_db->prepare( "INSERT INTO xref (accession,label,source_id,species_id) VALUES (?,?,?,?)" );


  # read stable IDs
  foreach my $type ('gene', 'transcript') {

    my $direct_xref_sth = $xref_db->prepare( "INSERT INTO ${type}_direct_xref (general_xref_id,ensembl_stable_id,linkage_xref) VALUES (?,?,?)" );
    if($verbose) { 
      print "Building xrefs from $type stable IDs\n";
    }

    my $wb_source_id = $self->get_source_id_for_source_name("wormbase_$type");

    my $sth = $db->prepare( "SELECT stable_id FROM ${type}" );
    $sth->execute();

    while(my @row = $sth->fetchrow_array()) {

      my $id = $row[0];
      # add an xref & a direct xref
      $xref_sth->execute($id, $id, $wb_source_id, $species_id);
      my $xref_id = $xref_sth->{'mysql_insertid'};
      $direct_xref_sth->execute($xref_id, $id, 'Stable ID direct xref');

    } # while fetch stable ID

  } # foreach type
  return 0;
}

1;

