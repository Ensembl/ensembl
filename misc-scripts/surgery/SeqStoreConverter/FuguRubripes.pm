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


use strict;
use warnings;

use SeqStoreConverter::BasicConverter;

package SeqStoreConverter::FuguRubripes;

use vars qw(@ISA);

@ISA = qw(SeqStoreConverter::BasicConverter);


sub create_coord_systems {
  my $self = shift;

  $self->debug("FuguRubripes Specific: creating scaffold coord system");

  #
  # Load the coord system table.
  # Fugu has only one coord system : scaffold
  #
  my $default_assembly = $self->get_default_assembly();
  my $target = $self->target();

  my $sth = $self->dbh()->prepare
   ("INSERT INTO $target.coord_system (name, version, attrib, rank)" .
    "VALUES (?,?,?,?)");
  $sth->execute('scaffold', $default_assembly,
                'default_version,sequence_level', 1);

  $sth->finish();
}



sub create_seq_regions {
  my $self = shift;

  $self->debug("FuguRubripes Specific: creating scaffolds");

  $self->contig_to_seq_region('scaffold');
}


sub create_assembly {
  my $self = shift;

  $self->debug("FuguRubripes Specific: no assembly data needed");
  #fugu has no assembly table
}


sub transfer_genes {
  my $self = shift;

  my $target = $self->target();
  my $source = $self->source();
  my $dbh    = $self->dbh();

  #
  # This is simple in fugu since all genes are actually in contig (now
  # renamed scaffold) coords.  We don't need joins to the assembly table
  # and we don't have to worry about stickies
  #

  $self->debug("FuguRubripes Specific: Building gene table " .
               "(no chromosomal conversion)");

  #
  # Transfer the gene table
  #

  $dbh->do
    ("INSERT INTO $target.gene " .
     "SELECT g.gene_id, g.type, g.analysis_id, e.contig_id, " .
     "       MIN(e.contig_start), MAX(e.contig_end), e.contig_strand, " .
     "       g.display_xref_id " .
     "FROM   $source.transcript t, $source.exon_transcript et, " .
     "       $source.exon e, $source.gene g " .
     "WHERE  t.transcript_id = et.transcript_id " .
     "AND    et.exon_id = e.exon_id " .
     "AND    g.gene_id = t.gene_id " .
     "GROUP BY g.gene_id");

  $self->debug("FuguRubripes Specific: Building transcript table " .
               "(no chromosomal conversion)");

  #
  # Transfer transcript table
  #

  $dbh->do
    ("INSERT INTO $target.transcript " .
     "SELECT t.transcript_id, t.gene_id, e.contig_id, " .
     "       MIN(e.contig_start), MAX(e.contig_end), e.contig_strand, " .
     "       t.display_xref_id " .
     "FROM   $source.transcript t, $source.exon_transcript et, " .
     "       $source.exon e " .
     "WHERE  t.transcript_id = et.transcript_id " .
     "AND    et.exon_id = e.exon_id " .
     "GROUP BY t.transcript_id");

  $self->debug("FuguRubripes Specific: Building exon table " .
               "(no chromosomal conversion)");


  #
  # Transfer exon table
  #

  $self->debug("FuguRubripes Specific: Building transcript table " .
               "(no chromosomal conversion)");

  $dbh->do
    ("INSERT INTO $target.exon " .
     "SELECT e.exon_id, e.contig_id, " .
     "       e.contig_start, e.contig_end, e.contig_strand, " .
     "       e.phase, e.end_phase " .
     "FROM   $source.exon e");


  # 
  # Transfer translation table
  # 

  $self->debug("Building translation table");

  $dbh->do
    ("INSERT INTO $target.translation " .
     "SELECT tl.translation_id, ts.transcript_id, tl.seq_start, " .
     "       tl.start_exon_id, tl.seq_end, tl.end_exon_id " .
     "FROM $source.transcript ts, $source.translation tl " .
     "WHERE ts.translation_id = tl.translation_id");

  return;
}


1;
