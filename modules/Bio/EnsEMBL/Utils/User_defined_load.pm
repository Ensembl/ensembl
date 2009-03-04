=head1 LICENSE

  Copyright (c) 1999-2009 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

use strict;
use warnings;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(warning throw);


=head2 Example of user defined loading

Add new databases by creating the load method and use it by setting
ENSEMBL_REGISTRY to your configuration file.

sub load_dupgene {
  my ( $class, @args ) = @_;

  my ($species) = rearrange( [qw(SPECIES)], @args );
  my $group = 'dupgene';

  push( @args, '-group' );
  push( @args, $group );

  my ( $spec, $gro ) =
    Bio::EnsEMBL::Registry->check_if_already_there(@args);
  if ($spec) {
    print STDERR "already defined returning\n";
    return;
  }

  my $dbc = new Bio::EnsEMBL::DBSQL::DBConnection(@args);

  my $dba = new_fast Bio::EnsEMBL::DBSQL::DBAdaptor( '-con' => $dbc );

  Bio::EnsEMBL::Registry->add_DBAdaptor( $species, $group, $dba );

  my %pairs = (
    'Analysis'        => 'Bio::EnsEMBL::DBSQL::AnalysisAdaptor',
    'ArchiveStableId' => 'Bio::EnsEMBL::DBSQL::ArchiveStableIdAdaptor',
    'Attribute'       => 'Bio::EnsEMBL::DBSQL::AttributeAdaptor',
    'AssemblyExceptionFeature' =>
      'Bio::EnsEMBL::DBSQL::AssemblyExceptionFeatureAdaptor',
    'AssemblyMapper' => 'Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor',
    'Blast'          => 'Bio::EnsEMBL::External::BlastAdaptor',
    'MetaContainer'  => 'Bio::EnsEMBL::DBSQL::MetaContainer',
    'CoordSystem'    => 'Bio::EnsEMBL::DBSQL::CoordSystemAdaptor',
    'CompressedSequence' =>
      'Bio::EnsEMBL::DBSQL::CompressedSequenceAdaptor',
    'DBEntry'         => 'Bio::EnsEMBL::DBSQL::DBEntryAdaptor',
    'DnaAlignFeature' => 'Bio::EnsEMBL::DBSQL::DnaAlignFeatureAdaptor',
    'DensityFeature'  => 'Bio::EnsEMBL::DBSQL::DensityFeatureAdaptor',
    'DensityType'     => 'Bio::EnsEMBL::DBSQL::DensityTypeAdaptor',
    'Exon'            => 'Bio::EnsEMBL::DBSQL::ExonAdaptor',
    'Gene'            => 'Bio::EnsEMBL::DBSQL::GeneAdaptor',
    'KaryotypeBand'   => 'Bio::EnsEMBL::DBSQL::KaryotypeBandAdaptor',
    'Marker'          => 'Bio::EnsEMBL::Map::DBSQL::MarkerAdaptor',
    'MarkerFeature' => 'Bio::EnsEMBL::Map::DBSQL::MarkerFeatureAdaptor',
    'MetaCoordContainer' => 'Bio::EnsEMBL::DBSQL::MetaCoordContainer',
    'MiscSet'            => 'Bio::EnsEMBL::DBSQL::MiscSetAdaptor',
    'MiscFeature'        => 'Bio::EnsEMBL::DBSQL::MiscFeatureAdaptor',
    'PredictionTranscript' =>
      'Bio::EnsEMBL::DBSQL::PredictionTranscriptAdaptor',
    'PredictionExon' => 'Bio::EnsEMBL::DBSQL::PredictionExonAdaptor',
    'ProteinFeature' => 'Bio::EnsEMBL::DBSQL::ProteinFeatureAdaptor',
    'ProteinAlignFeature' =>
      'Bio::EnsEMBL::DBSQL::ProteinAlignFeatureAdaptor',
    'ProxySNP'        => 'Bio::EnsEMBL::DBSQL::ProxySNPAdaptor',
    'QtlFeature'      => 'Bio::EnsEMBL::Map::DBSQL::QtlFeatureAdaptor',
    'Qtl'             => 'Bio::EnsEMBL::Map::DBSQL::QtlAdaptor',
    'RepeatConsensus' => 'Bio::EnsEMBL::DBSQL::RepeatConsensusAdaptor',
    'RepeatFeature'   => 'Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptor',
    'Sequence'        => 'Bio::EnsEMBL::DBSQL::SequenceAdaptor',
    'SimpleFeature'   => 'Bio::EnsEMBL::DBSQL::SimpleFeatureAdaptor',
    'Slice'           => 'Bio::EnsEMBL::DBSQL::SliceAdaptor',
    'SupportingFeature' =>
      'Bio::EnsEMBL::DBSQL::SupportingFeatureAdaptor',
    'Transcript'  => 'Bio::EnsEMBL::DBSQL::TranscriptAdaptor',
    'Translation' => 'Bio::EnsEMBL::DBSQL::TranslationAdaptor'
  );

  foreach my $key ( keys %pairs ) {
    Bio::EnsEMBL::Registry->add_adaptor( $species, $group, $key,
      $pairs{$key} );
  }

  #if dnadb has been set then for the follwing use it.
  foreach my $type (
    qw(Sequence AssemblyMapper KaryotypeBand RepeatFeature CoordSystem AssemblyExceptionFeature)
    )
  {
    Bio::EnsEMBL::Registry->set_get_via_dnadb_if_set( $species, $type );
  }
} ## end sub load_dupgene

=cut
