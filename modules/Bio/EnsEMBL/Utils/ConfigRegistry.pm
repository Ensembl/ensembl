#
# Ensembl module for Registry
#
# Copyright EMBL/EBI
##
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Utils::ConfigRegistry;

=head1 SYNOPSIS


Bio::EnsEMBL::Utils::ConfigRegistry->load_core( 
                      -species => "Homo Sapiens",
                      -host => 'kaka.sanger.ac.uk',
                      -user => 'anonymous',
                      -dbname => 'homo_sapiens_core_20_34c',
                      -port => '3306' );
 
Bio::EnsEMBL::Utils::ConfigRegistry->load_estgene(
                      -species => "Homo Sapiens",
                      -host => 'kaka.sanger.ac.uk',
                      -user => 'anonymous',
                      -dbname => 'homo_sapiens_estgene_20_34c',
                      -port => '3306' );

=head1 DESCRIPTION

The ConfigRegistry will "Register" a set of adaptors based on the type of
database that is being loaded.

=head1 CONTACT

Post questions to the Ensembl developer list: <ensembl-dev@ebi.ac.uk>


=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::Utils::ConfigRegistry;

#use Exporter;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::ProxySNPAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(warning throw  deprecate stack_trace_dump);
#use vars qw(@ISA @EXPORT_OK);
#@ISA = qw(Exporter);

#@EXPORT_OK = qw(&load_core &load_estgene);



=head2 load_core, load_estgene, load_vega, load_compara, load_pipeline, load_SNP, load_lite

  Arg [DBNAME] : string
                 The name of the database to connect to.
  Arg [HOST] : (optional) string
               The domain name of the database host to connect to.
               'localhost' by default.
  Arg [USER] : string
               The name of the database user to connect with
  Arg [PASS] : (optional) string
               The password to be used to connect to the database
  Arg [PORT] : int
               The port to use when connecting to the database
               3306 by default.
  Arg [DRIVER] : (optional) string
                 The type of database driver to use to connect to the DB
                 mysql by default.
  Arg [DISCONNECT_WHEN_INACTIVE]: (optional) boolean
                 If set to true, the database connection will be disconnected
                 everytime there are no active statement handles. This is
                 useful when running a lot of jobs on a compute farm
                 which would otherwise keep open a lot of connections to the
                 database.  Database connections are automatically reopened
                 when required.
  Arg [SPECIES] : string
                  The name of the species to be used in the registry.

  Example: Bio::EnsEMBL::Utils::ConfigRegistry->load_core(
                                  -species => "Homo Sapiens",
                                  -host    => 'kaka.sanger.ac.uk',
                                  -user    => 'anonymous',
                                  -dbname  => 'homo_sapiens_core_20_34c',
                                  -port    => '3306' );

  Description: Load all the necesary adaptors into the Register for the
               database.
  Returntype : none.
  Exceptions : thrown if USER or DBNAME are not specified, or if the
               database cannot be connected to.

=cut

sub load_core{
  my ($class, @args) = @_;

  my ($species) = rearrange([qw(SPECIES)],@args);

#  print STDERR "species is $species\n";
  my $group = 'core';

  push (@args, '-group');
  push (@args, $group);

  my ($spec,$gro) = Bio::EnsEMBL::Registry->check_if_already_there(@args);
  if($spec){
    return;
  }
  my $dbc = new Bio::EnsEMBL::DBSQL::DBConnection(@args);

  my $dba = new_fast  Bio::EnsEMBL::DBSQL::DBAdaptor('-con' => $dbc);

  Bio::EnsEMBL::Registry->add_DBAdaptor($species, $group, $dba);

  my %pairs =  ( 'Analysis'             => 'Bio::EnsEMBL::DBSQL::AnalysisAdaptor',
		 'ArchiveStableId'      => 'Bio::EnsEMBL::DBSQL::ArchiveStableIdAdaptor',
		 'Attribute'            => 'Bio::EnsEMBL::DBSQL::AttributeAdaptor',
		 'AssemblyExceptionFeature' => 'Bio::EnsEMBL::DBSQL::AssemblyExceptionFeatureAdaptor',
		 'AssemblyMapper'       => 'Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor',
		 'Blast'                => 'Bio::EnsEMBL::External::BlastAdaptor',
		 'MetaContainer'        => 'Bio::EnsEMBL::DBSQL::MetaContainer',
		 'CoordSystem'   => 'Bio::EnsEMBL::DBSQL::CoordSystemAdaptor',
		 'CompressedSequence' => 'Bio::EnsEMBL::DBSQL::CompressedSequenceAdaptor',
		 'DBEntry'              => 'Bio::EnsEMBL::DBSQL::DBEntryAdaptor',
		 'DnaAlignFeature'      => 'Bio::EnsEMBL::DBSQL::DnaAlignFeatureAdaptor',
		 'DensityFeature'       => 'Bio::EnsEMBL::DBSQL::DensityFeatureAdaptor',
		 'DensityType'          => 'Bio::EnsEMBL::DBSQL::DensityTypeAdaptor',
		 'Exon'                 => 'Bio::EnsEMBL::DBSQL::ExonAdaptor',
		 'Gene'                 => 'Bio::EnsEMBL::DBSQL::GeneAdaptor',
		 'KaryotypeBand'        => 'Bio::EnsEMBL::DBSQL::KaryotypeBandAdaptor',
		 'Marker'               => 'Bio::EnsEMBL::Map::DBSQL::MarkerAdaptor',
		 'MarkerFeature'        =>
		 'Bio::EnsEMBL::Map::DBSQL::MarkerFeatureAdaptor',
		 'MetaCoordContainer'   => 'Bio::EnsEMBL::DBSQL::MetaCoordContainer',
		 'MiscSet'              => 'Bio::EnsEMBL::DBSQL::MiscSetAdaptor',
		 'MiscFeature'          => 'Bio::EnsEMBL::DBSQL::MiscFeatureAdaptor',
		 'PredictionTranscript' => 'Bio::EnsEMBL::DBSQL::PredictionTranscriptAdaptor',
		 'PredictionExon'       => 'Bio::EnsEMBL::DBSQL::PredictionExonAdaptor',
		 'ProteinFeature'       => 'Bio::EnsEMBL::DBSQL::ProteinFeatureAdaptor',
		 'ProteinAlignFeature'  =>
		 'Bio::EnsEMBL::DBSQL::ProteinAlignFeatureAdaptor',
		 'ProxySNP'             => 'Bio::EnsEMBL::DBSQL::ProxySNPAdaptor',
		 'QtlFeature'           => 'Bio::EnsEMBL::Map::DBSQL::QtlFeatureAdaptor',
		 'Qtl'                  => 'Bio::EnsEMBL::Map::DBSQL::QtlAdaptor',
		 'RepeatConsensus'      => 'Bio::EnsEMBL::DBSQL::RepeatConsensusAdaptor',
		 'RepeatFeature'        => 'Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptor',
		 'Sequence'             => 'Bio::EnsEMBL::DBSQL::SequenceAdaptor',
		 'SimpleFeature'        => 'Bio::EnsEMBL::DBSQL::SimpleFeatureAdaptor',
		 'Slice'                => 'Bio::EnsEMBL::DBSQL::SliceAdaptor',
		 'SupportingFeature'    =>
		 'Bio::EnsEMBL::DBSQL::SupportingFeatureAdaptor',
		 'Transcript'           => 'Bio::EnsEMBL::DBSQL::TranscriptAdaptor',
		 'Translation'          => 'Bio::EnsEMBL::DBSQL::TranslationAdaptor');

  foreach my $key (keys %pairs){

    Bio::EnsEMBL::Registry->add_adaptor($species, $group, $key, $pairs{$key});
  }


  foreach my $type (qw(Sequence AssemblyMapper KaryotypeBand RepeatFeature CoordSystem AssemblyExceptionFeature)){
    Bio::EnsEMBL::Registry->set_get_via_dnadb_if_set($species,$type);
  }

}


sub load_lite{
  my ($class, @args) = @_;
  require Bio::EnsEMBL::Lite::DBAdaptor;

  my ($species) = rearrange([qw(SPECIES)],@args);
  my $group = 'lite';


  push (@args, '-group');
  push (@args, $group);

  my $dbc = new Bio::EnsEMBL::DBSQL::DBConnection(@args);

  my $dba = new_fast  Bio::EnsEMBL::Lite::DBAdaptor('-con' => $dbc);

  Bio::EnsEMBL::Registry->add_DBAdaptor($species, $group, $dba);

  eval "require Bio::EnsEMBL::Lite::SNPAdaptor";


  my $adap = Bio::EnsEMBL::Lite::SNPAdaptor->new($dba);

  my $prim_adap = Bio::EnsEMBL::DBSQL::ProxySNPAdaptor->new($dba,$adap);

  Bio::EnsEMBL::Registry->add_adaptor($species, $group, $group, $adap);
  Bio::EnsEMBL::Registry->add_adaptor($species, $group, "ProxySNP", $prim_adap);

}

sub load_SNP{
  my ($class, @args) = @_;
  require Bio::EnsEMBL::ExternalData::SNPSQL::DBAdaptor;
  require Bio::EnsEMBL::ExternalData::SNPSQL::SNPAdaptor;

  my ($species) = rearrange([qw(SPECIES)],@args);
  my $group = 'SNP';

  push (@args, '-group');
  push (@args, $group);

  my $dbc = new Bio::EnsEMBL::DBSQL::DBConnection(@args);

  my $dba = new_fast  Bio::EnsEMBL::ExternalData::SNPSQL::DBAdaptor('-con' => $dbc);

  Bio::EnsEMBL::Registry->add_DBAdaptor($species, $group, $dba);

  eval "require Bio::EnsEMBL::ExternalData::SNPSQL::DBAdaptor";
  my $adap = Bio::EnsEMBL::ExternalData::SNPSQL::SNPAdaptor->new($dba);
  my $prim_adap = Bio::EnsEMBL::DBSQL::ProxySNPAdaptor->new($dba,$adap);


  Bio::EnsEMBL::Registry->add_adaptor($species, $group, $group, $adap);

  Bio::EnsEMBL::Registry->add_adaptor($species, $group, "ProxySNP", $prim_adap);
}

sub attach_database{
  my ($class, $species, $core, $name1) = @_;

  my $first =  Bio::EnsEMBL::Registry->get_DBAdaptor($species,$name1);
  my $coredb = Bio::EnsEMBL::Registry->get_DBAdaptor($species,$core);

  Bio::EnsEMBL::Registry->add_db($coredb,$name1, $first);

}

sub attach_dna{
  my ($class, $species, $main, $attach) = @_; 

  my $no_seq =  Bio::EnsEMBL::Registry->get_DBAdaptor($species,$main);
  my $seq = Bio::EnsEMBL::Registry->get_DBAdaptor($species,$attach);

  Bio::EnsEMBL::Registry->add_DNAAdaptor($species,$no_seq->group,$seq);
}


sub load_blast{
  my ($class, @args) = @_;

  my ($species) = rearrange([qw(SPECIES)],@args);
  if(!defined($species)){
    $species = 'NONE'; # NOT species dependent.
    push (@args, '-species');
    push (@args, $species);
  }
  my $group = 'blast';

  push (@args, '-group');
  push (@args, $group);

  eval "require Bio::EnsEMBL::External::BlastAdaptor";

  my $dba = Bio::EnsEMBL::External::BlastAdaptor->new(@args);
  my $dbc = $dba->db();

  my $key = 'Blast';

  Bio::EnsEMBL::Registry->add_adaptor($species, $group, $key, $dba);
}


sub add_blast_link{
  my ($class, $species, $group) = @_;
  my $dba =undef;
  my $blast = undef;

  if($dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species, $group)){
  }
  else{
    throw("Cannot find $species $group in the registry\n");
  }

  if($blast = Bio::EnsEMBL::Registry->get_adaptor("NONE", "blast", "Blast")){
    $dba->add_db_adaptor("Blast",$blast);
  }
  else{
    throw("Sorry no Blast database has been set up to link to.\n");
  }
}



sub load_estgene{
  my ($class, @args) = @_;

  my ($species) = rearrange([qw(SPECIES)],@args);
  my $group = 'estgene';

  push (@args, '-group');
  push (@args, $group);

  my $dbc = new Bio::EnsEMBL::DBSQL::DBConnection(@args);

  my $dba = new_fast  Bio::EnsEMBL::DBSQL::DBAdaptor('-con' => $dbc);

  Bio::EnsEMBL::Registry->add_DBAdaptor($species, $group, $dba);

  my %pairs =  ( 'Analysis'             => 'Bio::EnsEMBL::DBSQL::AnalysisAdaptor',
		 'ArchiveStableId'      => 'Bio::EnsEMBL::DBSQL::ArchiveStableIdAdaptor',
		 'Attribute'            => 'Bio::EnsEMBL::DBSQL::AttributeAdaptor',
		 'AssemblyExceptionFeature' => 'Bio::EnsEMBL::DBSQL::AssemblyExceptionFeatureAdaptor',
		 'AssemblyMapper'       => 'Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor',
		 'Blast'                => 'Bio::EnsEMBL::External::BlastAdaptor',
		 'MetaContainer'        => 'Bio::EnsEMBL::DBSQL::MetaContainer',
		 'CoordSystem'   => 'Bio::EnsEMBL::DBSQL::CoordSystemAdaptor',
		 'CompressedSequence' => 'Bio::EnsEMBL::DBSQL::CompressedSequenceAdaptor',
		 'DBEntry'              => 'Bio::EnsEMBL::DBSQL::DBEntryAdaptor',
		 'DnaAlignFeature'      => 'Bio::EnsEMBL::DBSQL::DnaAlignFeatureAdaptor',
		 'DensityFeature'       => 'Bio::EnsEMBL::DBSQL::DensityFeatureAdaptor',
		 'DensityType'          => 'Bio::EnsEMBL::DBSQL::DensityTypeAdaptor',
		 'Exon'                 => 'Bio::EnsEMBL::DBSQL::ExonAdaptor',
		 'Gene'                 => 'Bio::EnsEMBL::DBSQL::GeneAdaptor',
		 'KaryotypeBand'        => 'Bio::EnsEMBL::DBSQL::KaryotypeBandAdaptor',
		 'Marker'               => 'Bio::EnsEMBL::Map::DBSQL::MarkerAdaptor',
		 'MarkerFeature'        =>
		 'Bio::EnsEMBL::Map::DBSQL::MarkerFeatureAdaptor',
		 'MetaCoordContainer'   => 'Bio::EnsEMBL::DBSQL::MetaCoordContainer',
		 'MiscSet'              => 'Bio::EnsEMBL::DBSQL::MiscSetAdaptor',
		 'MiscFeature'          => 'Bio::EnsEMBL::DBSQL::MiscFeatureAdaptor',
		 'PredictionTranscript' =>
		 'Bio::EnsEMBL::DBSQL::PredictionTranscriptAdaptor',
		 'PredictionExon'       => 'Bio::EnsEMBL::DBSQL::PredictionExonAdaptor',
		 'ProteinFeature'       => 'Bio::EnsEMBL::DBSQL::ProteinFeatureAdaptor',
		 'ProteinAlignFeature'  =>
		 'Bio::EnsEMBL::DBSQL::ProteinAlignFeatureAdaptor',
		 'ProxySNP'             => 'Bio::EnsEMBL::DBSQL::ProxySNPAdaptor',
		 'QtlFeature'           => 'Bio::EnsEMBL::Map::DBSQL::QtlFeatureAdaptor',
		 'Qtl'                  => 'Bio::EnsEMBL::Map::DBSQL::QtlAdaptor',
		 'RepeatConsensus'      => 'Bio::EnsEMBL::DBSQL::RepeatConsensusAdaptor',
		 'RepeatFeature'        => 'Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptor',
		 'Sequence'             => 'Bio::EnsEMBL::DBSQL::SequenceAdaptor',
		 'SimpleFeature'        => 'Bio::EnsEMBL::DBSQL::SimpleFeatureAdaptor',
		 'Slice'                => 'Bio::EnsEMBL::DBSQL::SliceAdaptor',
		 'SupportingFeature'    =>
		 'Bio::EnsEMBL::DBSQL::SupportingFeatureAdaptor',
		 'Transcript'           => 'Bio::EnsEMBL::DBSQL::TranscriptAdaptor',
		 'Translation'          => 'Bio::EnsEMBL::DBSQL::TranslationAdaptor' );


  foreach my $key (keys %pairs){
    Bio::EnsEMBL::Registry->add_adaptor($species, $group, $key, $pairs{$key});
  }

  #if dnadb has been set then for the follwing use it.
  foreach my $type (qw(Sequence AssemblyMapper KaryotypeBand RepeatFeature CoordSystem AssemblyExceptionFeature)){
    Bio::EnsEMBL::Registry->set_get_via_dnadb_if_set($species,$type);
  }

}


sub load_est{
  my ($class, @args) = @_;

  my ($species) = rearrange([qw(SPECIES)],@args);
  my $group = 'est';

  push (@args, '-group');
  push (@args, $group);

  my $dbc = new Bio::EnsEMBL::DBSQL::DBConnection(@args);

  my $dba = new_fast  Bio::EnsEMBL::DBSQL::DBAdaptor('-con' => $dbc);

  Bio::EnsEMBL::Registry->add_DBAdaptor($species, $group, $dba);

  my %pairs =  ( 'Analysis'             => 'Bio::EnsEMBL::DBSQL::AnalysisAdaptor',
		 'ArchiveStableId'      => 'Bio::EnsEMBL::DBSQL::ArchiveStableIdAdaptor',
		 'Attribute'            => 'Bio::EnsEMBL::DBSQL::AttributeAdaptor',
		 'AssemblyExceptionFeature' => 'Bio::EnsEMBL::DBSQL::AssemblyExceptionFeatureAdaptor',
		 'AssemblyMapper'       => 'Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor',
		 'Blast'                => 'Bio::EnsEMBL::External::BlastAdaptor',
		 'MetaContainer'        => 'Bio::EnsEMBL::DBSQL::MetaContainer',
		 'CoordSystem'   => 'Bio::EnsEMBL::DBSQL::CoordSystemAdaptor',
		 'CompressedSequence' => 'Bio::EnsEMBL::DBSQL::CompressedSequenceAdaptor',
		 'DBEntry'              => 'Bio::EnsEMBL::DBSQL::DBEntryAdaptor',
		 'DnaAlignFeature'      => 'Bio::EnsEMBL::DBSQL::DnaAlignFeatureAdaptor',
		 'DensityFeature'       => 'Bio::EnsEMBL::DBSQL::DensityFeatureAdaptor',
		 'DensityType'          => 'Bio::EnsEMBL::DBSQL::DensityTypeAdaptor',
		 'Exon'                 => 'Bio::EnsEMBL::DBSQL::ExonAdaptor',
		 'Gene'                 => 'Bio::EnsEMBL::DBSQL::GeneAdaptor',
		 'KaryotypeBand'        => 'Bio::EnsEMBL::DBSQL::KaryotypeBandAdaptor',
		 'Marker'               => 'Bio::EnsEMBL::Map::DBSQL::MarkerAdaptor',
		 'MarkerFeature'        =>
		 'Bio::EnsEMBL::Map::DBSQL::MarkerFeatureAdaptor',
		 'MetaCoordContainer'   => 'Bio::EnsEMBL::DBSQL::MetaCoordContainer',
		 'MiscSet'              => 'Bio::EnsEMBL::DBSQL::MiscSetAdaptor',
		 'MiscFeature'          => 'Bio::EnsEMBL::DBSQL::MiscFeatureAdaptor',
		 'PredictionTranscript' =>
		 'Bio::EnsEMBL::DBSQL::PredictionTranscriptAdaptor',
		 'PredictionExon'       => 'Bio::EnsEMBL::DBSQL::PredictionExonAdaptor',
		 'ProteinFeature'       => 'Bio::EnsEMBL::DBSQL::ProteinFeatureAdaptor',
		 'ProteinAlignFeature'  =>
		 'Bio::EnsEMBL::DBSQL::ProteinAlignFeatureAdaptor',
		 'ProxySNP'             => 'Bio::EnsEMBL::DBSQL::ProxySNPAdaptor',
		 'QtlFeature'           => 'Bio::EnsEMBL::Map::DBSQL::QtlFeatureAdaptor',
		 'Qtl'                  => 'Bio::EnsEMBL::Map::DBSQL::QtlAdaptor',
		 'RepeatConsensus'      => 'Bio::EnsEMBL::DBSQL::RepeatConsensusAdaptor',
		 'RepeatFeature'        => 'Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptor',
		 'Sequence'             => 'Bio::EnsEMBL::DBSQL::SequenceAdaptor',
		 'SimpleFeature'        => 'Bio::EnsEMBL::DBSQL::SimpleFeatureAdaptor',
		 'Slice'                => 'Bio::EnsEMBL::DBSQL::SliceAdaptor',
		 'SupportingFeature'    =>
		 'Bio::EnsEMBL::DBSQL::SupportingFeatureAdaptor',
		 'Transcript'           => 'Bio::EnsEMBL::DBSQL::TranscriptAdaptor',
		 'Translation'          => 'Bio::EnsEMBL::DBSQL::TranslationAdaptor' );


  foreach my $key (keys %pairs){
    Bio::EnsEMBL::Registry->add_adaptor($species, $group, $key, $pairs{$key});
  }

  #if dnadb has been set then for the follwing use it.
  foreach my $type (qw(Sequence AssemblyMapper KaryotypeBand RepeatFeature CoordSystem AssemblyExceptionFeature)){
    Bio::EnsEMBL::Registry->set_get_via_dnadb_if_set($species,$type);
  }
}


sub load_vega{
  my ($class, @args) = @_;

  my ($species) = rearrange([qw(SPECIES)],@args);
  my $group = 'vega';

  push (@args, '-group');
  push (@args, $group);

  my $dbc = new Bio::EnsEMBL::DBSQL::DBConnection(@args);

  my $dba = new_fast  Bio::EnsEMBL::DBSQL::DBAdaptor('-con' => $dbc);

  Bio::EnsEMBL::Registry->add_DBAdaptor($species, $group, $dba);

  my %pairs =  ( 'Analysis'             => 'Bio::EnsEMBL::DBSQL::AnalysisAdaptor',
		 'ArchiveStableId'      => 'Bio::EnsEMBL::DBSQL::ArchiveStableIdAdaptor',
		 'Attribute'            => 'Bio::EnsEMBL::DBSQL::AttributeAdaptor',
		 'AssemblyExceptionFeature' => 'Bio::EnsEMBL::DBSQL::AssemblyExceptionFeatureAdaptor',
		 'AssemblyMapper'       => 'Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor',
		 'Blast'                => 'Bio::EnsEMBL::External::BlastAdaptor',
		 'MetaContainer'        => 'Bio::EnsEMBL::DBSQL::MetaContainer',
		 'CoordSystem'   => 'Bio::EnsEMBL::DBSQL::CoordSystemAdaptor',
		 'CompressedSequence' => 'Bio::EnsEMBL::DBSQL::CompressedSequenceAdaptor',
		 'DBEntry'              => 'Bio::EnsEMBL::DBSQL::DBEntryAdaptor',
		 'DnaAlignFeature'      => 'Bio::EnsEMBL::DBSQL::DnaAlignFeatureAdaptor',
		 'DensityFeature'       => 'Bio::EnsEMBL::DBSQL::DensityFeatureAdaptor',
		 'DensityType'          => 'Bio::EnsEMBL::DBSQL::DensityTypeAdaptor',
		 'Exon'                 => 'Bio::EnsEMBL::DBSQL::ExonAdaptor',
		 'Gene'                 => 'Bio::EnsEMBL::DBSQL::GeneAdaptor',
		 'KaryotypeBand'        => 'Bio::EnsEMBL::DBSQL::KaryotypeBandAdaptor',
		 'Marker'               => 'Bio::EnsEMBL::Map::DBSQL::MarkerAdaptor',
		 'MarkerFeature'        =>
		 'Bio::EnsEMBL::Map::DBSQL::MarkerFeatureAdaptor',
		 'MetaCoordContainer'   => 'Bio::EnsEMBL::DBSQL::MetaCoordContainer',
		 'MiscSet'              => 'Bio::EnsEMBL::DBSQL::MiscSetAdaptor',
		 'MiscFeature'          => 'Bio::EnsEMBL::DBSQL::MiscFeatureAdaptor',
		 'PredictionTranscript' =>
		 'Bio::EnsEMBL::DBSQL::PredictionTranscriptAdaptor',
		 'PredictionExon'       => 'Bio::EnsEMBL::DBSQL::PredictionExonAdaptor',
		 'ProteinFeature'       => 'Bio::EnsEMBL::DBSQL::ProteinFeatureAdaptor',
		 'ProteinAlignFeature'  =>
		 'Bio::EnsEMBL::DBSQL::ProteinAlignFeatureAdaptor',
		 'ProxySNP'             => 'Bio::EnsEMBL::DBSQL::ProxySNPAdaptor',
		 'QtlFeature'           => 'Bio::EnsEMBL::Map::DBSQL::QtlFeatureAdaptor',
		 'Qtl'                  => 'Bio::EnsEMBL::Map::DBSQL::QtlAdaptor',
		 'RepeatConsensus'      => 'Bio::EnsEMBL::DBSQL::RepeatConsensusAdaptor',
		 'RepeatFeature'        => 'Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptor',
		 'Sequence'             => 'Bio::EnsEMBL::DBSQL::SequenceAdaptor',
		 'SimpleFeature'        => 'Bio::EnsEMBL::DBSQL::SimpleFeatureAdaptor',
		 'Slice'                => 'Bio::EnsEMBL::DBSQL::SliceAdaptor',
		 'SupportingFeature'    =>
		 'Bio::EnsEMBL::DBSQL::SupportingFeatureAdaptor',
		 'Transcript'           => 'Bio::EnsEMBL::DBSQL::TranscriptAdaptor',
		 'Translation'          => 'Bio::EnsEMBL::DBSQL::TranslationAdaptor' );
  foreach my $key (keys %pairs){
    Bio::EnsEMBL::Registry->add_adaptor($species, $group, $key, $pairs{$key});
  }
}


sub load_compara{
  my ($class, @args) = @_;
  require Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;

  my ($species) = rearrange([qw(SPECIES)],@args);
  my $group = 'compara';


  push (@args, '-group');
  push (@args, $group);

  my $dbc = new Bio::EnsEMBL::DBSQL::DBConnection(@args);

  my $dba = new_fast  Bio::EnsEMBL::Compara::DBSQL::DBAdaptor('-con' => $dbc);

  Bio::EnsEMBL::Registry->add_DBAdaptor($species, $group, $dba);

  my %pairs =  ( "MetaContainer" => "Bio::EnsEMBL::DBSQL::MetaContainer",
	      'SyntenyRegion'   => 'Bio::EnsEMBL::Compara::DBSQL::SyntenyRegionAdaptor',
	      "DnaAlignFeature" => "Bio::EnsEMBL::Compara::DBSQL::DnaAlignFeatureAdaptor",
	      "Synteny"         => "Bio::EnsEMBL::Compara::DBSQL::SyntenyAdaptor",
	      "GenomeDB"        => "Bio::EnsEMBL::Compara::DBSQL::GenomeDBAdaptor",
	      "DnaFrag" => "Bio::EnsEMBL::Compara::DBSQL::DnaFragAdaptor",
	      "GenomicAlign" => "Bio::EnsEMBL::Compara::DBSQL::GenomicAlignAdaptor",
	      "Homology" => "Bio::EnsEMBL::Compara::DBSQL::HomologyAdaptor",
	      "Family" => "Bio::EnsEMBL::Compara::DBSQL::FamilyAdaptor",
	      "Domain" => "Bio::EnsEMBL::Compara::DBSQL::DomainAdaptor",
	      "Subset" => "Bio::EnsEMBL::Compara::DBSQL::SubsetAdaptor",
	      "Member" => "Bio::EnsEMBL::Compara::DBSQL::MemberAdaptor",
	      "Attribute" => "Bio::EnsEMBL::Compara::DBSQL::AttributeAdaptor",
	      "Taxon" => "Bio::EnsEMBL::Compara::DBSQL::TaxonAdaptor",
	      "PeptideAlignFeature" => "Bio::EnsEMBL::Compara::DBSQL::PeptideAlignFeatureAdaptor",
	      "AnalysisAdaptor" => "Bio::EnsEMBL::DBSQL::AnalysisAdaptor"
        );

  foreach my $key (keys %pairs){

    Bio::EnsEMBL::Registry->add_adaptor($species, $group, $key, $pairs{$key});
  }

}

sub load_hive{
  my ($class, @args) = @_;
  require Bio::EnsEMBL::Hive::DBSQL::DBAdaptor;

  my ($species) = rearrange([qw(SPECIES)],@args);
  my $group = 'hive';


  push (@args, '-group');
  push (@args, $group);

  my $dbc = new Bio::EnsEMBL::DBSQL::DBConnection(@args);

  my $dba = new_fast  Bio::EnsEMBL::Hive::DBSQL::DBAdaptor('-con' => $dbc);

  Bio::EnsEMBL::Registry->add_DBAdaptor($species, $group, $dba);

  my %pairs =  (
	      "Analysis"         => "Bio::EnsEMBL::DBSQL::AnalysisAdaptor",
	      "Queen"            => "Bio::EnsEMBL::Hive::Queen",
	      "AnalysisJob"      => "Bio::EnsEMBL::Hive::DBSQL::AnalysisJobAdaptor",
	      "AnalysisStats"    => "Bio::EnsEMBL::Hive::DBSQL::AnalysisStatsAdaptor",
	      "AnalysisCtrlRule" => "Bio::EnsEMBL::Hive::DBSQL::AnalysisCtrlRuleAdaptor",
	      "DataflowRule"     => "Bio::EnsEMBL::Hive::DBSQL::DataflowRuleAdaptor",
	      "SimpleRule"       => "Bio::EnsEMBL::Hive::DBSQL::SimpleRuleAdaptor");

  foreach my $key (keys %pairs){

    Bio::EnsEMBL::Registry->add_adaptor($species, $group, $key, $pairs{$key});
  }

}

sub load_go{

 my ($class, @args) = @_;

 my ($species) = rearrange([qw(SPECIES)],@args);

 require Bio::EnsEMBL::ExternalData::GO::GOAdaptor;
 
 my $adap = new Bio::EnsEMBL::ExternalData::GO::GOAdaptor(@args);
 Bio::EnsEMBL::Registry->add_alias($species,$species);

 #not a dba adaptor but stored as one
 Bio::EnsEMBL::Registry->add_DBAdaptor($species,"go",$adap );
}


sub load_haplotype{
  my ($class, @args) = @_;

  require Bio::EnsEMBL::ExternalData::Haplotype::DBAdaptor;

  my ($species) = rearrange([qw(SPECIES)],@args);

  my $group = 'haplotype';

  push (@args, '-group');
  push (@args, $group);

  my $dbc = new Bio::EnsEMBL::DBSQL::DBConnection(@args);

  my $dba = new_fast Bio::EnsEMBL::ExternalData::Haplotype::DBAdaptor('-con' => $dbc);


  Bio::EnsEMBL::Registry->add_DBAdaptor($species, $group, $dba);

  my %pairs = ('Haplotype' => 'Bio::EnsEMBL::ExternalData::Haplotype::HaplotypeAdaptor');

  foreach my $key (keys %pairs){
    Bio::EnsEMBL::Registry->add_adaptor($species, $group, $key, $pairs{$key});
  }

}


sub load_pipeline{
  my ($class, @args) = @_;

  require Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;

  my ($species) = rearrange([qw(SPECIES)],@args);

  my $group = 'pipeline';

  push (@args, '-group');
  push (@args, $group);

  my $dbc = new Bio::EnsEMBL::DBSQL::DBConnection(@args);

  my $dba = new_fast Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor('-con' => $dbc);


  Bio::EnsEMBL::Registry->add_DBAdaptor($species, $group, $dba);

  my %pairs =  ( 'ArchiveStableId'      => 'Bio::EnsEMBL::DBSQL::ArchiveStableIdAdaptor',
		 'Attribute'            => 'Bio::EnsEMBL::DBSQL::AttributeAdaptor',
		 'AssemblyExceptionFeature' => 'Bio::EnsEMBL::DBSQL::AssemblyExceptionFeatureAdaptor',
		 'AssemblyMapper'       => 'Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor',
		 'MetaContainer'        => 'Bio::EnsEMBL::DBSQL::MetaContainer',
		 'CoordSystem'   => 'Bio::EnsEMBL::DBSQL::CoordSystemAdaptor',
		 'CompressedSequence' => 'Bio::EnsEMBL::DBSQL::CompressedSequenceAdaptor',
		 'DBEntry'              => 'Bio::EnsEMBL::DBSQL::DBEntryAdaptor',
		 'DnaAlignFeature'      => 'Bio::EnsEMBL::DBSQL::DnaAlignFeatureAdaptor',
		 'DensityFeature'       => 'Bio::EnsEMBL::DBSQL::DensityFeatureAdaptor',
		 'DensityType'          => 'Bio::EnsEMBL::DBSQL::DensityTypeAdaptor',
		 'Exon'                 => 'Bio::EnsEMBL::DBSQL::ExonAdaptor',
		 'Gene'                 => 'Bio::EnsEMBL::DBSQL::GeneAdaptor',
		 'KaryotypeBand'        => 'Bio::EnsEMBL::DBSQL::KaryotypeBandAdaptor',
		 'Marker'               => 'Bio::EnsEMBL::Map::DBSQL::MarkerAdaptor',
		 'MarkerFeature'        =>
		 'Bio::EnsEMBL::Map::DBSQL::MarkerFeatureAdaptor',
		 'MetaCoordContainer'   => 'Bio::EnsEMBL::DBSQL::MetaCoordContainer',
		 'MiscSet'              => 'Bio::EnsEMBL::DBSQL::MiscSetAdaptor',
		 'MiscFeature'          => 'Bio::EnsEMBL::DBSQL::MiscFeatureAdaptor',
		 'PredictionTranscript' => 'Bio::EnsEMBL::DBSQL::PredictionTranscriptAdaptor',
		 'PredictionExon'       => 'Bio::EnsEMBL::DBSQL::PredictionExonAdaptor',
		 'ProteinFeature'       => 'Bio::EnsEMBL::DBSQL::ProteinFeatureAdaptor',
		 'ProteinAlignFeature'  =>
		 'Bio::EnsEMBL::DBSQL::ProteinAlignFeatureAdaptor',
		 'ProxySNP'             => 'Bio::EnsEMBL::DBSQL::ProxySNPAdaptor',
		 'QtlFeature'           => 'Bio::EnsEMBL::Map::DBSQL::QtlFeatureAdaptor',
		 'Qtl'                  => 'Bio::EnsEMBL::Map::DBSQL::QtlAdaptor',
		 'RepeatConsensus'      => 'Bio::EnsEMBL::DBSQL::RepeatConsensusAdaptor',
		 'RepeatFeature'        => 'Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptor',
		 'Sequence'             => 'Bio::EnsEMBL::DBSQL::SequenceAdaptor',
		 'SimpleFeature'        => 'Bio::EnsEMBL::DBSQL::SimpleFeatureAdaptor',
		 'Slice'                => 'Bio::EnsEMBL::DBSQL::SliceAdaptor',
		 'SupportingFeature'    =>
		 'Bio::EnsEMBL::DBSQL::SupportingFeatureAdaptor',
		 'Transcript'           => 'Bio::EnsEMBL::DBSQL::TranscriptAdaptor',
		 'Translation'          => 'Bio::EnsEMBL::DBSQL::TranslationAdaptor',
		 'Analysis'           => 'Bio::EnsEMBL::Pipeline::DBSQL::AnalysisAdaptor',
		 'Job'                => 'Bio::EnsEMBL::Pipeline::DBSQL::JobAdaptor',
		 'PmatchFeature'      => 'Bio::EnsEMBL::Pipeline::DBSQL::PmatchFeatureAdaptor',
		 'Rule'               => 'Bio::EnsEMBL::Pipeline::DBSQL::RuleAdaptor',
		 'StateInfoContainer' => 'Bio::EnsEMBL::Pipeline::DBSQL::StateInfoContainer');

  foreach my $key (keys %pairs){
    Bio::EnsEMBL::Registry->add_adaptor($species, $group, $key, $pairs{$key});
  }

}


sub dnadb_add{
  my $class = shift;
  my ($dnaspecies, $dnagroup, $species, $group) = @_;

  my $dnadb =  Bio::EnsEMBL::Registry->get_DBAdaptor($dnaspecies, $dnagroup);
  my $featdb = Bio::EnsEMBL::Registry->get_DBAdaptor($species, $group);

  $featdb->dnadb($dnadb);
}


sub add_alias{
  my ($class, @args) = @_;
  my ($species, $aliases) = rearrange([qw(SPECIES ALIAS)],@args);

 #make sure it exists itself
  Bio::EnsEMBL::Registry->add_alias($species,$species);

  if($aliases){
    foreach my $ali (@$aliases){
      Bio::EnsEMBL::Registry->add_alias($species,$ali);
    }
  }

}


sub get_alias{
  my ($class, $key) = @_;
  return Bio::EnsEMBL::Registry->get_alias($key);
}

sub load_privgene{
  my ($class, @args) = @_;

  my ($species) = rearrange([qw(SPECIES)],@args);
  my $group = 'privgene';

  push (@args, '-group');
  push (@args, $group);

  my ($spec,$gro) = Bio::EnsEMBL::Registry->check_if_already_there(@args);
  if($spec){
    print STDERR "already defined returning??????\n";
    return;
  }

  my $dbc = new Bio::EnsEMBL::DBSQL::DBConnection(@args);

  my $dba = new_fast  Bio::EnsEMBL::DBSQL::DBAdaptor('-con' => $dbc);

  Bio::EnsEMBL::Registry->add_DBAdaptor($species, $group, $dba);

  my %pairs =  ( 'Analysis'             => 'Bio::EnsEMBL::DBSQL::AnalysisAdaptor',
		 'ArchiveStableId'      => 'Bio::EnsEMBL::DBSQL::ArchiveStableIdAdaptor',
		 'Attribute'            => 'Bio::EnsEMBL::DBSQL::AttributeAdaptor',
		 'AssemblyExceptionFeature' => 'Bio::EnsEMBL::DBSQL::AssemblyExceptionFeatureAdaptor',
		 'AssemblyMapper'       => 'Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor',
		 'Blast'                => 'Bio::EnsEMBL::External::BlastAdaptor',
		 'MetaContainer'        => 'Bio::EnsEMBL::DBSQL::MetaContainer',
		 'CoordSystem'   => 'Bio::EnsEMBL::DBSQL::CoordSystemAdaptor',
		 'CompressedSequence' => 'Bio::EnsEMBL::DBSQL::CompressedSequenceAdaptor',
		 'DBEntry'              => 'Bio::EnsEMBL::DBSQL::DBEntryAdaptor',
		 'DnaAlignFeature'      => 'Bio::EnsEMBL::DBSQL::DnaAlignFeatureAdaptor',
		 'DensityFeature'       => 'Bio::EnsEMBL::DBSQL::DensityFeatureAdaptor',
		 'DensityType'          => 'Bio::EnsEMBL::DBSQL::DensityTypeAdaptor',
		 'Exon'                 => 'Bio::EnsEMBL::DBSQL::ExonAdaptor',
		 'Gene'                 => 'Bio::EnsEMBL::DBSQL::GeneAdaptor',
		 'KaryotypeBand'        => 'Bio::EnsEMBL::DBSQL::KaryotypeBandAdaptor',
		 'Marker'               => 'Bio::EnsEMBL::Map::DBSQL::MarkerAdaptor',
		 'MarkerFeature'        =>
		 'Bio::EnsEMBL::Map::DBSQL::MarkerFeatureAdaptor',
		 'MetaCoordContainer'   => 'Bio::EnsEMBL::DBSQL::MetaCoordContainer',
		 'MiscSet'              => 'Bio::EnsEMBL::DBSQL::MiscSetAdaptor',
		 'MiscFeature'          => 'Bio::EnsEMBL::DBSQL::MiscFeatureAdaptor',
		 'PredictionTranscript' =>
		 'Bio::EnsEMBL::DBSQL::PredictionTranscriptAdaptor',
		 'PredictionExon'       => 'Bio::EnsEMBL::DBSQL::PredictionExonAdaptor',
		 'ProteinFeature'       => 'Bio::EnsEMBL::DBSQL::ProteinFeatureAdaptor',
		 'ProteinAlignFeature'  =>
		 'Bio::EnsEMBL::DBSQL::ProteinAlignFeatureAdaptor',
		 'ProxySNP'             => 'Bio::EnsEMBL::DBSQL::ProxySNPAdaptor',
		 'QtlFeature'           => 'Bio::EnsEMBL::Map::DBSQL::QtlFeatureAdaptor',
		 'Qtl'                  => 'Bio::EnsEMBL::Map::DBSQL::QtlAdaptor',
		 'RepeatConsensus'      => 'Bio::EnsEMBL::DBSQL::RepeatConsensusAdaptor',
		 'RepeatFeature'        => 'Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptor',
		 'Sequence'             => 'Bio::EnsEMBL::DBSQL::SequenceAdaptor',
		 'SimpleFeature'        => 'Bio::EnsEMBL::DBSQL::SimpleFeatureAdaptor',
		 'Slice'                => 'Bio::EnsEMBL::DBSQL::SliceAdaptor',
		 'SupportingFeature'    =>
		 'Bio::EnsEMBL::DBSQL::SupportingFeatureAdaptor',
		 'Transcript'           => 'Bio::EnsEMBL::DBSQL::TranscriptAdaptor',
		 'Translation'          => 'Bio::EnsEMBL::DBSQL::TranslationAdaptor' );

  foreach my $key (keys %pairs){
    Bio::EnsEMBL::Registry->add_adaptor($species, $group, $key, $pairs{$key});
  }

  #if dnadb has been set then for the follwing use it.
  foreach my $type (qw(Sequence AssemblyMapper KaryotypeBand RepeatFeature CoordSystem AssemblyExceptionFeature)){
    Bio::EnsEMBL::Registry->set_get_via_dnadb_if_set($species,$type);
  }
}

#
# overwrite/load new types. Done this way to enable no changes to CVS for
# external users. External users should add there own "GROUPS" in the file
# User_defined_load.
#

eval{ require Bio::EnsEMBL::Utils::User_defined_load };
if ($@){ print STDERR  "No user defined loads\n"; }

1;
