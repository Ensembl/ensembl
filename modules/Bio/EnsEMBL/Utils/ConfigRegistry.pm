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


Bio::EnsEMBL::Utils::ConfigRegistry->load_core($dba );
 

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
my $reg = "Bio::EnsEMBL::Registry";
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

use Bio::EnsEMBL::Utils::Exception qw(warning throw  deprecate stack_trace_dump);
#use vars qw(@ISA @EXPORT_OK);
#@ISA = qw(Exporter);

#@EXPORT_OK = qw(&load_core &load_estgene);



=head2 load_core, load_estgene, load_vega, load_compara, load_pipeline, load_SNP, load_lite
  Arg [1]    : DBAdaptor with DBConnection alredy attached
  Returntype : DBAdaptor;
  Exceptions : none

=cut

sub gen_load{
  my ($dba) = @_;
  my $config_sub;

# at some point we hope to set the group in the DBadaptor
# hence this long check etc should be simpler

  if($dba->isa('Bio::EnsEMBL::Compara::DBSQL::DBAdaptor')){
    if(!defined($dba->group())){
      $dba->group('compara');
    }
    $config_sub =  \&Bio::EnsEMBL::Utils::ConfigRegistry::load_compara; 
  }
  elsif($dba->isa('Bio::EnsEMBL::Lite::DBAdaptor')){
    if(!defined($dba->group())){
      $dba->group('lite');
    }
    $config_sub =  \&Bio::EnsEMBL::Utils::ConfigRegistry::load_lite;
  }
  elsif($dba->isa('Bio::EnsEMBL::External::BlastAdaptor')){
    if(!defined($dba->group())){
      $dba->group('blast');
    }
    $config_sub =  \&Bio::EnsEMBL::Utils::ConfigRegistry::load_blast;
  }
  elsif($dba->isa('Bio::EnsEMBL::ExternalData::SNPSQL::DBAdaptor')){
    if(!defined($dba->group())){
      $dba->group('SNP');
    } 
    $config_sub =  \&Bio::EnsEMBL::Utils::ConfigRegistry::load_SNP;
  }
  elsif($dba->isa('Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor')){
    if(!defined($dba->group())){
      $dba->group('pipeline');
    }
    $config_sub =  \&Bio::EnsEMBL::Utils::ConfigRegistry::load_pipeline;
  }
  elsif($dba->isa('Bio::EnsEMBL::Hive::DBSQL::DBAdaptor')){
    if(!defined($dba->group())){
      $dba->group('hive');
    }
    $config_sub =  \&Bio::EnsEMBL::Utils::ConfigRegistry::load_hive;
  }
  elsif($dba->isa('Bio::EnsEMBL::ExternalData::Haplotype::DBAdaptor')){
    if(!defined($dba->group())){
      $dba->group('haplotype');
    }
    $config_sub =  \&Bio::EnsEMBL::Utils::ConfigRegistry::load_haplotype;    
  }
  elsif($dba->isa('Bio::EnsEMBL::DBSQL::DBAdaptor')){
    #vega uses the core DBAdaptor so test if vega is in the dbname
    if($dba->dbc->dbname() =~ /vega/){
      if(!defined($dba->group())){
	$dba->group('vega');
      }
      $config_sub =  \&Bio::EnsEMBL::Utils::ConfigRegistry::load_vega;
     }
    else{
      if(!defined($dba->group())){
	$dba->group('core');
      }
      $config_sub =  \&Bio::EnsEMBL::Utils::ConfigRegistry::load_core;
    }
  }
  else{
    throw("Unknown DBAdaptor type $dba\n");
  }


  # return if the connection and species, group are the same
  my $db_reg;
  if($reg->get_alias($dba->species,"no throw")){
    $db_reg = $reg->get_DBAdaptor($dba->species,$dba->group);
  }
  if(defined($db_reg)){
    if($dba->dbc->equals($db_reg->dbc)){
      return $db_reg;
    }
    else{ # diff conn details
      $dba->species(find_unique_species($dba->species,$dba->group));
    }
  }

  Bio::EnsEMBL::Registry->add_DBAdaptor($dba->species(), $dba->group(), $dba);

  #call the loading subroutine. (add the adaptors to the DBAdaptor)
  &{$config_sub}($dba);
   
  return $dba;
}



sub find_unique_species{
  my ($species, $group) = @_;

  $reg->add_alias($species,$species);

  my $i = 0;
  my $free =0;
  while(!$free){
    if($i == 0){
      if(!defined($reg->get_DBAdaptor($species, $group))){
	$free =1;
	$i ="";
      }
      else{
	$i = 1;
      }
    }
    else{
      $reg->add_alias($species.$i,$species.$i); #set needed self alias
      if(!defined($reg->get_DBAdaptor($species.$i, $group))){
	$free =1;
      }
      else{
	$i++;
      }
    }
  }
  
  $species .= $i;
  return ($species);
}

sub load_core{
  my ($dba) = @_;


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
		 'SNP'             => 'Bio::EnsEMBL::DBSQL::ProxySNPAdaptor',
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

    Bio::EnsEMBL::Registry->add_adaptor($dba->species, $dba->group, $key, $pairs{$key});
  }


  foreach my $type (qw(Sequence AssemblyMapper KaryotypeBand RepeatFeature CoordSystem AssemblyExceptionFeature)){
    Bio::EnsEMBL::Registry->set_get_via_dnadb_if_set($dba->species,$type);
  }

}


sub load_lite{
    my ($dba) = @_;
  
  Bio::EnsEMBL::Registry->add_DBAdaptor($dba->species, $dba->group, $dba);
    
  Bio::EnsEMBL::Registry->add_adaptor($dba->species, $dba->group, "SNP", "Bio::EnsEMBL::Lite::SNPAdaptor" );
}

sub load_SNP{
  my ($dba) = @_;
  Bio::EnsEMBL::Registry->add_DBAdaptor($dba->species, $dba->group, $dba);
    
  Bio::EnsEMBL::Registry->add_adaptor($dba->species, $dba->group, "SNP", 
				      "Bio::EnsEMBL::ExternalData::SNPSQL::SNPAdaptor" );
}

sub attach_database{
  my ($class, $species, $core, $name1) = @_;

  print STDERR "attach_databse called ".caller(). "\n";

  my $first =  Bio::EnsEMBL::Registry->get_DBAdaptor($species,$name1);
  my $coredb = Bio::EnsEMBL::Registry->get_DBAdaptor($species,$core);

  Bio::EnsEMBL::Registry->add_db($coredb,$name1, $first);

}

sub attach_dna{
  my ($class, $species, $main, $attach) = @_; 
  print STDERR "attach_dna called ".caller(). "\n";

  my $no_seq =  Bio::EnsEMBL::Registry->get_DBAdaptor($species,$main);
  my $seq = Bio::EnsEMBL::Registry->get_DBAdaptor($species,$attach);

  Bio::EnsEMBL::Registry->add_DNAAdaptor($species,$no_seq->group,$seq);
}


sub load_blast{
    my ($dba) = @_;

  Bio::EnsEMBL::Registry->add_DBAdaptor($dba->species, $dba->group, $dba);
  Bio::EnsEMBL::Registry->add_adaptor($dba->species, $dba->group, 'Blast', 
				      "Bio::EnsEMBL::External::BlastAdaptor");
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
  my ($dba) = @_;

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
    Bio::EnsEMBL::Registry->add_adaptor($dba->species, $dba->group, $key, $pairs{$key});
  }

  #if dnadb has been set then for the follwing use it.
  foreach my $type (qw(Sequence AssemblyMapper KaryotypeBand RepeatFeature CoordSystem AssemblyExceptionFeature)){
    Bio::EnsEMBL::Registry->set_get_via_dnadb_if_set($dba->species,$type);
  }

}


sub load_est{
  my ($dba) = @_;

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
    Bio::EnsEMBL::Registry->add_adaptor($dba->species, $dba->group, $key, $pairs{$key});
  }

  #if dnadb has been set then for the follwing use it.
  foreach my $type (qw(Sequence AssemblyMapper KaryotypeBand RepeatFeature CoordSystem AssemblyExceptionFeature)){
    Bio::EnsEMBL::Registry->set_get_via_dnadb_if_set($dba->species,$type);
  }
}


sub load_vega{
  my ($dba) = @_;

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
    Bio::EnsEMBL::Registry->add_adaptor($dba->species, $dba->group, $key, $pairs{$key});
  }
}


sub load_compara{
  my ($dba) = @_;

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
	      "Analysis" => "Bio::EnsEMBL::DBSQL::AnalysisAdaptor"
        );

  foreach my $key (keys %pairs){

    Bio::EnsEMBL::Registry->add_adaptor($dba->species, $dba->group, $key, $pairs{$key});
  }

}

sub load_hive{
  my ($dba) = @_;

  my %pairs =  (
        "MetaContainer"    => 'Bio::EnsEMBL::DBSQL::MetaContainer',
	      "Analysis"         => "Bio::EnsEMBL::DBSQL::AnalysisAdaptor",
	      "Queen"            => "Bio::EnsEMBL::Hive::Queen",
	      "AnalysisJob"      => "Bio::EnsEMBL::Hive::DBSQL::AnalysisJobAdaptor",
	      "AnalysisStats"    => "Bio::EnsEMBL::Hive::DBSQL::AnalysisStatsAdaptor",
	      "AnalysisCtrlRule" => "Bio::EnsEMBL::Hive::DBSQL::AnalysisCtrlRuleAdaptor",
	      "DataflowRule"     => "Bio::EnsEMBL::Hive::DBSQL::DataflowRuleAdaptor",
	      "SimpleRule"       => "Bio::EnsEMBL::Hive::DBSQL::SimpleRuleAdaptor");

  foreach my $key (keys %pairs){

    Bio::EnsEMBL::Registry->add_adaptor($dba->species, $dba->group, $key, $pairs{$key});
  }

}

sub load_go{

 my ($class, $dba) = @_;

 # shouldnt go into the registry ... sorry
}


sub load_haplotype{
  my ($dba) = @_;

  require Bio::EnsEMBL::ExternalData::Haplotype::DBAdaptor;

  my %pairs = ('Haplotype' => 'Bio::EnsEMBL::ExternalData::Haplotype::HaplotypeAdaptor');

  foreach my $key (keys %pairs){
    Bio::EnsEMBL::Registry->add_adaptor($dba->species, $dba->group, $key, $pairs{$key});
  }

}


sub load_pipeline{
  my ($dba) = @_;

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
    Bio::EnsEMBL::Registry->add_adaptor($dba->species, $dba->group, $key, $pairs{$key});
  }

}


sub dnadb_add{
  my $class = shift;
  my ($dnaspecies, $dnagroup, $species, $group) = @_;
  print STDERR "dnadb_add called ".caller(). "\n";

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


#
# overwrite/load new types. Done this way to enable no changes to CVS for
# external users. External users should add there own "GROUPS" in the file
# User_defined_load.
#

eval{ require Bio::EnsEMBL::Utils::User_defined_load };
#if ($@){ print STDERR  "No user defined loads\n"; }

1;
