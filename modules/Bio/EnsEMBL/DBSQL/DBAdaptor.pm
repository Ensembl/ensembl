
=head1 NAME - Bio::EnsEMBL::DBSQL::DBAdaptor

=head1 SYNOPSIS

    $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -user   => 'root',
        -dbname => 'pog',
        -host   => 'caldy',
        -driver => 'mysql'
        );

    $gene_adaptor = $db->get_GeneAdaptor();

    $gene = $gene_adaptor()->fetch_by_stable_id($stable_id);

    $slice = $db->get_SliceAdaptor()->fetch_by_chr_start_end('X', 1, 10000);

=head1 DESCRIPTION

This is the primary interface to an EnsEMBL database. It maintains an active
connection to the database and allows for the retrieval of ObjectAdaptors,
via a set of get_XxxxAdaptor methods (where Xxxx is the type of adaptor).

ObjectAdaptors can then be used to obtain objects and actual information
from the database.


=head1 CONTACT

Post questions to the EnsEMBL development list <ensembl-dev@ebi.ac.uk>

=head1 METHODS

The rest of the documentation details each of the object methods. Internal
methods are usually preceded with a _

=cut


package Bio::EnsEMBL::DBSQL::DBAdaptor;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;

use Bio::EnsEMBL::Utils::Exception qw(throw warning deprecate);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

@ISA = qw(Bio::EnsEMBL::DBSQL::DBConnection);



#Override constructor inherited by Bio::EnsEMBL::DBSQL::DBConnection


=head2 new

  Arg [-DNADB]: (optional) Bio::EnsEMBL::DBSQL::DBAdaptor DNADB 
               All sequence, assembly, contig information etc, will be
               retrieved from this database instead.
  Arg [..]   : Other args are passed to superclass
               Bio::EnsEMBL::DBSQL::DBConnection
  Example    : $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
						    -user   => 'root',
						    -dbname => 'pog',
						    -host   => 'caldy',
						    -driver => 'mysql' );
  Description: Constructor for DBAdaptor.
  Returntype : Bio::EnsEMBL::DBSQL::DBAdaptor
  Exceptions : none
  Caller     : general

=cut

sub new {
  my($class, @args) = @_;

  #call superclass constructor
  my $self = $class->SUPER::new(@args);

  my ( $dnadb ) = rearrange([qw(DNADB)],@args);

  if(defined $dnadb) {
    $self->dnadb($dnadb);
  }

	# $self here is actually a Container object
	# so need to call _obj to get the DBAdaptor
	$self->_obj->{'default_module'} =
    { 'Analysis'             => 'Bio::EnsEMBL::DBSQL::AnalysisAdaptor',
      'ArchiveStableId'      => 'Bio::EnsEMBL::DBSQL::ArchiveStableIdAdaptor',
      'AssemblyMapper'       => 'Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor',
      'Blast'                => 'Bio::EnsEMBL::External::BlastAdaptor',
      'Chromosome'           => 'Bio::EnsEMBL::DBSQL::ChromosomeAdaptor',
      'Clone'                => 'Bio::EnsEMBL::DBSQL::CloneAdaptor',
      'CoordSystem'   => 'Bio::EnsEMBL::DBSQL::CoordSystemAdaptor',
      'DBEntry'              => 'Bio::EnsEMBL::DBSQL::DBEntryAdaptor',
      'DnaAlignFeature'      => 'Bio::EnsEMBL::DBSQL::DnaAlignFeatureAdaptor',
      'Exon'                 => 'Bio::EnsEMBL::DBSQL::ExonAdaptor',
      'Gene'                 => 'Bio::EnsEMBL::DBSQL::GeneAdaptor',
      'KaryotypeBand'        => 'Bio::EnsEMBL::DBSQL::KaryotypeBandAdaptor',
      'MapFrag'              => 'Bio::EnsEMBL::DBSQL::MapFragAdaptor',
      'Marker'               => 'Bio::EnsEMBL::Map::DBSQL::MarkerAdaptor',
      'MarkerFeature'        =>
          'Bio::EnsEMBL::Map::DBSQL::MarkerFeatureAdaptor',
      'MetaContainer'        => 'Bio::EnsEMBL::DBSQL::MetaContainer',
      'MiscSet'              => 'Bio::EnsEMBL::DBSQL::MiscSetAdaptor',
      'MiscFeature'          => 'Bio::EnsEMBL::DBSQL::MiscFeatureAdaptor',
      'PredictionTranscript' =>
           'Bio::EnsEMBL::DBSQL::PredictionTranscriptAdaptor',
      'PredictionExon'       => 'Bio::EnsEMBL::DBSQL::PredictionExonAdaptor',
      'ProteinFeature'       => 'Bio::EnsEMBL::DBSQL::ProteinFeatureAdaptor',
      'Protein'              => 'Bio::EnsEMBL::DBSQL::ProteinAdaptor',
      'ProteinAlignFeature'  =>
           'Bio::EnsEMBL::DBSQL::ProteinAlignFeatureAdaptor',
      'ProxySNP'             => 'Bio::EnsEMBL::DBSQL::ProxySNPAdaptor',
      'ProxyGene'            => 'Bio::EnsEMBL::DBSQL::ProxyGeneAdaptor',
      'ProxyRepeatFeature'   =>
          'Bio::EnsEMBL::DBSQL::ProxyRepeatFeatureAdaptor',
      'QtlFeature'           => 'Bio::EnsEMBL::Map::DBSQL::QtlFeatureAdaptor',
      'Qtl'                  => 'Bio::EnsEMBL::Map::DBSQL::QtlAdaptor',
      'RawContig'            => 'Bio::EnsEMBL::DBSQL::RawContigAdaptor',
      'RepeatConsensus'      => 'Bio::EnsEMBL::DBSQL::RepeatConsensusAdaptor',
      'RepeatFeature'        => 'Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptor',
      'Sequence'             => 'Bio::EnsEMBL::DBSQL::SequenceAdaptor',
      'SimpleFeature'        => 'Bio::EnsEMBL::DBSQL::SimpleFeatureAdaptor',
      'Slice'                => 'Bio::EnsEMBL::DBSQL::SliceAdaptor',
      'SupportingFeature'    =>
          'Bio::EnsEMBL::DBSQL::SupportingFeatureAdaptor',
      'Transcript'           => 'Bio::EnsEMBL::DBSQL::TranscriptAdaptor',
      'Translation'          => 'Bio::EnsEMBL::DBSQL::TranslationAdaptor',
      'AssemblyExceptionFeature' => 'Bio::EnsEMBL::DBSQL::AssemblyExceptionFeatureAdaptor'
    };


	# initialise storage for hash of names of current modules
	%{$self->_obj->{'current_module'}} = %{$self->_obj->{'default_module'}};

	# keep a hash of objects representing objects of each adaptor type
	# instantiated as required in get adaptor
	$self->_obj->{'current_objects'} = {};

	# initialise generic feature adaptor storage
	$self->_obj->{'generic_feature_adaptors'} = {};

  return $self;
}



=head2 get_ArchiveStableIdAdaptor

  Args       : none 
  Example    : none
  Description: ...
  Returntype : Bio::EnsEMBL::DBSQL::ArchiveStableIdAdaptor
  Exceptions : none
  Caller     : general

=cut

sub get_ArchiveStableIdAdaptor {
    my( $self ) = @_;

    return 
      $self->get_adaptor("ArchiveStableId");
}


=head2 get_QtlFeatureAdaptor

  Args       : none 
  Example    : none
  Description: ...
  Returntype : Bio::EnsEMBL::Map::DBSQL::QtlFeatureAdaptor
  Exceptions : none
  Caller     : general

=cut

sub get_QtlFeatureAdaptor {
  my $self = shift;
  return $self->get_adaptor("QtlFeature");
}

=head2 get_QtlAdaptor

  Args       : none 
  Example    : none
  Description: ...
  Returntype : Bio::EnsEMBL::Map::DBSQL::QtlAdaptor
  Exceptions : none
  Caller     : general

=cut

sub get_QtlAdaptor {
    my $self  = shift;
    return $self->get_adaptor("Qtl");
}

=head2 get_MetaContainer

  Args       : none
  Example    : $meta_container = $db_adaptor->get_MetaContainer(); 
  Description: Gets a MetaContainer object for this database
  Returntype : Bio::EnsEMBL::DBSQL::MetaContainer
  Exceptions : none
  Caller     : general

=cut

sub get_MetaContainer {
    my $self = shift;
    return $self->get_adaptor('MetaContainer');
}


=head2 get_ProteinFeatureAdaptor

  Args       : none 
  Example    : $pfa = $database_adaptor->get_ProteinFeatureAdaptor();  
  Description: Gets a ProteinFeatureAdaptor for this database.  
               Formerly named get_Protfeat_Adaptor()
  Returntype : Bio::EnsEMBL::DBSQL::ProteinFeatureAdaptor 
  Exceptions : none
  Caller     : general

=cut

sub get_ProteinFeatureAdaptor {
    my $self = shift;
    return $self->get_adaptor("ProteinFeature");
}


=head2 get_SNPAdaptor

  Args       : none 
  Example    : $snp_adaptor = $db_adaptor->get_SNPAdaptor();
  Description: Gets a SNPAdaptor for this database
  Returntype : Bio::EnsEMBL::DBSQL::SNPAdaptor
  Exceptions : none
  Caller     : general

=cut

sub get_SNPAdaptor {
  my ($self)  = @_;

  my $lite = $self->get_db_adaptor('lite');
  my $primary_adaptor;

  if($lite) {
    $primary_adaptor = $lite->get_SNPAdaptor();
  } else {
    my $snp = $self->get_db_adaptor('SNP');

    unless($snp) {
      warning("No lite or SNP database, cannot get snp adaptor\n");
      return undef;
    }

    $primary_adaptor = $snp->get_SNPAdaptor();
    $primary_adaptor->ensembl_db( $self );
  }

  #return a proxy adaptor which can use the lite or the core database
  return $self->get_adaptor("ProxySNP",
                            $primary_adaptor);
}


=head2 get_BlastAdaptor

  Args       : none 
  Example    : $blast_adaptor = $db_adaptor->get_BlastAdaptor();
  Description: Gets a BlastAdaptor for retrieving stored blast hits
  Returntype : Bio::EnsEMBL::External::BlastAdaptor
  Exceptions : none
  Caller     : general

=cut

sub get_BlastAdaptor {
  my ($self)  = @_;

  my $db_apt = $self->get_db_adaptor('blast');

  return $self->get_adaptor("Blast",
			     $db_apt);
}


=head2 get_PredictionTranscriptAdaptor

  Args       : none 
  Example    : $pta = $db_adaptor->get_PredictionTranscriptAdaptor();
  Description: Gets a PredictionTranscriptAdaptor for this database
  Returntype : Bio::EnsEMBL::DBSQL::PredictionTranscriptAdaptor
  Exceptions : none
  Caller     : general

=cut

sub get_PredictionTranscriptAdaptor {
  my ($self) = @_;

  return $self->get_adaptor("PredictionTranscript");
}


=head2 get_PredictionExonAdaptor

  Args       : none
  Example    : $pea = $db_adaptor->get_PredictionExonAdaptor();
  Description: Gets a PredictionExonAdaptor for this database
  Returntype : Bio::EnsEMBL::DBSQL::PredictionExonAdaptor
  Exceptions : none
  Caller     : general

=cut

sub get_PredictionExonAdaptor {
  my ($self) = @_;

  return $self->get_adaptor("PredictionExon");
}


=head2 get_SequenceAdaptor

  Args       : none 
  Example    : $sequence_adaptor = $db_adaptor->get_SequenceAdaptor();
  Description: Gets a SequenceAdaptor for this database
  Returntype : Bio::EnsEMBL::DBSQL::SequenceAdaptor
  Exceptions : none
  Caller     : general

=cut

sub get_SequenceAdaptor {
   my $self = shift;

   #return the sequence adaptor for the dnadb (which may be this db)
   return $self->dnadb->get_adaptor("Sequence");
}


=head2 get_GeneAdaptor

  Args       : none 
  Example    : $gene_adaptor = $db_adaptor->get_GeneAdaptor();
  Description: Gets a GeneAdaptor for this database
  Returntype : Bio::EnsEMBL::DBSQL::GeneAdaptor
  Exceptions : none
  Caller     : general

=cut

sub get_GeneAdaptor {
  my( $self ) = @_;
  #get a core db adaptor
  my $core_adaptor = $self->get_adaptor("Gene");
  
  #use a proxy gene adaptor, capable of making decisions with regards to the
  #database that it uses, passing in the core adaptor as a constructor arg
  return $self->get_adaptor("ProxyGene",
			     $core_adaptor);
}


=head2 get_ExonAdaptor

  Args       : none 
  Example    : $exon_adaptor = $db_adaptor->get_ExonAdaptor();
  Description: Gets an ExonAdaptor for this database
  Returntype : Bio::EnsEMBL::DBSQL::ExonAdaptor
  Exceptions : none
  Caller     : general

=cut

sub get_ExonAdaptor {
  my( $self ) = @_;
  
  return $self->get_adaptor("Exon");
}


=head2 get_TranscriptAdaptor

  Args       : none 
  Example    : $transcript_adaptor = $db_adaptor->get_TranscriptAdaptor();
  Description: Gets a TranscriptAdaptor for this database
  Returntype : Bio::EnsEMBL::DBSQL::TranscriptAdaptor
  Exceptions : none
  Caller     : general

=cut

sub get_TranscriptAdaptor {
  my( $self ) = @_;

  return $self->get_adaptor("Transcript");
}


=head2 get_TranslationAdaptor

  Args       : none 
  Example    : $translation_adaptor = $db_adaptor->get_TranslationAdaptor();
  Description: Gets a TranslationAdaptor for this database
  Returntype : Bio::EnsEMBL::DBSQL::TranslationAdaptor
  Exceptions : none
  Caller     : general

=cut

sub get_TranslationAdaptor {
    my( $self ) = @_;

    return $self->get_adaptor("Translation");
}



=head2 get_SliceAdaptor

  Args       : none 
  Example    : $slice_adaptor = $db_adaptor->get_SliceAdaptor();
  Description: Gets a SliceAdaptor for this database
  Returntype : Bio::EnsEMBL::DBSQL::SliceAdaptor
  Exceptions : none
  Caller     : general

=cut

sub get_SliceAdaptor {
  my( $self ) = @_;
  
  return $self->get_adaptor("Slice");
}


=head2 get_AnalysisAdaptor

  Args       : none 
  Example    : $analysis_adaptor = $db_adaptor->get_AnalysisAdaptor();
  Description: Gets an AnalysisAdaptor for this database
  Returntype : Bio::EnsEMBL::DBSQL::AnalysisAdaptor
  Exceptions : none
  Caller     : general

=cut

sub get_AnalysisAdaptor {
    my( $self ) = @_;

    return $self->get_adaptor("Analysis");
}


=head2 get_SimpleFeatureAdaptor

  Args       : none 
  Example    : $sfa = $db_adaptor->get_SimpleFeatureAdaptor();
  Description: Gets a SimpleFeatureAdaptor for this database
  Returntype : Bio::EnsEMBL::DBSQL::SimpleFeatureAdaptor
  Exceptions : none
  Caller     : general

=cut

sub get_SimpleFeatureAdaptor {
  my( $self ) = @_;
  
  return $self->get_adaptor("SimpleFeature");
}


=head2 get_RepeatConsensusAdaptor

  Args       : none 
  Example    : $rca = $db_adaptor->get_RepeatConsensusAdaptor();
  Description: Gets a RepeatConsensusAdaptor for this database
  Returntype : Bio::EnsEMBL::DBSQL::RepeatConsensusAdaptor
  Exceptions : none
  Caller     : general

=cut

sub get_RepeatConsensusAdaptor {
  my( $self ) = @_;
  
  return $self->get_adaptor("RepeatConsensus");
}


=head2 get_RepeatFeatureAdaptor

  Args       : none 
  Example    : $rfa = $db_adaptor->get_RepeatFeatureAdaptor();
  Description: Gets a RepeatFeatureAdaptor for this database
  Returntype : Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptor
  Exceptions : none
  Caller     : general

=cut

sub get_RepeatFeatureAdaptor {
  my( $self ) = @_;
  
  my $core_adaptor = 
    $self->dnadb->get_adaptor("RepeatFeature");
  
  #create a proxy adaptor, using a core RepeatFeatureAdaptor as constructor arg
  return $self->get_adaptor("ProxyRepeatFeature",
			    $core_adaptor);
}


=head2 get_ProteinAlignFeatureAdaptor

  Args       : none 
  Example    : $pafa = $db_adaptor->get_ProteinAlignFeatureAdaptor();
  Description: Gets a ProteinAlignFeatureAdaptor for this database
  Returntype : Bio::EnsEMBL::DBSQL::ProteinAlignFeatureAdaptor
  Exceptions : none
  Caller     : general

=cut

sub get_ProteinAlignFeatureAdaptor {
  my( $self ) = @_;
    
  return 
    $self->get_adaptor("ProteinAlignFeature");
}

  
=head2 get_DnaAlignFeatureAdaptor

  Args       : none 
  Example    : $dafa = $db_adaptor->get_DnaAlignFeatureAdaptor();
  Description: Gets a DnaAlignFeatureAdaptor for this database
  Returntype : Bio::EnsEMBL::DBSQL::DnaAlignFeatureAdaptor
  Exceptions : none
  Caller     : general

=cut

sub get_DnaAlignFeatureAdaptor {
  my $self = shift;

  return $self->get_adaptor("DnaAlignFeature");
}


=head2 get_AssemblyMapperAdaptor

  Args       : none 
  Example    : $asma = $db_adaptor->get_AssemblyMapperAdaptor();
  Description: Gets an AsemblyMapperAdaptor for this database
  Returntype : Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor
  Exceptions : none
  Caller     : general

=cut

sub get_AssemblyMapperAdaptor {
  my( $self ) = @_;

  return 
    $self->dnadb->get_adaptor("AssemblyMapper");
}


=head2 get_DBEntryAdaptor

  Args       : none 
  Example    : $dbentry_adaptor = $db_adaptor->get_DBEntryAdaptor();
  Description: Gets a DBEntryAdaptor for this database
  Returntype : Bio::EnsEMBL::DBSQL::DBEntryAdaptor
  Exceptions : none
  Caller     : general

=cut

sub get_DBEntryAdaptor {
    my( $self ) = @_;

    return $self->get_adaptor("DBEntry");
}



=head2 get_KaryotypeBandAdaptor

  Args       : none 
  Example    : $kba = $db_adaptor->get_KaryotypeBandAdaptor();
  Description: Gets a KaryotypeBandAdaptor for this database
  Returntype : Bio::EnsEMBL::DBSQL::KaryotypeBandAdaptor
  Exceptions : none
  Caller     : general

=cut

sub get_KaryotypeBandAdaptor {
    my( $self ) = @_;

    return 
      $self->dnadb->get_adaptor("KaryotypeBand");
}



=head2 get_SupportingFeatureAdaptor

  Arg [1]    : none
  Example    : $sfa = $db_adaptor->get_SupportingFeatureAdaptor();
  Description: Gets a SupportingFeatreAdaptor for this database
  Returntype : Bio::EnsEMBL::DBSQL::SupportingFeatureAdaptor
  Exceptions : none
  Caller     : Bio::EnsEMBL::Exon

=cut

sub get_SupportingFeatureAdaptor {
  my $self = shift;

  return $self->get_adaptor("SupportingFeature");
}



=head2 get_MarkerFeatureAdaptor

  Arg [1]    : none
  Example    : $mfa = $db_adaptor->get_MarkerFeatureAdaptor;
  Description: Gets a MarkerFeatureAdaptor for this database
  Returntype : Bio::EnsEMBL::Map::DBSQL::MarkerFeatureAdaptor
  Exceptions : none
  Caller     : general

=cut

sub get_MarkerFeatureAdaptor {
  my $self = shift;

  return $self->get_adaptor('MarkerFeature');
}


=head2 get_MarkerAdaptor

  Arg [1]    : none
  Example    : $ma = $db_adaptor->get_MarkerAdaptor;
  Description: Gets a MarkerAdaptor for this database
  Returntype : Bio::EnsEMBL::Map::DBSQL::MarkerAdaptor
  Exceptions : none
  Caller     : general

=cut

sub get_MarkerAdaptor {
  my $self = shift;

  return $self->get_adaptor('Marker');
}



=head2 get_CoordSystemAdaptor

  Arg [1]    : none
  Example    : $csa = $db_adaptor->get_CoordSystemAdaptor();
  Description: Gets a CoordSystemAdaptor for this database
  Returntype : Bio::EnsEMBL::DBSQL::CoordSystemAdaptor
  Exceptions : none
  Caller     : general

=cut

sub get_CoordSystemAdaptor {
  my $self = shift;

  return $self->get_adaptor('CoordSystem');
}




=head2 get_MiscSetAdaptor

  Arg [1]    : none
  Example    : $msa = $db_adaptor->get_MiscSetAdaptor();
  Description: Gets a MiscSet adaptor for this database
  Returntype : Bio::EnsEMBL::DBSQL::MiscSetAdaptor
  Exceptions : none
  Caller     : general

=cut

sub get_MiscSetAdaptor {
  my $self = shift;

  return $self->get_adaptor('MiscSet');
}



=head2 get_MiscFeatureAdaptor

  Arg [1]    : none
  Example    : $mfa = $db_adaptor->get_MiscFeatureAdaptor();
  Description: Gets a MiscFeature adaptor for this database
  Returntype : Bio::EnsEMBL::DBSQL::MiscFeatureAdaptor
  Exceptions : none
  Caller     : general

=cut

sub get_MiscFeatureAdaptor {
  my $self = shift;

  return $self->get_adaptor('MiscFeature');
}


=head2 get_AssemblyExceptionFeatureAdaptor

  Arg [1]    : none
  Example    : $aefa = $db_adaptor->get_AssebmyExceptionFeatureAdaptor();
  Description: Gets a AssemblyExceptionFeature adaptor for this database
  Returntype : Bio::EnsEMBL::DBSQL::AssemblyExceptionFeatureAdaptor
  Exceptions : none
  Caller     : general

=cut

sub get_AssemblyExceptionFeatureAdaptor {
  my $self = shift;

  return $self->get_adaptor('AssemblyExceptionFeature');
}


=head2 dnadb

 Title   : dnadb
 Usage   : my $dnadb = $db->dnadb;
 Function: returns the database adaptor where the dna lives
           Useful if you only want to keep one copy of the dna
           on disk but have other databases with genes and features in
 Returns : dna database adaptor
  Args    : Bio::EnsEMBL::DBSQL::DBAdaptor

=cut

sub dnadb {
  my ($self,$arg) = @_;

  if($arg) {
    if(ref $arg && (($arg->isa('Bio::EnsEMBL::Container') && 
		     $arg->_obj == $self) || $arg == $self)) {
      #we don't want to store a circular reference to our self
      return;
    }

    #if this a container, we don't care, hang onto it
    $self->{'dnadb'} = $arg;
  }

  return $self->{'dnadb'} || $self;
}



=head2 deleteObj

  Arg [1]    : none
  Example    : none
  Description: Cleans up circular reference loops so proper garbage collection
               can occur.
  Returntype : none
  Exceptions : none
  Caller     : DBAdaptorContainer::DESTROY

=cut

sub deleteObj {
  my $self = shift;
  #print "called deleteObj on DBAdaptor\n";

  #clean up external feature adaptor references
  if(exists $self->{'_xf_adaptors'}) {
    foreach my $key (keys %{$self->{'_xf_adaptors'}}) {
      delete $self->{'_xf_adaptors'}->{$key};
    }
  }

  if(exists $self->{'current_objects'}) {
    foreach my $adaptor_name (keys %{$self->{'current_objects'}}) {
      my $adaptor = $self->{'current_objects'}->{$adaptor_name};

      if($adaptor && $adaptor->can('deleteObj')) {
        $adaptor->deleteObj();
      }

      delete $self->{'current_objects'}->{$adaptor_name};
    }
  }

  delete $self->{'_meta_container'};
  delete $self->{'dnadb'};

  #call the superclass deleteObj method
  $self->SUPER::deleteObj;
}


###########################################################
#
# Support for DAS
#
###########################################################

=head2 add_DASFeatureFactory

  Arg [1]    : Bio::EnsEMBL::ExternalFeatureFactory $value 
  Example    : none
  Description: Attaches a DAS Feature Factory to this method.  
               ExternalFeatureFactory objects are not really used right now.
               They may be reintroduced or taken out completely.  The fate
               of this function is unknown (although it is presently needed).
  Returntype : none
  Exceptions : none
  Caller     : EnsWeb

=cut

sub add_DASFeatureFactory{
 
 my ($self,$value) = @_;
  
  push(@{$self->{'_das_ff'}},$value);
}


=head2 _each_DASFeatureFactory

  Args       : none
  Example    : none
  Description: Not sure if this is used, or if it should be removed.  It 
               does not seem to be used at the moment
  Returntype : Bio::EnsEMBL::ExternalFeatureFactory
  Exceptions : none
  Caller     : ??

=cut

sub _each_DASFeatureFactory{
   my ($self) = @_;

   return @{$self->{'_das_ff'}||[]}
}


################################################################## 
# 
# SUPPORT FOR EXTERNAL FEATURE FACTORIES 
# 
##################################################################



=head2 add_ExternalFeatureAdaptor

  Arg [1]    : Bio::EnsEMBL::External::ExternalFeatureAdaptor
  Example    : $db_adaptor->add_ExternalFeatureAdaptor($xfa);
  Description: Adds an external feature adaptor to this database adaptor.
               Adding the external adaptor in this way allows external
               features to be obtained from Slices and from RawContigs.

               The external feature adaptor which is passed to this method
               will have its db attribuite set to this DBAdaptor object via 
               the db accessor method. 

               ExternalFeatureAdaptors passed to this method are stored 
               internally in a hash keyed on the string returned by the 
               ExternalFeatureAdaptors track_name method.
               
               If the track name method is not implemented then the 
               a default key named 'External features' is assigned.  In the
               event of duplicate key names, a number is appended to the
               key name, and incremented for each subsequent adaptor with the
               same track name.  For example, if no track_names are specified 
               then the the external feature adaptors will be stored under the
               keys 'External features', 'External features2' 
               'External features3' etc.
  Returntype : none
  Exceptions : none
  Caller     : general

=cut

sub add_ExternalFeatureAdaptor {
  my ($self, $adaptor) = @_;

  unless($adaptor && ref $adaptor && 
	 $adaptor->isa('Bio::EnsEMBL::External::ExternalFeatureAdaptor')) {
     throw("[$adaptor] is not a " .
           "Bio::EnsEMBL::External::ExternalFeatureAdaptor");
  }

  unless(exists $self->{'_xf_adaptors'}) {
    $self->{'_xf_adaptors'} = {};
  }

  my $track_name = $adaptor->{'_track_name'};

  #use a generic track name if one hasn't been defined
  unless(defined $track_name) {
    $track_name = "External features";
  }

  #if the track name exists add numbers to the end until a free name is found
  if(exists $self->{'_xf_adaptors'}->{"$track_name"}) {
    my $num = 2;
    $num++ while(exists $self->{'_xf_adaptors'}->{"$track_name$num"});
    $self->{'_xf_adaptors'}->{"$track_name$num"} = $adaptor;
  } else {
    $self->{'_xf_adaptors'}->{"$track_name"} = $adaptor;
  }

  $adaptor->db($self);
}



=head2 get_ExternalFeatureAdaptors

  Arg [1]    : none
  Example    : @xfas = values %{$db_adaptor->get_ExternalFeatureAdaptors}; 
  Description: Retrieves all of the ExternalFeatureAdaptors which have been
               added to this DBAdaptor.  The ExternalFeatureAdaptors are 
               returned in a reference to a hash keyed on the track names
               of the external adaptors
  Returntype : Reference to a hash of ExternalFeatureAdaptors keyed on 
               their track names.
  Exceptions : none
  Caller     : general

=cut

sub get_ExternalFeatureAdaptors {
  my $self = shift;

  return $self->{'_xf_adaptors'};
}


=head2 add_ExternalFeatureFactory

  Arg [1]    : Bio::EnsEMBL::DB::ExternalFeatureFactoryI $value
  Example    : $db_adaptor->add_ExternalFeatureFactory
  Description: It is recommended that add_ExternalFeatureAdaptor be used 
               instead.  See documentation for 
               Bio::EnsEMBL::External::ExternalFeatureAdaptor

               Adds an external feature factory to the core database
               so that features from external sources can be displayed in 
               ensembl. This method is still available mainly for legacy
               support for external EnsEMBL installations.
  Returntype : none
  Exceptions : none
  Caller     : external

=cut

sub add_ExternalFeatureFactory{
   my ($self,$value) = @_;

   $self->add_ExternalFeatureAdaptor($value);
}

#
# OVERWRITABLE STANDARD ADAPTORS
#

=head2 get_adaptor

  Arg [1]    : Canonical data type for which an adaptor is required.
  Example    : $db_adaptor->get_adaptor("Protein")
  Description: Gets an adaptor object for a standard data type.
  Returntype : Adaptor Object of arbitrary type
  Exceptions : thrown if there is no associated module
  Caller     : external

=cut

sub get_adaptor() {
	my ($self, $canonical_name, @other_args) = @_;

  if ($self->isa('Bio::EnsEMBL::Container')) {
    $self = $self->_obj;
  }

	# throw if module for $canonical_name does not exist
	throw("No such data type $canonical_name") 
    if (!exists($self->{'current_module'}->{$canonical_name}));

	# get module name for $canonical_name
	my $module_name = $self->{'default_module'}->{$canonical_name};

	# create and store a new one if necessary
	if (!exists($self->{'current_objects'}->{$canonical_name})) {
	  $self->{'current_objects'}->{$canonical_name} =
      $self->_get_adaptor($module_name, @other_args);
	}

	return $self->{'current_objects'}->{$canonical_name};

}


=head2 set_adaptor

  Arg [1]    : Canonical data type for new adaptor.
	Arg [2]    : Object defining the adaptor for arg1.
  Example    : $pa = Bio::EnsEMBL::DBSQL::ProteinAdaptor->new($db_adaptor);
             : $db_adaptor->set_adaptor("Protein", $pa)
  Description: Stores the object which represents the adaptor for the
               arg1 data type.
  Returntype : none
  Exceptions : If arg2 is not a subclass of the default module for this
               data type.
  Caller     : external

=cut

sub set_adaptor() {
	my ($self, $canonical_name, $new_object) = @_;

  if ($self->isa('Bio::EnsEMBL::Container')) {
    $self = $self->_obj;
  }

  # throw if an unrecognised canonical_name is used
	throw("No such data type $canonical_name") 
    if(!exists($self->{'default_module'}->{$canonical_name}));

	my $default_module = $self->{'default_module'}->{$canonical_name};
	
	# Check that $new_module is a subclass of $default_module	
	if (!$new_object->isa($default_module)) {  # polymorphism should work
		throw("ref($new_object) is not a subclass of $default_module");
	}

	# set the value in current_module
	$self->{'current_objects'}->{$canonical_name} = $new_object;
}

#
# GENERIC FEATURE ADAPTORS
#

=head2 get_GenericFeatureAdaptors

  Arg [1]    : List of names of feature adaptors to get. If no
               adaptor names are given, all the defined adaptors are returned.
  Example    : $db->get_GenericFeature("SomeFeature", "SomeOtherFeature")
  Description: Returns a hash containing the named feature adaptors (or
               all feature adaptors).
  Returntype : reference to a Hash containing the named
               feature adaptors (or all feature adaptors).
  Exceptions : If any of the the named generic feature adaptors do not exist.
  Caller     : external

=cut

sub get_GenericFeatureAdaptors() {

  my ($self, @names) = @_;

  my %adaptors = ();

  if (!@names) {
    %adaptors = %{$self->{'generic_feature_adaptors'}};
  } else {
    foreach my $name (@names) {
      if (!exists($self->{'generic_feature_adaptors'}->{$name})) {
        throw("No generic feature adaptor has been defined for $name" );
      }
      $adaptors{$name} = $self->{'generic_feature_adaptors'}->{$name};
    }
  }

  return \%adaptors;
}


=head2 add_GenericFeatureAdaptor

  Arg [1]    : The name of the feature.
  Arg [2]    : Adaptor object for a generic feature.
  Example    : $db->add_GenericFeatureAdaptor("SomeFeature",
                              "Bio::EnsEMBL::DBSQL::SomeFeatureAdaptor")
  Description: Stores the object which represents the adaptor for the
               named feature type.
  Returntype : none
  Exceptions :
  Caller     : external

=cut

sub add_GenericFeatureAdaptor() {

	my ($self, $name, $adaptor_obj) = @_;
	
	# check that $adaptor is an object that subclasses BaseFeatureAdaptor	
	if (!$adaptor_obj->isa("Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor")) {
	  throw("$name is a " . ref($adaptor_obj) . "which is not a " .
          "subclass of Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor" );
	}

	$self->{'generic_feature_adaptors'}->{$name} = $adaptor_obj;
}




#########################
# sub DEPRECATED METHODS
#########################



=head2 get_MapFragAdaptor

  Description: MapFragAdaptor is deprecated. Use MiscFeatureAdaptor instead.

=cut

sub get_MapFragAdaptor {
  my $self = shift;

  return $self->get_adaptor( "MapFrag" );
}

=head2 get_CloneAdaptor

  Description: CloneAdaptor is deprecated.  Use SliceAdaptor instead.

=cut

sub get_CloneAdaptor {
  my( $self ) = @_;

  return $self->dnadb->get_adaptor("Clone");
}

=head2 get_ChromosomeAdaptor

  Description: ChromosomeAdaptor is deprecated.  Use SliceAdaptor instead.

=cut

sub get_ChromosomeAdaptor {
    my( $self ) = @_;

    return 
      $self->dnadb->get_adaptor("Chromosome");
}

=head2 get_RawContigAdaptor

  Description: RawContigAdaptor is deprecated. Use SliceAdaptor instead.

=cut

sub get_RawContigAdaptor {
    my( $self ) = @_;

    return $self->dnadb->get_adaptor("RawContig");
}


=head2 get_ProteinAdaptor

  Description: ProteinAdaptor is deprecated. Use TranslationAdaptor instead

=cut

sub get_ProteinAdaptor {
    my $self  = shift;
    deprecate("The ProteinAdaptor is deprecated. Use the TranslationAdaptor " .
              "instead of the ProteinAdaptor and Translation instead of " .
              "Protein.");
    return $self->get_adaptor("Protein");
}


sub source {
  deprecate('Do not use - this method does nothing');
}


=head2 assembly_type

  Description: DEPRECATED - Use CoordSystemAdaptor to obtain default coordinate
               system instead.

=cut

sub assembly_type{
  my $self = shift;

  deprecate('Use CoordSystemAdaptor::fetch_top_level instead');

  my $csa = $self->get_CoordSystemAdaptor();

  #get the default top-level coord system
  my ($dbID,$name,$version) = $csa->fetch_top_level();

  return $version;
}



=head2 list_supported_assemblies

  Description: DEPRECATED - Use CoordSystemAdaptor to obtain list of top-level
               coordinate systems instead

=cut

sub list_supported_assemblies {
  my($self) = @_;
  deprecate('Use CoordSystemAdaptor::fetch_all_top_level instead');

  my $csa = $self->get_CoordSystemAdaptor();
  my @versions;
  foreach my $cs (@{$csa->fetch_all_top_level()}) {
    my ($dbID, $name, $version) = @$csa;
    push @versions, $version;
  }
  return @versions;
}


1;
