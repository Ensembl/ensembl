
=Head1 NAME - Bio::EnsEMBL::DBSQL::DBAdaptor

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

This object represents a database that is implemented somehow (you shouldn\'t
care much as long as you can get the object). Once created you can retrieve
database adaptors specific to various database objects that allow the
retrieval and creation of objects from the database,

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::DBSQL::DBAdaptor;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DBSQL::DBConnection;

@ISA = qw(Bio::EnsEMBL::DBSQL::DBConnection);

#Override constructor inherited by Bio::EnsEMBL::DBSQL::DBConnection


=head2 new

  Arg [1]    : string SOURCE 
               The source name of the database.  This may be removed soon.
  Arg [2]    : Bio::EnsEMBL::DBSQL::DBAdaptor DNADB 
               The dna database to be attached to this database.  This may also
               be changed.
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
  
  my ( $source, $dnadb ) = $self->_rearrange([qw(SOURCE DNADB)],@args);  

  if(defined $source) {
    $self->source($source);
  }

  if(defined $dnadb) {
    $self->dnadb($dnadb);
  }
 
  return $self; # success - we hope!
}


=head2 source

  Arg [1]    : (optional) string $source
               The source of info in the database connected to (e.g. 'embl') 
  Example    : $db_adaptor->source('sanger');
  Description: Sets/Gets the source or human readable name of the genes in 
               the connected database. For example for the sanger db the source
               would be 'sanger'. 
  Returntype : string
  Exceptions : none 
  Caller     : Bio::EnsEMBL::GeneAdaptor  Bio::EnsEMBL::LiteGeneAdaptor EnsWeb 

=cut

sub source {
  my ($self, $source) = @_;

  if(defined $source) {
    $self->{'_source'} = $source;
  }

  return $self->{'_source'};
}


=head2 get_MetaContainer

  Args       : none
  Example    : $meta_container = $db_adaptor->get_MetaContainer(); 
  Description: Gets a MetaContainer object for this database
  Returntype : Bio::EnsEMBL::DBSQL::MetaContainer
  Exceptions : none
  Caller     : ?

=cut

sub get_MetaContainer {
    my( $self ) = @_;
    
    my( $mc );
    unless ($mc = $self->{'_meta_container'}) {
        require Bio::EnsEMBL::DBSQL::MetaContainer;
        $mc = Bio::EnsEMBL::DBSQL::MetaContainer->new($self);
        $self->{'_meta_container'} = $mc;
    }
    return $mc;
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
    my( $self ) = @_;
    
    return 
      $self->_get_adaptor("Bio::EnsEMBL::DBSQL::ProteinFeatureAdaptor");
}


=head2 get_ProteinAdaptor

  Args       : none 
  Example    : $pa = $database_adaptor->get_ProteinAdaptor();  
  Description: Gets a ProteinAdaptor for this database.  
               Formerly named get_Protein_Adaptor() 
  Returntype : Bio::EnsEMBL::DBSQL::ProteinAdaptor 
  Exceptions : none 
  Caller     : general

=cut

sub get_ProteinAdaptor {
    my( $self ) = @_;
 
    return $self->_get_adaptor("Bio::EnsEMBL::DBSQL::ProteinAdaptor");
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

  #return a proxy adaptor which can use the lite or the core database
  return $self->_get_adaptor("Bio::EnsEMBL::DBSQL::ProxySNPAdaptor");
}


=head2 get_MapFragAdaptor

  Args       : none 
  Example    : $map_frag_adaptor = $db_adaptor->get_MapFragAdaptor();
  Description: Gets a MapFragAdaptor for this database
  Returntype : Bio::EnsEMBL::DBSQL::MapFragAdaptor
  Exceptions : none
  Caller     : general

=cut

sub get_MapFragAdaptor {
  my $self = shift;

  return $self->_get_adaptor( "Bio::EnsEMBL::DBSQL::MapFragAdaptor" );
}


=head2 get_CloneAdaptor

  Args       : none 
  Example    : $clone_adaptor = $db_adaptor->get_CloneAdaptor();
  Description: Gets a CloneAdaptor for this database
  Returntype : Bio::EnsEMBL::DBSQL::CloneAdaptor
  Exceptions : none
  Caller     : general

=cut

sub get_CloneAdaptor {
  my( $self ) = @_;
  
  return $self->_get_adaptor("Bio::EnsEMBL::DBSQL::CloneAdaptor");
}


=head2 get_LandmarkMarkerAdaptor

  Args       : none 
  Example    : $marker_adaptor = $db_adaptor->get_LandmarkMarkerAdaptor();
  Description: Gets a LandmarkMarkerAdaptor which currently is always retrieved
               from the lite database attached to this database.
  Returntype : Bio::EnsEMBL::DBSQL::LandmarkMarkerAdaptor
  Exceptions : none
  Caller     : general

=cut

sub get_LandmarkMarkerAdaptor {
  my $self = shift;

  if( defined $self->lite_DBAdaptor() ) {
    return $self->lite_DBAdaptor()->get_LandmarkMarkerAdaptor();
  } else {
    return undef;
  }
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

  return 
    $self->_get_adaptor("Bio::EnsEMBL::DBSQL::PredictionTranscriptAdaptor");
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

   return $self->_get_adaptor("Bio::EnsEMBL::DBSQL::SequenceAdaptor");
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
  my $core_adaptor = $self->_get_adaptor("Bio::EnsEMBL::DBSQL::GeneAdaptor");
  
  #use a proxy gene adaptor, capable of making decisions with regards to the
  #database that it uses, passing in the core adaptor as a constructor arg
  return $self->_get_adaptor("Bio::EnsEMBL::DBSQL::ProxyGeneAdaptor",
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
  
  return $self->_get_adaptor("Bio::EnsEMBL::DBSQL::ExonAdaptor");
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
  
  return $self->_get_adaptor("Bio::EnsEMBL::DBSQL::TranscriptAdaptor");
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

    return $self->_get_adaptor("Bio::EnsEMBL::DBSQL::TranslationAdaptor");
}


=head2 get_FeatureAdaptor

  Args       : none 
  Example    : $feature_adaptor = $db_adaptor->get_FeatureAdaptor();
  Description: Gets a FeatureAdaptor for this database
  Returntype : Bio::EnsEMBL::DBSQL::FeatureAdaptor
  Exceptions : none
  Caller     : general

=cut

sub get_FeatureAdaptor {
    my( $self ) = @_;

    return $self->_get_adaptor("Bio::EnsEMBL::DBSQL::FeatureAdaptor");
}


=head2 get_RawContigAdaptor

  Args       : none 
  Example    : $rca = $db_adaptor->get_RawContigAdaptor();
  Description: Gets a RawContigAdaptor for this database
  Returntype : Bio::EnsEMBL::DBSQL::RawContigAdaptor
  Exceptions : none
  Caller     : general

=cut

sub get_RawContigAdaptor {
    my( $self ) = @_;

    return $self->_get_adaptor("Bio::EnsEMBL::DBSQL::RawContigAdaptor");
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
  
  return $self->_get_adaptor("Bio::EnsEMBL::DBSQL::SliceAdaptor");
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

    return $self->_get_adaptor("Bio::EnsEMBL::DBSQL::AnalysisAdaptor");
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
  
  return $self->_get_adaptor("Bio::EnsEMBL::DBSQL::SimpleFeatureAdaptor");
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
  
  return $self->_get_adaptor("Bio::EnsEMBL::DBSQL::RepeatConsensusAdaptor");
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
    $self->_get_adaptor("Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptor");
  
  #create a proxy adaptor, using a core RepeatFeatureAdaptor as constructor arg
  return $self->_get_adaptor("Bio::EnsEMBL::DBSQL::ProxyRepeatFeatureAdaptor",
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
    $self->_get_adaptor("Bio::EnsEMBL::DBSQL::ProteinAlignFeatureAdaptor");
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
  
  my $core_adaptor = 
    $self->_get_adaptor("Bio::EnsEMBL::DBSQL::DnaAlignFeatureAdaptor");

  #return a proxy adaptor which can choose between the core and est DBs
  return 
    $self->_get_adaptor("Bio::EnsEMBL::DBSQL::ProxyDnaAlignFeatureAdaptor",
			$core_adaptor);
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
    
  return $self->_get_adaptor("Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor");
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

    return $self->_get_adaptor("Bio::EnsEMBL::DBSQL::DBEntryAdaptor");
}


=head2 get_StaticGoldenPathAdaptor

  Args       : none 
  Example    : $sgpa = $db_adaptor->get_StaticGoldenPathAdaptor();
  Description: Gets a StaticGoldenPathAdaptor for this database
               Use of the StaticGoldenPathAdaptor is not recommended.  It is
               being phased out, and is largly deprecated already.  The 
               SliceAdaptor or AssemblyMapperAdaptor may be a viable 
               alternatives.
  Returntype : Bio::EnsEMBL::DBSQL::AnalysisAdaptor
  Exceptions : none
  Caller     : general

=cut

sub get_StaticGoldenPathAdaptor{
  my( $self ) = @_;
  
  return $self->_get_adaptor("Bio::EnsEMBL::DBSQL::StaticGoldenPathAdaptor");
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

    return $self->_get_adaptor("Bio::EnsEMBL::DBSQL::KaryotypeBandAdaptor");
}


=head2 get_ChromosomeAdaptor

  Args       : none 
  Example    : $ca = $db_adaptor->get_ChromosomeAdaptor();
  Description: Gets a ChromosomeAdaptor for this database
  Returntype : Bio::EnsEMBL::DBSQL::ChromosomeAdaptor
  Exceptions : none
  Caller     : general

=cut

sub get_ChromosomeAdaptor {
    my( $self ) = @_;

    return $self->_get_adaptor("Bio::EnsEMBL::DBSQL::ChromosomeAdaptor");
}



=head2 list_supported_assemblies

  Args       : none
  Example    : my @supported = $database_adaptor->list_supported_assemblies);
  Description: Returns a list of assemblies supported by this database.
  Returntype : list of strings
  Exceptions : thrown if SQL query fails
  Caller     : ?

=cut

sub list_supported_assemblies {
    my($self) = @_;
    my @out;

    my $query = q{
        SELECT distinct type
        FROM   assembly
    };

    my $sth = $self->prepare($query) ||
     $self->throw("Error in list_supported_assemblies");
    my $res = $sth->execute ||
     $self->throw("Error in list_supported_assemblies");

    while (my($type) = $sth->fetchrow_array) {
       push(@out, $type);
    }
    return @out;
}

	
=head2 deleteObj

  Args       : none
  Example    : $db_adaptor->deleteObj();
  Description: Explicitly destroys this object and objects referenced by 
               this object.  This method should only be called if you know
               what you are doing, and is only needed for object destruction
               when circular references are present (these will prevent 
               perls automatic garbage collection).
  Returntype : none
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBSQL::DBAdaptor

=cut		     

sub deleteObj {

  my  $self=shift;
  my $dummy;
  $self->DESTROY;
  
  foreach my $name ( keys %{$self} ) {
    eval {
      $dummy = $self->{$name}; 
      $self->{$name}  = undef;
      $dummy->deleteObj;
    };
  }
}


=head2 assembly_type

  Arg [1]    : (optional) string $value
                the new assembly type value
  Example    : $db_adaptor->assembly_type($newval);
  Description: Getter / Setter for the type of assembly used by
               this database.  If the value is not set then the 
               default value is obtained from the MetaContainer
  Returntype : string
  Exceptions : thrown if there is no defined assembly type, and the default
               assembly type cannot be obtained from the MetaContainer 
  Caller     : ?

=cut

sub assembly_type{
   my ($obj,$value) = @_;
   if($value) {
      $obj->{'assembly'} = $value;
    }
    if (! defined $obj->{'assembly'}) {
      print STDERR "DBAdaptor.pm: using default assembly type\n";
      my $ass;
      eval {
        $ass = $obj->get_MetaContainer()->get_default_assembly();
      };
      if ( $@ ) {
        $obj->throw("*** get_MetaContainer->get_default_assembly failed:\n$@\n"
          ."assembly type must be set with assembly_type() first");
      } elsif (! $ass) {
        $obj->throw("No default assembly defined"
          . " - must set with assembly_type() first");
      }
      $obj->{'assembly'} = $ass;
    }
    return $obj->{'assembly'};

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

  if (defined($arg)) {
    if (! $arg->isa("Bio::EnsEMBL::DBSQL::DBAdaptor")) {
      $self->throw("[$arg] is not a Bio::EnsEMBL::DBSQL::DBAdaptor");
    }
    $self->{_dnadb} = $arg;
  }
  return $self->{_dnadb} || $self;
}


=head2 lite_DBAdaptor

  Arg [1]    : (optional) Bio::EnsEMBL::Lite::DBAdaptor $liteDBConnection
               A denormalized Lite database to attach to this database
  Example    : $lite_db = $db_adaptor->lite_DBAdaptor();
  Description: A Getter/Setter for the lite database adaptor attached to this
               database
  Returntype : Bio::EnsEMBL::Lite::DBAdaptor
  Exceptions : none
  Caller     : EnsWeb

=cut

sub lite_DBAdaptor {
  my ($self, $arg ) = @_;
  if ( defined $arg ) {
    $self->{_liteDB} = $arg;

  }

  return $self->{_liteDB};
}


=head2 SNP_DBAdaptor

  Arg [1]    : (optional) Bio::EnsEMBL::ExternalData::SNPSQL::DBAdaptor $arg
               An external SNP database to be attached to this database
  Example    : $db_adaptor->SNP_DBAdaptor($snp_database_adaptor);
  Description: A Getter/Setter for the external SNP database adaptor 
               attached to this database
  Returntype : Bio::EnsEMBL::ExternalData::SNPSQL::DBAdaptor
  Exceptions : none
  Caller     : EnsWeb

=cut

sub SNP_DBAdaptor {
  my ($self, $arg) = @_;

  if(defined $arg) {
    $self->{_SNP_db} = $arg;
  }

  return $self->{_SNP_db};
}


=head2 map_DBAdaptor

  Arg [1]    : (optional) Bio::EnsEMBL::Map::DBSQL::DBAdaptor $arg
               A Map database to attach to this database
  Example    : $db_adaptor->lite_db($map_db_adaptor);
  Description: A Getter/Setter for the Map database adaptor attached to this
               database
  Returntype : Bio::EnsEMBL::Map::DBSQL::DBAdaptor
  Exceptions : none
  Caller     : EnsWeb

=cut

sub map_DBAdaptor {
  my ($self, $arg ) = @_;
  if ( defined $arg ) {
    $self->{_mapDB} = $arg;
  }

  return $self->{_mapDB};
}


=head2 est_DBAdaptor

  Arg [1]    : (optional) Bio::EnsEMBL::ExternalData::ESTSQL::DBAdaptor $arg
               An external EST database to attach to this database
  Example      $db_adaptor->est_DBAdaptor();
  Description: A Getter/Setter for the EST database adaptor attached to this
               database
  Returntype : Bio::EnsEMBL::ExternalData::ESTSQL::DBAdator
  Exceptions : none
  Caller     : EnsWeb

=cut

sub est_DBAdaptor {
  my ($self, $arg) = @_;
  
  if(defined $arg) {
    $self->{_estDB} = $arg;
  }

  return $self->{_estDB};
}

=head2 add_DASFeatureFactory

  Arg [1]    : Bio::EnsEMBL::DB::ExternalFeatureFactoryI $value 
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
  
  unless( ref $value && 
	  $value->isa('Bio::EnsEMBL::DB::ExternalFeatureFactoryI') ) {
    $self->throw("[$value] is not a Bio::EnsEMBL::DB::ExternalFeatureFactoryI"
		 . " but it should be!");
  }
  
  push(@{$self->{'_das_ff'}},$value);
}



=head2 _each_DASFeatureFactory

  Args       : none
  Example    : none
  Description: Not sure if this is used, or if it should be removed.  It 
               does not seem to be used at the moment
  Returntype : Bio::EnsEMBL::DB::ExternalFeatureFactoryI
  Exceptions : none
  Caller     : ??

=cut

sub _each_DASFeatureFactory{
   my ($self) = @_;

   return @{$self->{'_das_ff'}}
}






################################################################## 
# 
# SUPPORT FOR EXTERNAL ADAPTORS & FEATURE FACTORIES 
# 
# These are not implemented on the new main trunk and at this
# point it is not clear if they ever will be.
#
##################################################################




=head2 extension_tables

  Arg [1]    : none
  Example    : none
  Description: NOT IMPLEMENTED
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub extension_tables{
   my $obj = shift;

   $obj->warn("Extension tables is not implemented, and will either be ". 
	      "deprecated or implemented in the near future\n");

   return undef;

#   if( @_ ) {
#      my $value = shift;
#      $obj->{'extension_tables'} = $value;
#    }
#    return $obj->{'extension_tables'};

}


=head2 list_ExternalAdaptors

  Arg [1]    : none
  Example    : none
  Description: NOT CURRENTLY USED
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub list_ExternalAdaptors {
    my ($self) = @_;

    $self->warn("DBAdaptor::list_ExternalAdaptors is not implmented. It will "
		. "either be implemented or deprecated at a later date\n");

    return undef;

    #return keys % {$self->{_ext_adaptors}};
}


=head2 add_ExternalAdaptor

  Arg [1]    : none
  Example    : none
  Description: NOT CURRENTLY USED
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub add_ExternalAdaptor {
    my ($self, $adtor_name, $adtor_obj) = @_;

    $self->warn("DBAdaptor::add_ExternalAdaptor is not implemented. It will "
		. "either be implemented or deprecated at a later date\n");

    #$self->_ext_adaptor($adtor_name, $adtor_obj);
    #undef;
}


=head2 get_ExternalAdaptor

  Arg [1]    : none
  Example    : none
  Description: NOT CURRENTLY USED
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_ExternalAdaptor {
    my ($self, $adtor_name) = @_;
    
    $self->warn("DBAdaptor::get_ExternalAdaptor is not implemented. It will "
		. "either be implemented or deprecated at a later date\n");

    return undef;

    #$self->_ext_adaptor($adtor_name);
}


=head2 remove_ExternalAdaptor

  Arg [1]    : none
  Example    : none
  Description: NOT CURRENTLY USED
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub remove_ExternalAdaptor {
    my ($self, $adtor_name) = @_;

    $self->warn("remove_ExternalAdaptor is not implemented.  It will either "
               . "be implemented, or deprecated at a later date\n");

    return undef;

    #$self->_ext_adaptor($adtor_name, 'DELETE');
    #undef;
}


=head2 add_ExternalFeatureFactory

  Arg [1]    : none
  Example    : none
  Description: NOT CURRENTLY USED
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub add_ExternalFeatureFactory{
   my ($self,$value) = @_;

   $self->warn("DBAdaptor::add_ExternalFeatureFactory is not yet implemented."
	    . " it will either be implemented or deprecated at a later date");

#   unless( ref $value && $value->isa('Bio::EnsEMBL::DB::ExternalFeatureFactoryI') ) {
#       $self->throw("[$value] is not a Bio::EnsEMBL::DB::ExternalFeatureFactoryI but it should be!");
#   }

#   push(@{$self->{'_external_ff'}},$value);
}


=head2 _each_ExternalFeatureFactory

  Arg [1]    : none
  Example    : none
  Description: NOT CURRENTLY USED
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub _each_ExternalFeatureFactory{
   my ($self) = @_;


   $self->warn("DBAdaptor::_each_ExternalFeatureFactory is not yet implemented"
	     . "it will either be implemented or deprecated at a later date");

   return undef;

   #return @{$self->{'_external_ff'}}
}


## internal stuff for external adaptors


=head2 _ext_adaptor

  Arg [1]    : none
  Example    : none
  Description: NOT CURRENTLY USED
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub _ext_adaptor {
    my ($self, $adtor_name, $adtor_obj) = @_;
    
    $self->throw("No adaptor name given") unless $adtor_name;
    
    if( $adtor_obj ) {
        if ($adtor_obj eq 'DELETE') { 
            delete $adtor_obj->{'_ext_adaptors'}{$adtor_name};
        } else {
            $self->{'_ext_adaptors'}{$adtor_name} = $adtor_obj;
        }
    }
    return $self->{'_ext_adaptors'}{$adtor_name};
}


=head2 add_db_adaptor

  Arg [1]    : none
  Example    : none
  Description: NOT CURRENTLY USED
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub add_db_adaptor {
  my ($self, $name, $adaptor) = @_;

  $self->{'_db_adaptors'}->{$name} = $adaptor;
}

=head2 remove_db_adaptor

  Arg [1]    : none
  Example    : none
  Description: NOT CURRENTLY USED
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub remove_db_adaptor {
  my ($self, $name) = @_;

  my $adaptor = $self->{'_db_adaptors'}->{$name};
  delete $self->{'_db_adaptors'}->{$name};

  return $adaptor;
}

=head2 get_all_db_adaptors

  Arg [1]    : none
  Example    : none
  Description: NOT CURRENTLY USED
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_all_db_adaptors {
  my ($self, $name) = @_;   

  return $self->{'_db_adaptors'};
}

=head2 get_db_adaptor

  Arg [1]    : none
  Example    : none
  Description: NOT CURRENTLY USED
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_db_adaptor {
  my ($self, $name) = @_;

  return $self->{'_db_adaptors'}->{$name};
}



##############################################################
###################DEPRECATED METHODS#########################
##                                                          ##
##  All the methods below are deprecated methods,           ##
##  only kept here to allow old scripts to work             ##
##  They all send a warning and call the new method instead ##
##                                                          ##
##############################################################
##############################################################


=head1 Old Deprecated Functions 

Functions which are completely deprecated 

=cut


=head2 get_LiteAdaptor

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use Bio::EnsEMBL::Lite::DBAdaptor instead or
               lite_DBAdaptor
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_LiteAdaptor {
    my( $self ) = @_;
    
    $self->throw("The lite adaptor is deprecated. Use the " .
		"Bio::EnsEMBL::Lite::DBAdaptor instead.\n" .
		 "This may be attached to the core DBAdaptor using the" .
		"lite_DBAdaptor method\n");

    return undef;
}


=head2 feature_Obj

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub feature_Obj {
    my $self = shift;
    $self->throw("No more Feature Objs!");

}


=head2 get_Gene

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED 
               use Bio::EnsEMBL::DBSQL::GenAdaptor::fetch_by_stable_id instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_Gene {
   my ($self,$geneid, $supporting) = @_;

   $self->warn("DBAdaptor::get_Gene is deprecated\n" .
	   "use Bio::EnsEMBL::DBSQL::GeneAdaptor::fetch_by_stable_id instead");

   return $self->get_GeneAdaptor->fetch_by_stable_id($geneid,$supporting);
}


=head2 get_Gene_by_Transcript_id

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_Gene_by_Transcript_id {
   my ($self,$tid, $supporting) = @_;

   $self->throw("Call to deprecated method " .
		"Bio::EnsEMBL::DBSQL::DBAdaptor::get_Gene_by_Transcript_id\n");
	
   return undef;
}


=head2 get_Gene_by_DBLink

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_Gene_by_DBLink {
   my ($self,$ext_id, $supporting) = @_;

   $self->throw("Call to deprecated method " .
		"Bio::EnsEMBL::DBSQL::DBAdaptor::get_Gene_by_DBLink\n");

   return undef;
}



=head2 get_Gene_array

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED 
              use Bio::EnsEMBL::DBSQL::GeneAdaptor::fetch_by_stable_id instead 
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_Gene_array {
    my ($self,@geneid) = @_;

    $self->throw("Call to deprecated method get_Gene_array");

    return undef;
}


=head2 get_Gene_array_supporting

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED 
               use Bio::EnsEMBL::GeneAdaptor::fetch_by_stable_id instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_Gene_array_supporting {
    my ($self,$supporting,@geneid) = @_;
    
    $self->throw("Call to deprecated method get_Gene_array_supporting");

    return undef;
}


=head2 get_Virtual_Contig_by_Transcript_id

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use 
               Bio::EnsEMBL::DBSQL::SliceAdaptor->fetch_by_transcript_stable_id
               instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_Virtual_Contig_by_Transcript_id {
   my ($self,$tid, $maxlen) = @_;

   $self->warn("get_Virtual_contig is deprecated. Use " .
	     "Bio::EnsEMBL::DBSQL::SliceAdaptor->fetch_by_transcript_stable_id"
	     . " instead.");

   return $self->get_SliceAdaptor->fetch_by_transcript_stable_id($tid,$maxlen);
}


=head2 get_Transcript_in_VC_coordinates

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED 
               use Bio::EnsEMBL::DBSQL::TranscriptAdaptor->fetch_by_stable_id
               and Transcript::transform(Bio::EnsEMBL::Slice) instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_Transcript_in_VC_coordinates{
   my ($self,$tid) = @_;

  $self->throw("call to deprecated method get_Transcript_in_VC_coordinates\n");

   return undef;
}


=head2 get_donor_locator

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_donor_locator {
    my ($self) = @_;

    $self->throw("call to deprecated method get_donor_locator\n");

    return undef;
}


=head2 get_last_update_offset

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_last_update_offset{
    my ($self) = @_;

    $self->throw("call to deprecated method get_donor_locator\n");

    return undef;
}    



=head2 get_last_update

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_last_update{
    my ($self) = @_;

    $self->throw("call to deprecated method get_last_update_offset\n");
    
    return undef;
}     


=head2 get_now_offset

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_now_offset{
    my ($self) = @_;

    $self->throw("call to deprecated method get_now_offset\n");

    return undef;
}



=head2 get_offset

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_offset{
    my ($self) = @_;

    $self->throw("call to deprecated method get_offset\n");

    return undef;
}
    


=head2 get_Protein_annseq

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_Protein_annseq{
    my ($self,$ENSP) = @_;

    $self->throw("call to deprecated method get_Protein_annseq\n");
    
    return undef;
} 


=head2 get_Transcript

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use 
               Bio::EnsEMBL::DBSQL::TranscriptAdaptor::fetch_by_dbID instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut
    
sub get_Transcript{
    my ($self,$transid) = @_;

    $self->warn("call to deprecated method get_Transcript " .
		"use Bio::EnsEMBL::DBSQL::TranscriptAdaptor::fetch_by_dbID " .
		"instead\n");
 
    return $self->get_TranscriptAdaptor->fetch_by_dbID($transid);
}



=head2 get_Translation

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use 
               Bio::EnsEMBL::DBSQL::TranslationAdaptor::fetch_by_dbID instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_Translation{
   my ($self,$translation_id) = @_;

   $self->warn("call to deprecated method get_Translation " .
	       "use Bio::EnsEMBL::DBSQL::TranslationAdaptor::fetch_by_dbID " .
	       "instead\n");
   
   return $self->get_TranslationAdaptor->fetch_by_dbID($translation_id);
}


=head2 get_Exon

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED 
               use Bio::EnsEMBL::DBSQL::ExonAdaptor::fetch_by_dbID instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_Exon{
   my ($self,$exonid) = @_;

   $self->warn("call to deprecated method get_Exon " .
	       "use Bio::EnsEMBL::DBSQL::ExonAdaptor::fetch_by_dbID instead");

   return $self->get_ExonAdaptor->fetch_by_dbID($exonid);
}


=head2 get_all_Gene_id

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_all_Gene_id{
   my ($self) = @_;

   $self->throw("call to deprecated method get_all_Gene_id\n");

   return undef;
}


=head2 get_all_Transcript_id

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_all_Transcript_id{
   my ($self) = @_;

   $self->throw("Call to deprecated method get_all_Transcript_id\n");

   return undef;
}



=head2 delete_Exon

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use Bio::EnsEMBL::DBSQL::ExonAdaptor::remove instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub delete_Exon{
    my ($self,$exon_id) = @_;

    $self->throw("call to deprecated method delete exon. use " .
		 "Bio::EnsEMBL::DBSQL::ExonAdaptor::remove instead");

    return undef;
}


=head2 delete_Supporting_Evidence

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub delete_Supporting_Evidence {
    my ($self,$exon_id) = @_;
 
    $self->throw("call to deprecated method delete_Supporting_Evidence\n");

    return undef;
}


=head2 delete_Features

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub delete_Features {
    my ($self,$contig) = @_;

    $self->throw("call to deprecated method delete_Features\n");

    return undef;
} 


=head2 delete_Gene

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub delete_Gene {
  my ($self, @args) = @_;

  $self->throw("call to deprecated method delete_Gene");

  return undef;
}


=head2 geneid_to_cloneid

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub geneid_to_cloneid {
  my ($self, @args) = @_;

  $self->throw("call to deprecated method geneid_to_cloneid");

  return undef;
}


=head2 write_Gene

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub write_Gene {
  my ($self, @args) = @_;

  $self->throw("call to deprecated method write_Gene");

  return undef;
}


=head2 write_all_Protein_features

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub write_all_Protein_features {
  my ($self, @args) = @_;

  $self->throw("call to deprecated method write_all_Protein_features");

  return undef;
}


=head2 write_Protein_feature

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub write_Protein_feature {
  my ($self, @args) = @_;

  $self->throw("call to deprecated method write_Protein_feature");

  return undef;
}


=head2 write_Feature

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub write_Feature {
  my ($self, @args) = @_;

  $self->throw("call to deprecated method write_Feature");

  return undef;
}


=head2 write_supporting_evidence

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub write_supporting_evidence {
  my ($self, @args) = @_;

  $self->throw("call to deprecated method write_supporting_evidence");

  return undef;
}


=head2 get_supporting_evidence

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_supporting_evidence {
  my ($self, @args) = @_;

  $self->throw("call to deprecated method get_supporting_evidence");

  return undef;
}


=head2 write_Analysis

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub write_Analysis {
  my ($self, @args) = @_;

  $self->throw("call to deprecated method write_Analysis");

  return undef;
}


=head2 exists_Homol_Feature

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub exists_Homol_Feature {
  my ($self, @args) = @_;

  $self->throw("call to deprecated method exists_Homol_Feature");

  return undef;
}


=head2 get_Analysis

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED 
               use Bio::EnsEMBL::DBSQL::AnalysisAdaptor::fetch_by_dbID
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_Analysis {
    my ($self,$id) = @_;

    $self->warn("call to deprecated method get_Analysis. Use " .
		"Bio::EnsEMBL::DBSQL::AnalysisAdaptor::fetch_by_DBID instead");
    
    $self->get_AnalysisAdaptor()->fetch_by_dbID($id);
} 


=head2 exists_Analysis

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub exists_Analysis {
  my ($self, @args) = @_;

  $self->throw("call to deprecated method exists_Analysis");

  return undef;
}


=head2 write_Transcript

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub write_Transcript {
  my ($self, @args) = @_;

  $self->throw("call to deprecated method write_Transcript");

  return undef;
}


=head2 write_Translation

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub write_Translation {
  my ($self, @args) = @_;

  $self->throw("call to deprecated method write_Translation");

  return undef;
}


=head2 write_Exon

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub write_Exon {
  my ($self, @args) = @_;

  $self->throw("call to deprecated method write_Exon");

  return undef;
}


=head2 get_PredictionFeature_as_Transcript

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_PredictionFeature_as_Transcript {
  my ($self, @args) = @_;

  $self->throw("call to deprecated method get_PredictionFeature_as_Transcript");

  return undef;
}


=head2 write_Clone

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub write_Clone {
  my ($self, @args) = @_;

  $self->throw("call to deprecated method write_Clone");

  return undef;

#    $clone_ad->store($clone);
   # my $clone_id = $clone->id;

#    $clone || $self->throw("Trying to write a clone without a clone object!\n");
#    if( !$clone->isa('Bio::EnsEMBL::DB::CloneI') ) {
#	$self->throw("Clone '$clone' is not a 'Bio::EnsEMBL::DB::CloneI'");
#    }
    
#    my @sql;
#    my $sql = "insert into clone(name, embl_acc, version, embl_version, htg_phase, created, modified) values('$clone_id', '".$clone->embl_id."', ".$clone->version.",".$clone->embl_version.", ".$clone->htg_phase.", FROM_UNIXTIME(".$clone->created."), FROM_UNIXTIME(".$clone->modified."))";
#    my $sth = $self->prepare($sql);
#    #my $sth = $self->prepare('insert into clone (clone_id, name,  embl_acc, version, embl_version, htg_phase, created, modified) values(?, ?, ?, ?, ?, ?, FROM_UNIXTIME(?), FROM_UNIXTIME(?)'); 
#    my $rv = $sth->execute();
        
#    $self->throw("Failed to insert clone $clone_id") unless $rv;
#    $sth = $self->prepare("select last_insert_id()");
#    my $res = $sth->execute;
#    my $row = $sth->fetchrow_hashref;
#    $sth->finish;
#    my $id  = $row->{'last_insert_id()'};
#    #print(STDERR "Clone $clone_id - $id\n");
    
#    foreach my $contig ( $clone->get_all_Contigs() ) {        
#        $self->write_Contig($contig,$id);
#    }
}


=head2 write_Contig

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED 
               use Bio::EnsEMBL::DBSQL::RawContigAdaptor::store instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub write_Contig {
    my($self, $contig, $clone)  = @_;
       

    $self->warn("call to deprecated method write_Contig. " .
	      "use Bio::EnsEMBL::DBSQL::RawContigAdaptor->store() instead.\n");
    
    my $rca = $self->get_RawContigAdaptor();

    $rca->store($contig, $clone);
    #Why do we have $clone if contig->cloneid is ok?
     
   # $self->throw("$contig is not a Bio::EnsEMBL::DB::ContigI - cannot insert contig for clone $clone")
#        unless $contig->isa('Bio::EnsEMBL::DB::ContigI');   
#    my $dna = $contig->primary_seq  || $self->throw("No sequence in contig object");
#    $dna->id                        || $self->throw("No contig id entered.");
#    $clone                          || $self->throw("No clone entered.");
    
##   (defined($contig->species)    && $contig->species   ->isa("Bio::EnsEMBL::Species"))    || $self->throw("No species object defined");
##    (defined($contig->chromosome) && $contig->chromosome->isa("Bio::EnsEMBL::Chromosome")) 
##                                    || $self->throw("No chromosomeobject defined");
                                    
##   my $species_id    = $self->write_Species   ($contig->species);
##   my $chromosome_id = $self->write_Chromosome($contig->chromosome,$species_id);    
#    my $contigid      = $contig->id;
#    my $len           = $dna   ->length;
#    my $seqstr        = $dna   ->seq;
#    my $offset        = $contig->embl_offset();
#    my $corder         = $contig->order();
#    #my $chromosome_id = $contig->chromosome->get_db_id;
#    my  $international_name = $contig->international_name();

    # Insert the sequence into the dna table
    #$self->_insertSequence($seqstr, $contig->seq_date);
    
 #   my @sql;
    
#    my $sth = $self->prepare("
#        insert into contig(name, dna_id, length, clone_id, offset, corder, international_name ) 
#        values(?, LAST_INSERT_ID(), ?, ?, ?, ?, ?)
#        "); 
#    #print STDERR "contig name = ",$contigid,"\n";
#    my $rv = $sth->execute(
#        $contigid,			   
#        $len,
#        $clone,
#        $offset,
#        $corder,
#        $international_name
#        );  
          
#    $self->throw("Failed to insert contig $contigid") unless $rv;
       
    
#    $sth = $self->prepare("select last_insert_id()");
#    $sth->execute;
#    my ($id) = $sth->fetchrow
#        or $self->throw("Failed to get last insert id");
#    #can no longer do this as get_all_SeqFeatures no longer exists
#    #if a contig is written to the database
#    # this is a nasty hack. We should have a cleaner way to do this.
#    #my @features = $contig->get_all_SeqFeatures;
#    #print(STDERR "Contig $contigid - $id\n"); 
#    # write sequence features. We write all of them together as it
#    # is more efficient
#    #$self->get_Feature_Obj->write($contig, @features);
    
    return 1;
}



=head2 _insertSequence

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub _insertSequence {
  my ($self, @args) = @_;

  $self->throw("call to deprecated method _insertSequence");

  return undef;

    #if ($self->dnadb ne $self) {
#      $self->throw("ERROR: Trying to write to a remote dna database");
#    } 
    
#    my $statement = $self->prepare("
#        insert into dna(sequence,created) 
#        values(?, FROM_UNIXTIME(?))
#        "); 
        
#    my $rv = $statement->execute($sequence, $date); 
    
#    $self->throw("Failed to insert dna $sequence") unless $rv;    
}


=head2 write_Species

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub write_Species {
  my ($self, @args) = @_;

  $self->throw("call to deprecated method write_Species");

  return undef;

    #if (!defined($species)) {
#	$self->throw("No species argument input");
#    }
#    if (!$species->isa("Bio::EnsEMBL::Species")) {
#	$self->throw("[$species] is not a Bio::EnsEMBL::Species object");
#    }

#    my $query = "select species_id " .
#	        "from   species " .
#		"where  nickname    = '" . $species->nickname . "' " . 
#		"and    taxonomy_id = "  . $species->taxonomy_id;

#    my $sth = $self->prepare($query);
#    my $res = $sth->execute;

#    if ($sth->rows == 1) {
#	my $rowhash    = $sth->fetchrow_hashref;
#	my $species_id = $rowhash->{species_id};
#	return $species_id;
#    } 

#    $query =  "insert into species(species_id,nickname,taxonomy_id) " . 
#	      "            values(null,'" . $species->nickname . "'," . $species->taxonomy_id . ")";
	
    
#    $sth = $self->prepare($query);
#    $res = $sth->execute;

#    $sth = $self->prepare("select last_insert_id()");
#    $res = $sth->execute;

#    my $rowhash = $sth->fetchrow_hashref;
#    my $species_id = $rowhash->{'last_insert_id()'};
   
#    return $species_id;
}


=head2 get_Update_Obj

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_Update_Obj {
  my ($self, @args) = @_;

  $self->throw("call to deprecated method get_Update_Obj");

  return undef;

 #   my( $update_obj );
#    unless ($update_obj = $self->{'_update_obj'}) {
#        require Bio::EnsEMBL::DBSQL::Update_Obj;
#        $update_obj = Bio::EnsEMBL::DBSQL::Update_Obj->new($self);
#        $self->{'_update_obj'} = $update_obj;
#    }
#    return $update_obj;
}



=head2 get_all_chr_ids

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_all_chr_ids {
   my ($self, $type) = @_;

   $self->warn("call to deprecated method get_all_chr_ids " . 
	    "Use Bio::EnsEMBL::DBSQL::ChromosomeAdaptor::fetch_all instead\n");

   return $self->get_ChromosomeAdaptor()->fetch_all();

#   $self->throw("no static_gold_path given") unless defined $type;
#   my @out;

#   my $q= "SELECT DISTINCT chromosome_id 
#           FROM assembly
#           WHERE type = '$type'";
#   my $sth = $self->prepare($q) || $self->throw("can't prepare: $q");
#   my $res = $sth->execute || $self->throw("can't prepare: $q");

#   while( my ($id) = $sth->fetchrow_array) {
#       push(@out, $id);
#   }
#   return @out;
}


=head2 get_all_fpcctg_ids

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_all_fpcctg_ids {
   my ($self, $type) = @_;

   $self->warn("DBAdaptor->get_all_fpcctg_ids is deprecated. \n" .
	       'Use $dba->get_StaticGoldenContigAdaptor->get_all_fpc_ids($id)'.
	       "instead");

   return $self->get_StaticGoldenPathAdaptor()->get_all_fpc_ids($type);

#  $self->throw("no static_gold_path given") unless defined $type;
#   my @out;

#   my $q= "SELECT DISTINCT superctg_name 
#           FROM assembly 
#           WHERE type = '$type'";
#   my $sth = $self->prepare($q) || $self->throw("can't prepare: $q");
#   my $res = $sth->execute || $self->throw("can't prepare: $q");

#   while( my ($id) = $sth->fetchrow_array) {
#       push(@out, $id);
#   }
#   return @out;
}



=head2 get_object_by_wildcard

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_object_by_wildcard {
  my ($self, @args) = @_;

  $self->throw("call to deprecated method get_object_by_wildcard");

  return undef;

 #  print STDERR "Got type: $type and string: $string\n";
#   my @ids;
#   my $sth = $self->prepare("select id from $type where id like \'$string\'");
#   print STDERR "mysql: select id from $type where id like \'$string\'\n";
#   my $res = $sth->execute || $self->throw("Could not get any ids!");
#   while( my $rowhash = $sth->fetchrow_hashref) {
#       push(@ids,$rowhash->{'id'});
#   }
   
#   if ($type eq 'gene') {
#       return $self->gene_Obj->get_array_supporting('without',@ids);
#   }
#   if ($type eq 'transcript') {
#       my @trans;
#       foreach my $id (@ids) {
#	   push @trans, $self->gene_Obj->get_Transcript($id);
#       }
#       return @trans;
#   }
#   if ($type eq 'exon') {
#       my @exons;
#       foreach my $id (@ids) {
#	   push @exons, $self->gene_Obj->get_Exon($id);
#       }
#       return @exons;
#   }
#   if ($type eq 'clone') {
#       my @clones;
#       foreach my $id (@ids) {
#	   push @clones, $self->get_Clone($id);
#       }
#       return @clones;
#   }
#   else {
#       $self->throw("Type $type not supported, only gene, transcript, exon and clone\n");
#   }
#   return;
}



=head2 write_Chromosome

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub write_Chromosome {
  my ($self, @args) = @_;

  $self->throw("call to deprecated method write_Chromosome");

  return undef;

#  $self->throw("No chromosome argument input") unless defined($chromosome);
   
  
#  if (!$chromosome->isa("Bio::EnsEMBL::Chromosome")) {
#    $self->throw("[$chromosome] is not a Bio::EnsEMBL::Chromosome object");
#  }
#  if(!$length){
#    $length = 0;
#  }
#  if(!$known_genes){
#    $known_genes = 0;
#  }
#  if(!$unknown_genes){
#    $unknown_genes = 0;
#  }
#  if(!$snps){
#    $snps = 0;
#  }
  
#  my $query = "select chromosome_id " .
#              "from   chromosome " .
#              "where  name       = '" . $chromosome->chr_name . "' " .
#	      " and    known_genes = "  . $known_genes . 
#	      " and    unknown_genes = ".$unknown_genes .
#	      " and    snps = ".$snps.
#	      " and    length = ".$length;
  
#    my $sth = $self->prepare($query);
#    my $res = $sth->execute;

#    if ($sth->rows == 1) {
#	my $rowhash       = $sth->fetchrow_hashref;
#	my $chromosome_id = $rowhash->{chromosome_id};
#	return $chromosome_id;
#    } 

#    $query =  "insert into chromosome(chromosome_id,name,known_genes,unknown_genes,snps,length) values(null,'" . $chromosome->chr_name . "',".$known_genes.",".$unknown_genes.",".$snps.",".$length.")";
	
#  print $query."\n";
#    $sth = $self->prepare($query);
#    $res = $sth->execute;

#    $sth = $self->prepare("select last_insert_id()");
#    $res = $sth->execute;

#    my $rowhash       = $sth->fetchrow_hashref;
#    my $chromosome_id = $rowhash->{'last_insert_id()'};
   
#    return $chromosome_id;
  }



=head2 _analysis_cache

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub _analysis_cache {
  my ($self, @args) = @_;

  $self->throw("call to deprecated method _analysis_cache");

  return undef;
#   if( @_ ) {
#      my $value = shift;
#      $obj->{'_analysis_cache'} = $value;
#    }
#    return $obj->{'_analysis_cache'};
}


=head2 _contig_seq_cache

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub _contig_seq_cache {
  my ($self, @args) = @_;

  $self->throw("call to deprecated method _contig_seq_cache");

  return undef;

#   if( $seq ) {
       
#       #
#       # Every 100 hits, flush the cache
#       #
#       if( $self->{'_contig_seq_cnt'} > 100 ) {
#	   $self->_flush_seq_cache;
#	   $self->{'_contig_seq_cnt'} = 0;
#       }

#       $self->{'_contig_seq_cnt'}++;
#       $self->{'_contig_seq_cache'}->{$id} = $seq;
#   }

#   return $self->{'_contig_seq_cache'}->{$id};
}


=head2 _flush_seq_cache

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub _flush_seq_cache {
  my ($self, @args) = @_;

  $self->throw("call to deprecated method _flush_seq_cache");

  return undef;
#   $self->{'_contig_seq_cache'} = {};
}


=head2 _lock_tables

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub _lock_tables {
 my( $self, @tables) = @_;
 $self->warn("call to deprecated method _lock_tables");   
 
 my $state;
 foreach my $table ( @tables ) {
   if( $self->{'_lock_table_hash'}->{$table} == 1 ) {
     $self->warn("$table already locked. Relock request ignored");
   } else {
     if( $state ) { $state .= ","; } 
     $state .= "$table write";
     $self->{'_lock_table_hash'}->{$table} = 1;
   }
 }
 
 my $sth = $self->prepare("lock tables $state");
 my $rv = $sth->execute();
 $self->throw("Failed to lock tables $state") unless $rv;

}


=head2 _unlock_tables

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub _unlock_tables {
   my ($self,@tables) = @_;
   
   $self->throw("call to deprecated method _unlock_tables");

   my $sth = $self->prepare("unlock tables");
   my $rv  = $sth->execute();
   $self->throw("Failed to unlock tables") unless $rv;
   %{$self->{'_lock_table_hash'}} = ();
}


=head2 get_Clone

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED Use 
               Bio::EnsEMBL::DBSQL::CloneAdaptor::fetch_by_accession_version 
               or Bio::EnsEMBL::DBSQL::CloneAdaptor::fetch_by_accession instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_Clone { 
    my( $self, $accession ) = @_;

    $self->warn('DBAdaptor->get_Clone() is deprecated.\n' .
	  'Use \$db->get_CloneAdaptor()->fetch_by_accession_version(\$acc_ver)'
	. ' or $db->get_CloneAdaptor()->fetch_by_accession(\$acc) instead');

    my $ca = $self->get_CloneAdaptor;

    if ($accession =~ /(.+?)\.(\d+)/) {
	$accession = $1;
	my $version   = $2;
	return $ca->fetch_by_accession_version($accession, $version);
    }
    else {
	return $ca->fetch_by_accession($accession);
    }

}


=head2 list_embl_version_by_Clone

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use
               Bio::EnsEMBL::DBSQL::CloneAdaptor::list_embl_versions_by_Clone 
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub list_embl_version_by_Clone { 
    my( $self, $accession ) = @_;

    $self->warn('DBAdaptor->list_embl_version_by_Clone() is deprecated.\n' .
	'Use \$db->get_CloneAdaptor()->list_embl_version_by_accession(\$acc)');

    my $ca = $self->get_CloneAdaptor;

    return $ca->list_embl_version_by_accession($accession);
}


=head2 get_Clone_by_Version

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use 
               Bio::EnsEMBL::DBSQL::CloneAdaptor::fetch_by_accession_version
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_Clone_by_version { 
    my ($self,$accession,$ver) = @_;

    $self->warn('DBAdaptor->get_CLone_by_version() is deprecated.\n' .
      'Use \$db->get_CloneAdaptor()->fetch_by_accession_version(\$acc,\$ver)');

    my $ca = $self->get_CloneAdaptor;

    return $ca->fetch_by_accession_version($accession,$ver);
}
  


=head2 get_Contig

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED 
               use Bio::EnsEMBL::DBSQL::RawContigAdaptor::fetch_by_name instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_Contig {
    my ($self,$id) = @_;

    $self->warn("DBAdaptor->get_Contig() is a deprecated.\n" .
	       "Use \$db->get_RawContigAdaptor()->fetch_by_name(\$id)"); 

    return $self->get_RawContigAdaptor->fetch_by_name($id);
}


=head2 get_Contig_by_internal_id

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED 
               use Bio::EnsEMBL::DBSQL::RawContigAdaptor::fetch_by_dbID instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_Contig_by_internal_id {
  my ($self,$id) = @_;

  $self->warn("DBAdaptor get_Contig_by_internal_id is deprecated.\n" .
	      "Use \$db->get_RawContigAdaptor()->fetch_by_dbID(\$id) "); 

  return $self->get_RawContigAdaptor()->fetch_by_dbID($id);

#  if (!defined($id)) {
#    $self->throw("No id defined\n");
#  }
#  my $query = "select id from contig where internal_id = $id";

#  my $sth = $self->prepare($query);
#  my $res = $sth->execute;

#  if ($sth->rows < 1) {
#    $self->throw("No contig available with id $id\n");
#  }
#  my $ref = $sth->fetchrow_hashref;
#  my $contigid = $ref->{'id'};

#  return $self->get_Contig($contigid);
}


=head2 get_contig_by_international_id

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_contig_by_international_id {
  my ($self, @args) = @_;

  $self->throw("call to deprecated method get_contig_by_international_id");

  return undef;
#   my $sth=$self->prepare("select id from contig where international_id = '$int_id'");
#   $sth->execute;
#   my $row = $sth->fetchrow_hashref;
#   my $id  = $row->{'id'};

#   return $self->get_Contig($id);
}


=head2 get_Contigs_by_Chromosome

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_Contigs_by_Chromosome {
  my ($self, @args) = @_;

  $self->throw("call to deprecated method get_Contigs_by_Chromosome");

  return undef;
}



=head2 get_all_Clone_id

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_all_Clone_id {
  my ($self, @args) = @_;

  $self->throw("call to deprecated method get_all_Clone_id");

  return undef;
}



=head2 get_all_Contig_id

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_all_Contig_id {
  my ($self, @args) = @_;

  $self->throw("call to deprecated method get_all_Contig_id");

  return undef;
}


=head2 perl_only_sequences

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub perl_only_sequences{
   my $self = shift;

   $self->warn('DBAdaptor->perl_only_sequences() is deprecated.\n' .
	       'no replacement has been written');

#   if( @_ ) {
#      my $value = shift;
#      $self->{'perl_only_sequences'} = $value;
#    }
#    return $obj->{'perl_only_sequences'};

}


=head2 perl_only_contigs

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub perl_only_contigs{
   my $self = shift;

   $self->warn('DBAdaptor->perl_only_contigs() is deprecated.\n' .
	       'no replacement has been written');

#   if( @_ ) {
#      my $value = shift;
#      $obj->{'perl_only_contigs'} = $value;
#    }
#    return $obj->{'perl_only_contigs'};

}


=head2 _crossdb

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub _crossdb {
   my $self = shift;

   $self->throw('DBAdaptor->_crossdb is deprecated.\n' .
		'No replacement has been written');

#   if( @_ ) {
#      my $value = shift;
#      $obj->{'_crossdb'} = $value;
#    }
#    return $obj->{'_crossdb'};
}


=head2 create_tables

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub create_tables { 
  my $self = shift;

  $self->throw('DBAdaptor->create_tables is deprecated.\n' .
	       'No replacement has been written');

  # get all adaptors ...
  # call create_tables on them

  # create some tables without adaptors
  # (which should disappear once)
}


=head2 get_FamilyAdaptor

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use 
               Bio::EnsEMBL::ExternalData::Family::DBSQL::DBAdaptor instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

## following is not part of core EnsEMBL, so maybe doesn't belong here and
## has to be moved elsehwere (e.g. as part of a more dynamical
## add_external_adaptor scheme). For convenience I have it here, now,
## though. It will break things for people who don't have ensembl-external
## checked out ...

sub get_FamilyAdaptor {
    my( $self ) = @_;
    
    $self->throw("DBAdaptor->get_FamilyAdaptor is deprecated.\n" .
"Family db has now its own DBAdaptor. Use it to connect a family database,\n" .
"Keep it cached using add_ExternalAdaptor from the core DBAdaptor\n" .
"From the family DBAdaptor you can then get_FamilyAdaptor\n" .
"For more info, see perldoc Bio::EnsEMBL::ExternalData::Family::DBSQL::DBAdaptor\n");
    return undef;

#    my( $fa );
#    unless ($fa = $self->{'_externaldata_family_familyadaptor'}) {
#        eval{
#            require Bio::EnsEMBL::ExternalData::Family::FamilyAdaptor;
#        };
#        if ($@) {
#            $self->throw(
#                "Unable to load 'Bio::EnsEMBL::ExternalData::Family::FamilyAdaptor'\n"
#                . "It is not part of the core Ensembl distribution.\n"
#                . "Have you installed it?");
#        }
#        $fa = Bio::EnsEMBL::ExternalData::Family::FamilyAdaptor->new($self);
#        $self->{'_externaldata_family_familyadaptor'} = $fa;
#    }
#    return $fa;
}


=head2 find_GenomeHits

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub find_GenomeHits {
  my ($self, @args) = @_;

  $self->throw("call to deprecated method find_GenomeHits");

  return undef;
}


=head2 release_number

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub release_number {
  my ($self, @args) = @_;

  $self->throw("call to deprecated method release_number");

  return undef;
}



1;
