
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

Post questions to the EnsEMBL development list <ensembl-dev@ebi.ac.uk>

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

  return $self;
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

  my $lite = $self->get_db_adaptor('lite');
  my $primary_adaptor;

  if($lite) {
    $primary_adaptor = $lite->get_SNPAdaptor();
  } else {
    my $snp = $self->get_db_adaptor('SNP');
    
    unless($snp) {
      warn("No lite or SNP database, cannot get snp adaptor\n");
      return undef;
    }

    $primary_adaptor = $snp->get_SNPAdaptor();
    $primary_adaptor->ensembl_db( $self );
  }
  
  #return a proxy adaptor which can use the lite or the core database
  return $self->_get_adaptor("Bio::EnsEMBL::DBSQL::ProxySNPAdaptor",
			     $primary_adaptor);
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

  my $lite = $self->get_db_adaptor('lite');

  if( defined $lite ) {
    return $lite->get_LandmarkMarkerAdaptor();
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

   #return the sequence adaptor for the dnadb (which may be this db)
   return $self->dnadb->_get_adaptor("Bio::EnsEMBL::DBSQL::SequenceAdaptor");
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
    $self->dnadb->_get_adaptor("Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptor");
  
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

  return $self->dnadb->_get_adaptor("Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor");
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

  return $self->_get_adaptor("Bio::EnsEMBL::DBSQL::SupportingFeatureAdaptor");
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

  return $self->_get_adaptor('Bio::EnsEMBL::Map::DBSQL::MarkerFeatureAdaptor');
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

  return $self->_get_adaptor('Bio::EnsEMBL::Map::DBSQL::MarkerAdaptor');
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
      my $ass;
      eval {
        $ass = $obj->dnadb->get_MetaContainer()->get_default_assembly();
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
      $self->{'_xf_adaptors'}->{$key} = undef;
    }
  }


  $self->{'dnadb'} = undef;

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
    $self->throw("[$adaptor] is not a " .
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


1;
