
=Head1 NAME - Bio::EnsEMBL::DBSQL::DBAdaptor

=head1 SYNOPSIS

    $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -user   => 'root',
        -dbname => 'pog',
        -host   => 'caldy',
        -driver => 'mysql',
        );

    $gene_adaptor = $db->get_GeneAdaptor();

    $gene = $gene_adaptor()->fetch_by_stable_id($stable_id);

    $slice = $db->get_SliceAdaptor()->fetch_by_dbID($slice_id);

    
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
  Arg  1    : scalar string $source
              The source of info in the database connected to (e.g. 'embl')
  Function  : Sets/Gets the source or human readable name of the genes in 
              the connected database. For example for the sanger db the source
              would be 'sanger'.
  Returntype: scalar string
  Exceptions: none
  Caller    : Bio::EnsEMBL::GeneAdaptor  Bio::EnsEMBL::LiteGeneAdaptor EnsWeb

=cut

sub source {
  my ($self, $source) = @_;

  if(defined $source) {
    $self->{'_source'} = $source;
  }

  return $self->{'_source'};
}


=head2 get_MetaContainer

 Title   : get_Meta_Container
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


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







=head2 add_DASFeatureFactory

 Title   : add_DASFeatureFactory
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_DASFeatureFactory{
   my ($self,$value) = @_;

   unless( ref $value && $value->isa('Bio::EnsEMBL::DB::ExternalFeatureFactoryI') ) {
       $self->throw("[$value] is not a Bio::EnsEMBL::DB::ExternalFeatureFactoryI but it should be!");
   }

   push(@{$self->{'_das_ff'}},$value);
}

=head2 _each_DASFeatureFactory

 Title   : _each_DASFeatureFactory
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _each_DASFeatureFactory{
   my ($self) = @_;

   return @{$self->{'_das_ff'}}
}





=head2 get_ProteinFeatureAdaptor

 Title   : get_ProteinFeatureAdaptor
 Usage   : my $pfa = $dba->get_ProteinFeatureAdaptor(); 
 Function: Gets a ProteinFeatureAdaptor for this database.  
           Formerly named get_Protfeat_Adaptor()
 Example : my $protfeat = $dba->get_ProteinFeatureAdaptor->fetch_by_dbID($id);
 Returns : Bio::EnsEMBL::DBSQL::ProteinFeatureAdaptor
 Args    : none

=cut

sub get_ProteinFeatureAdaptor {
    my( $self ) = @_;
    
    return 
      $self->_get_adaptor("Bio::EnsEMBL::DBSQL::ProteinFeatureAdaptor");
}


=head2 get_ProteinAdaptor

 Title   : get_ProteinAdaptor
 Usage   : my $pa = $dba->get_ProteinAdaptor(); 
 Function: Gets a ProteinAdaptor for this database.  
           Formerly named get_Protein_Adaptor()
 Example : my $prot = $dba->get_ProteinAdaptor->fetch_by_dbID($id);
 Returns : Bio::EnsEMBL::DBSQL::ProteinAdaptor
 Args    : none

=cut

sub get_ProteinAdaptor {
    my( $self ) = @_;
 
    return $self->_get_adaptor("Bio::EnsEMBL::DBSQL::ProteinAdaptor");
}


sub get_SNPAdaptor {
  my ($self)  = @_;

  #return a proxy adaptor which can use the lite or the core database
  return $self->_get_adaptor("Bio::EnsEMBL::DBSQL::ProxySNPAdaptor");
}


sub get_MapFragAdaptor {
  my $self = shift;

  return $self->_get_adaptor( "Bio::EnsEMBL::DBSQL::MapFragAdaptor" );
}



=head2 get_CloneAdaptor

    my $ca = $dba->get_CloneAdaptor;

Returns a B<Bio::EnsEMBL::DBSQL::CloneAdaptor>
object, which is used for reading and writing
B<Clone> objects from and to the SQL database.

=cut 

sub get_CloneAdaptor {
  my( $self ) = @_;
  
  return $self->_get_adaptor("Bio::EnsEMBL::DBSQL::CloneAdaptor");
}


sub get_LandmarkMarkerAdaptor {
  my $self = shift;

  if( defined $self->lite_DBAdaptor() ) {
    return $self->lite_DBAdaptor()->get_LandmarkMarkerAdaptor();
  } else {
    return undef;
  }
}


=head2 get_PredictionTranscriptAdaptor

  Args      : none
  Function  : PredictionTranscript Adaptor for this database.
  Returntype: Bio::EnsEMBL::DBSQL::PredictionTranscriptAdaptor
  Exceptions: none
  Caller    : general

=cut

sub get_PredictionTranscriptAdaptor {
   my ($self) = @_;

   return $self->_get_adaptor("Bio::EnsEMBL::DBSQL::PredictionTranscriptAdaptor");
 }

=head2 get_SequenceAdaptor

  Args      : none
  Function  : The sequence producing adaptor. Could be hooked to a different
              database than the rest for example.
  Returntype: Bio::EnsEMBL::DBSQL::SequenceAdaptor
  Exceptions: none
  Caller    : general, Bio::EnsEMBL::Slice, Bio::EnsEMBL::RawContig

=cut

sub get_SequenceAdaptor {
   my $self = shift;

   return $self->_get_adaptor("Bio::EnsEMBL::DBSQL::SequenceAdaptor");
}

=head2 lite_DBAdaptor

  Arg    [1]: Bio::EnsEMBL::Lite::DBAdaptor $liteDBConnection
              From there you get the denormalized GeneAdaptor.
  Function  : The liteDB is set in EnsWEB to get denormalized access to
              EnsEMBL data. It provides a GeneAdaptor that makes half
              filled Genes.
  Returntype: Bio::EnsEMBL::Lite::DBAdaptor
  Exceptions: 
  Caller    : set in EnsWEB, get internal

=cut


sub lite_DBAdaptor {
  my ($self, $arg ) = @_;
  if ( defined $arg ) {
    $self->{_liteDB} = $arg;

  }

  return $self->{_liteDB};
}


sub SNP_DBAdaptor {
  my ($self, $arg) = @_;

  if(defined $arg) {
    $self->{_SNP_db} = $arg;
  }

  return $self->{_SNP_db};
}



sub map_DBAdaptor {
  my ($self, $arg ) = @_;
  if ( defined $arg ) {
    $self->{_mapDB} = $arg;
  }

  return $self->{_mapDB};
}


sub est_DBAdaptor {
  my ($self, $arg) = @_;
  
  if(defined $arg) {
    $self->{_estDB} = $arg;
  }

  return $self->{_estDB};
}

sub get_GeneAdaptor {
    my( $self ) = @_;
    #get a core db adaptor
    my $core_adaptor = $self->_get_adaptor("Bio::EnsEMBL::DBSQL::GeneAdaptor");

    #use a proxy gene adaptor, capable of making decisions with regards to the
    #database that it uses, passing in the core adaptor as a constructor arg
    return $self->_get_adaptor("Bio::EnsEMBL::DBSQL::ProxyGeneAdaptor",
			     $core_adaptor);
  }



sub get_ExonAdaptor {
  my( $self ) = @_;
  
  return $self->_get_adaptor("Bio::EnsEMBL::DBSQL::ExonAdaptor");
}

sub get_TranscriptAdaptor {
  my( $self ) = @_;
  
  return $self->_get_adaptor("Bio::EnsEMBL::DBSQL::TranscriptAdaptor");
}

sub get_TranslationAdaptor {
    my( $self ) = @_;

    return $self->_get_adaptor("Bio::EnsEMBL::DBSQL::TranslationAdaptor");
}

sub get_FeatureAdaptor {
    my( $self ) = @_;

    return $self->_get_adaptor("Bio::EnsEMBL::DBSQL::FeatureAdaptor");
}

sub get_RawContigAdaptor {
    my( $self ) = @_;

    return $self->_get_adaptor("Bio::EnsEMBL::DBSQL::RawContigAdaptor");
}

sub get_SliceAdaptor {
  my( $self ) = @_;
  
  return $self->_get_adaptor("Bio::EnsEMBL::DBSQL::SliceAdaptor");
}

=head2 get_AnalysisAdaptor

 Title   : get_AnalysisAdaptor
 Usage   : $analysisAdaptor = $dbObj->get_AnalysisAdaptor;
 Function: gives the adaptor to fetch/store Analysis objects.
 Example :
 Returns : the adaptor
 Args    :

=cut

sub get_AnalysisAdaptor {
    my( $self ) = @_;

#    print "Getting an analysis adaptor from" . $self->dbname() . "\n";

    return $self->_get_adaptor("Bio::EnsEMBL::DBSQL::AnalysisAdaptor");
}

=head2 get_SimpleFeatureAdaptor

 Title   : get_SimpleFeatureAdaptor
 Usage   : $sa = $db->get_SimpleFeatureAdaptor;
 Example :
 Returns : the adaptor
 Args    : none

=cut

sub get_SimpleFeatureAdaptor {
  my( $self ) = @_;
  
  return $self->_get_adaptor("Bio::EnsEMBL::DBSQL::SimpleFeatureAdaptor");
}


=head2 get_RepeatConsensusAdaptor

 Title   : get_SimpleFeatureAdaptor
 Usage   : $sa = $db->get_SimpleFeatureAdaptor;
 Example :
 Returns : the adaptor
 Args    : none

=cut

sub get_RepeatConsensusAdaptor {
  my( $self ) = @_;
  
  return $self->_get_adaptor("Bio::EnsEMBL::DBSQL::RepeatConsensusAdaptor");
}

sub get_RepeatFeatureAdaptor {
  my( $self ) = @_;
  
  my $core_adaptor = 
    $self->_get_adaptor("Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptor");
  
  #create a proxy adaptor, using a core RepeatFeatureAdaptor as constructor arg
  
  return $self->_get_adaptor("Bio::EnsEMBL::DBSQL::ProxyRepeatFeatureAdaptor",
			    $core_adaptor);
}

=head2 get_ProteinAlignFeatureAdaptor

 Title   : get_ProteinAlignFeatureAdaptor
 Usage   : $sa = $db->get_ProteinAlignFeatureAdaptor;
 Returns : the adaptor


=cut

sub get_ProteinAlignFeatureAdaptor {
  my( $self ) = @_;
    
  return $self->_get_adaptor("Bio::EnsEMBL::DBSQL::ProteinAlignFeatureAdaptor");
}


=head2 get_DnaAlignFeatureAdaptor

 Args      : none
 Function  : Returns the DnaAlignFeatuire Adaptor which is connected to
             this database. There should only be one around.
 Returntype: Bio::EnsEMBL::DBSQL::DnaAlignFeatureAdaptor
 Exceptions: none
 Caller    : FeatureAdaptor, general

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

 Title   : get_AssemblyMapperAdaptor
 Usage   : $sa = $db->get_AssemblyMapperAdaptor;
 Example :
 Returns : the adaptor
 Args    :

=cut

sub get_AssemblyMapperAdaptor {
  my( $self ) = @_;
    
  return $self->_get_adaptor("Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor");
}


sub get_DBEntryAdaptor {
    my( $self ) = @_;

    return $self->_get_adaptor("Bio::EnsEMBL::DBSQL::DBEntryAdaptor");
}


=head2 get_StaticGoldenPathAdaptor

 Title   : get_StaticGoldenPathAdaptor
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_StaticGoldenPathAdaptor{
  my( $self ) = @_;
  
  return $self->_get_adaptor("Bio::EnsEMBL::DBSQL::StaticGoldenPathAdaptor");
}


=head2 get_KaryotypeBandAdaptor

 Title   : get_KaryotypeBandAdaptor
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_KaryotypeBandAdaptor {
    my( $self ) = @_;

    return $self->_get_adaptor("Bio::EnsEMBL::DBSQL::KaryotypeBandAdaptor");
}

=head2 get_ChromosomeAdaptor

 Title   : get_ChromosomeAdaptor
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_ChromosomeAdaptor {
    my( $self ) = @_;

    return $self->_get_adaptor("Bio::EnsEMBL::DBSQL::ChromosomeAdaptor");
}

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

    Title   : deleteObj
    Usage   : $dbObj->deleteObj
    Function: Call when you are done with this object. Breaks links between objects. Necessary to clean up memory.
    Example : -
    Returns : -
    Args    : -


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

 Title   : assembly_type
 Usage   : $obj->assembly_type($newval)
 Function: 
 Example : 
 Returns : name of assembly
 Args    : newvalue (optional)


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




################################################################## 
# 
# SUPPORT FOR EXTERNAL ADAPTORS & FEATURE FACTORIES 
# 
# These are not implemented on the new main trunk and at this
# point it is not clear if they ever will be.
#
##################################################################



=head2 extension_tables

 Title   : extension_tables
 Usage   : $obj->extension_tables($newval)
 Function: 
 Returns : value of extension_tables
 Args    : newvalue (optional)


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

 Title   : list_ExternalAdaptors
 Usage   : $obj->list_ExternalAdaptors
 Function: returns all the names of installed external adaptors
 Returns : a (possibly empty) list of name of external adaptors
 Args    : none

=cut

sub list_ExternalAdaptors {
    my ($self) = @_;

    $self->warn("DBAdaptor::list_ExternalAdaptors is not implmented. It will "
		. "either be implemented or deprecated at a later date\n");

    return undef;

    #return keys % {$self->{_ext_adaptors}};
}

=head2 add_ExternalAdaptor

 Title   : add_ExternalAdaptor
 Usage   : $obj->add_ExternalAdaptor('family', $famAdaptorObj);
 Function: adds the external adaptor the internal hash of known 
           external adaptors. If an adaptor of the same name is installed, 
           it will be overwritten.
 Returns : undef
 Args    : a name and a adaptor object. 

=cut

sub add_ExternalAdaptor {
    my ($self, $adtor_name, $adtor_obj) = @_;

    $self->warn("DBAdaptor::add_ExternalAdaptor is not implemented. It will "
		. "either be implemented or deprecated at a later date\n");

    #$self->_ext_adaptor($adtor_name, $adtor_obj);
    #undef;
}

=head2 get_ExternalAdaptor

 Title   : get_ExternalAdaptor
 Usage   : $obj->get_ExternalAdaptor('family');
 Function: retrieve external adaptor by name
 Returns : an adaptor (sub-type of BaseAdaptor) or undef
 Args    : the name 

=cut

sub get_ExternalAdaptor {
    my ($self, $adtor_name) = @_;
    
    $self->warn("DBAdaptor::get_ExternalAdaptor is not implemented. It will "
		. "either be implemented or deprecated at a later date\n");

    return undef;

    #$self->_ext_adaptor($adtor_name);
}


=head2 remove_ExternalAdaptor

 Title   : remove_ExternalAdaptor
 Usage   : $obj->remove_ExternalAdaptor('family')
 Function: removes the named external adaptor from the internal hash of known 
           external adaptors. If the adaptor name is not known, nothing 
           happens. 
 Returns : undef
 Args    : a name

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

 Title   : add_ExternalFeatureFactory
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


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

 Title   : _each_ExternalFeatureFactory
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


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

 Title   : _ext_adaptor
 Usage   : $obj->_ext_adaptor('family' [, $famAdaptorObj] )
 Function: 
 Returns : an adaptor or undef
 Args    : a name and a adaptor object. 

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
    
 Title   : DEPRECATED get_LiteAdaptor
 Usage   : DEPRECATED my $la = $db->get_LiteAdaptor;
 Function: DEPRECATED Returns the lite database object handler
 Example : DEPRECATED
 Returns : DEPRECATED Bio::EnsEMBL::DBSQL::LiteAdaptor
 Args    : DEPRECATED

=cut

sub get_LiteAdaptor {
    my( $self ) = @_;
    
    $self->throw("The lite adaptor is deprecated. Use the " .
		"Bio::EnsEMBL::Lite::DBAdaptor instead.\n" .
		 "This may be attached to the core DBAdaptor using the" .
		."lite_DBAdaptor method\n");

    return undef;
}


=head2 feature_Obj
    
 Title   : feature_Obj
 Usage   : my $featureobj = $db->feature_Obj
 Function: Returns the feature object database handle
 Example : 
 Returns : Bio::EnsEMBL::DB::Feature_ObjI
 Args    : 

=cut

sub feature_Obj {
    my $self = shift;
    $self->throw("No more Feature Objs!");

}




sub add_db_adaptor {
  my ($self, $name, $adaptor) = @_;

  $self->{'_db_adaptors'}->{$name} = $adaptor;
}

sub remove_db_adaptor {
  my ($self, $name) = @_;

  my $adaptor = $self->{'_db_adaptors'}->{$name};
  delete $self->{'_db_adaptors'}->{$name};

  return $adaptor;
}

sub get_all_db_adaptors {
  my ($self, $name) = @_;   

  return $self->{'_db_adaptors'};
}

sub get_db_adaptor {
  my ($self, $name) = @_;

  return $self->{'_db_adaptors'}->{$name};
}


=head2 get_Gene

 Title   : get_Gene
 Usage   : $obj->get_Gene($geneid, $supporting)
 Function: gets one gene out of the db with or without supporting evidence
 Example : $obj->get_Gene('ENSG00000009151','evidence')
 Returns : gene object (with transcripts, exons and supp.evidence if wanted)
 Args    : gene id and supporting tag (if latter not specified, assumes without
	   Note that it is much faster to get genes without supp.evidence!


=cut

sub get_Gene {
   my ($self,$geneid, $supporting) = @_;

   $self->warn("Obj->get_Gene is a deprecated method!\nCalling gene_Obj->get instead!");

   return $self->gene_Obj->get($geneid,$supporting);
}

=head2 get_Gene_by_Transcript_id

 Title   : get_Gene_by_Transcript_id
 Usage   : $obj->get_Gene_by_Transcript_id($transid, $supporting)
 Function: gets one gene out of the db with or without supporting evidence
 Example : $obj->get_Gene_by_Transcript_id('ENST00000009151')
 Returns : gene object (with transcripts, exons and supp.evidence if wanted)
 Args    : transcript id and supporting tag (if latter not specified, assumes without
	   Note that it is much faster to get genes without supp.evidence!


=cut

sub get_Gene_by_Transcript_id {
   my ($self,$tid, $supporting) = @_;

   $self->warn("Obj->get_Gene_by_Transcript_id is a deprecated method! 
Calling gene_Obj->get instead!");

   return $self->gene_Obj->get_Gene_by_Transcript_id($tid,$supporting);
}



=head2 get_Gene_by_DBLink

 Title   : get_Gene_by_DBLink
 Usage   : $obj->get_Gene_by_DBLink($ext_id, $supporting)
 Function: gets one gene out of the db with or without supporting evidence
 Example : $obj->get_Gene_by_DBLink( 'MC1R')
 Returns : gene object (with transcripts, exons and supp.evidence if wanted)
 Args    : transcript id and supporting tag (if latter not specified, assumes without
	   Note that it is much faster to get genes without supp.evidence!


=cut

sub get_Gene_by_DBLink {
   my ($self,$ext_id, $supporting) = @_;

  # $self->warn("Obj->get_Gene_by_DBLink is a deprecated method! Calling gene_Obj->get instead!");

   return $self->gene_Obj->get_Gene_by_DBLink($ext_id,$supporting);
}


=head2 get_Gene_array

 Title   : get_Gene_array
 Usage   :
 Function: old deprecated method, points to new method
           get_gene_array_supporting without asking for supp.evidence
 Example :
 Returns : 
 Args    :


=cut

sub get_Gene_array {
    my ($self,@geneid) = @_;

    $self->throw("Very deprecated method, should call methods with supporting evidence and from
gene_Obj!");
}

=head2 get_Gene_array_supporting

 Title   : get_Gene_array_supporting
 Usage   : $obj->get_Gene_array_supporting($supporting,@geneid)
 Function: Gets an array of genes, with transcripts and exons. If $supporting
           equal to 'evidence' the supporting evidence for each exon is also read
           from the supporting evidence table
 Example : $obj->get_Gene_array_supporting ('evidence',@geneid)
 Returns : an array of gene objects
 Args    : 'evidence' and gene id array

=cut

sub get_Gene_array_supporting {
    my ($self,$supporting,@geneid) = @_;

    $self->warn("Obj->get_Gene_array_supporting is a deprecated method!
Calling gene_Obj->get_array_supporting instead!");

    return $self->gene_Obj->get_array_supporting($supporting,@geneid);
}



=head2 get_Virtual_Contig
    
 Title   : get_Virtual_Contig
 Usage   : $obj->get_Virtual_Contig($transcript,$max_length)
 Function: Gets a Bio::EnsEMBL::DB::Virtual Contig object which 
           spans the whole sequence on which the given 
           Bio::EnsEMBL::Transcript object lies, as long 
           as its length does not exceed max_length. If max_length
           is exceeded, undef is returned instead.
 Example : $obj->get_Virtual_Contig($transcript,50000)
 Returns : VirtualContig Object (or undef)
 Args    : Bio::EnsEMBL::Transcript object and max_length int variable


=cut



sub get_Virtual_Contig_by_Transcript_id {
   my ($self,$tid, $maxlen) = @_;

#   $self->warn("Obj->get_Virtual_contig is a deprecated method! 
#Calling gene_Obj->get_Virtual_contig instead!");

   return my $vc =$self->gene_Obj->get_Virtual_Contig($tid,$maxlen);
}





=head2 get_Transcript_in_VC_coordinates
    
 Title   : get_Transcript_in_VC_coordinates
 Usage   : $obj->get_Transcript_in_VC_coordinates($transcript_id)
 Function: Gets a Bio::EnsEMBL::Transcript object in vc coordinates
 Example : $obj->get_Virtual_Contig($transcript_id)
 Returns : Bio::EnsEMBL::Transcript
 Args    : transcript id


=cut




sub get_Transcript_in_VC_coordinates{
   my ($self,$tid) = @_;

   $self->warn("Obj->get_Transcript_in_VC_coordinates is a deprecated method! 
Calling gene_Obj->get_Virtual_contig instead!");

   return my $transcript =$self->gene_Obj->get_Transcript_in_VC_coordinates($tid);
}

=head2 donor_locator
    
 Title   : get_donor_locator
 Usage   : $obj->get_donor_locator; 
 Function: Reads the meta table of the database to get the donor_database_locator
 Example : get_donor_locator
 Returns : locator string
 Args    : none


=cut

sub get_donor_locator {
    my ($self) = @_;

    $self->warn("Obj->get_donor_locator is a deprecated method! 
Calling Update_Obj->get_donor_locator instead!");
    
    return $self->get_Update_Obj->get_donor_locator();
}

=head2 get_last_update_offset

 Title   : get_last_update_offset
 Usage   : $obj->get_last_update_offset; 
 Function: Reads the meta table of the database to get the last_update time - offset time
 Example : get_last_update_offset
 Returns : UNIX TIME of last update - offset time
 Args    : none

=cut

sub get_last_update_offset{
    my ($self) = @_;

    $self->warn("Obj->get_last_update_offset is a deprecated method! 
Calling Update_Obj->get_last_update_offset instead!");
 
    return $self->get_Update_Obj->get_last_update_offset();
}    

=head2 get_last_update

 Title   : get_last_update
 Usage   : $obj->get_last_update; 
 Function: Reads the db_update table of the database to get the finishing time of the
           last complete update
 Example : get_last_update
 Returns : UNIX TIME of last update
 Args    : none

=cut

sub get_last_update{
    my ($self) = @_;
    
    $self->warn("Obj->get_last_update is a deprecated method! 
Calling Update_Obj->get_last_update_offset instead!");
    
    return $self->get_Update_Obj->get_last_update_offset();
}     

=head2 get_now_offset

 Title   : get_now_offset
 Usage   : $obj->get_now_minus_offset; 
 Function: Gets the current time from the point of view of the database, substracts the
           offset time found in the meta table and gives back unix time of now-offset
 Example : get_now_offset
 Returns : UNIX TIME of now - offset_time
 Args    : none


=cut

sub get_now_offset{
    my ($self) = @_;

    $self->warn("Obj->get_now_offset is a deprecated method! 
Calling Update_Obj->get_now_offset instead!");
   
    return $self->get_Update_Obj->get_now_offset();
}

=head2 get_offset

 Title   : get_offset
 Usage   : $obj->get_offset; 
 Function: Gets the offset time found in the meta table
 Example : get_offset
 Returns : UNIX TIME of offset_time
 Args    : none


=cut

sub get_offset{
    my ($self) = @_;

    $self->throw("Obj->get_offset should not be needed any more!"); 
}
    

=head2 get_Protein_annseq

 Title   : get_Protein_annseq
 Usage   : get_Protein_annseq ($ENSP); 
 Function: Creates an annseq object for a particular peptide, storing the peptide
           sequence in $annseq->primary_seq, and adding all the protein features as generic
           Seqfeatures
 Example : 
 Returns : $annseq
 Args    : $ENSP


=cut

sub get_Protein_annseq{
    my ($self,$ENSP) = @_;

    $self->warn("Obj->get_Protein_annseq is a deprecated method! 
Calling Feature_Obj->get_Protein_annseq instead!");
    
    $self->get_Feature_Obj->get_Protein_annseq($ENSP);
} 

=head2 get_Transcript
    
 Title   : get_Transcript
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut
    
sub get_Transcript{
    my ($self,$transid) = @_;
 
    $self->warn("Obj->get_Transcript is a deprecated method! 
Calling gene_Obj->get_Transcript instead!");

    return $self->gene_Obj->get_Transcript($transid);
}

=head2 get_Translation

 Title   : get_Translation
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Translation{
   my ($self,$translation_id) = @_;

   $self->warn("Obj->get_Translation is a deprecated method! 
Calling gene_Obj->get_Translation instead!");

   return $self->gene_Obj->get_Translation($translation_id);
}

=head2 get_Exon

 Title   : get_Exon
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Exon{
   my ($self,$exonid) = @_;

   $self->warn("Obj->get_Exon is a deprecated method! 
Calling gene_Obj->get_Exon instead!");

   return $self->gene_Obj->get_Exon($exonid);
}

=head2 get_all_Gene_id

 Title   : get_all_Gene_id
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_Gene_id{
   my ($self) = @_;

   $self->warn("Obj->get_all_Gene_id is a deprecated method! 
Calling gene_Obj->get_all_Gene_id instead!");

   return $self->gene_Obj->get_all_Gene_id();
}



=head2 get_all_Transcript_id

 Title   : get_all_Transcript_id
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_Transcript_id{
   my ($self) = @_;

   $self->warn("Obj->get_all_Transcript_id is a deprecated method! 
Calling gene_Obj->get_all_Gene_id instead!");

   return $self->gene_Obj->get_all_Transcript_id();
}

=head2 delete_Exon

 Title   : delete_Exon
 Usage   : $obj->delete_Exon($exon_id)
 Function: Deletes exon, including exon_transcript rows
 Example : $obj->delete_Exon(ENSE000034)
 Returns : nothing
 Args    : $exon_id


=cut

sub delete_Exon{
    my ($self,$exon_id) = @_;

    $self->warn("Obj->delete_Exon is a deprecated method
Calling gene_Obj->delete_Exon instead!");

    return $self->gene_Obj->delete_Exon($exon_id);
}

=head2 delete_Supporting_Evidence

 Title   : delete_Supporting_Evidence
 Usage   : $obj->delete_Supporting_Evidence($exon_id)
 Function: Deletes exon\'s supporting evidence entries
 Example : $obj->delete_Supporting_Evidence(ENSE000034)
 Returns : nothing
 Args    : $exon_id


=cut

sub delete_Supporting_Evidence {
    my ($self,$exon_id) = @_;
 
    $self->warn("Obj->delete_Supporting_Evidence is a deprecated method
Calling gene_Obj->delete_Supporting_Evidence instead!");

    return $self->gene_Obj->delete_Supporting_Evidence($exon_id);
}

=head2 delete_Features

 Title   : delete_Features
 Usage   :
 Function: deletes all features from a contig;
 Example :
 Returns : 
 Args    :


=cut

sub delete_Features {
    my ($self,$contig) = @_;

    $self->warn("Obj->delete_Features is a deprecated method! 
Calling Feature_Obj->delete instead!");

   $self->get_Feature_Obj->delete($contig);
} 

=head2 delete_Gene

 Title   : delete_Gene
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub delete_Gene{
   my ($self,$geneid) = @_;

   $self->warn("Obj->delete_Gene is a deprecated method! 
Calling gene_Obj->delete instead!");

   return $self->gene_Obj->delete($geneid);
}

=head2 geneid_to_cloneid

 Title   : geneid_to_cloneid
 Usage   : @cloneid = $db->geneid_to_cloneid($geneid);
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub geneid_to_cloneid{
    my ($self,$geneid) = @_;
    
    $self->throw("Obj->geneid_to_cloneid is a deprecated method, called gene_Obj->each_cloneid instead!
All the gene, transcript, and exon methods are now to be found in gene_Obj");
    return $self->gene_Obj->each_cloneid($geneid);
}

=head2 write_Gene

 Title   : write_Gene
 Usage   : $obj->write_Gene($gene)
 Function: writes a particular gene into the database
           
 Example :
 Returns : 
 Args    :


=cut


sub write_Gene{
   my ($self,$gene) = @_;

   my ( $p, $f, $l ) = caller;
   $self->warn("Obj->write_Gene is a deprecated method! Calling from $p::$f::$l\n" );

   return $self->get_GeneAdaptor()->store( $gene );
}

=head2 write_all_Protein_features

 Title   : write_all_Protein_features
 Usage   : $obj->write_all_Protein_features($ENSP)
 Function: writes all protein features of a particular peptide into the database          
 Example :
 Returns : 
 Args    :


=cut

sub write_all_Protein_features {
    my ($self,$prot_annseq,$ENSP) = @_;

    $self->warn("Obj->write_all_Protein_features is a deprecated method! 
Calling Feature_Obj->write_all_Protein_features instead!");
    
    $self->get_Feature_Obj->write_all_Protein_features($prot_annseq,$ENSP);
} 


=head2 write_Protein_feature

 Title   : write_Protein_feature
 Usage   : $obj->write_Protein_feature($ENSP, $feature)
 Function: writes a protein feature object of a particular peptide into the database          
 Example :
 Returns : 
 Args    :


=cut

sub write_Protein_feature {
    my ($self,$ENSP,$feature) = @_;
 
    $self->warn("Obj->write_Protein_feature is a deprecated method! 
Calling Feature_Obj->write_Protein_feature instead!");
    
    $self->get_Feature_Obj->write_Protein_feature($ENSP,$feature);
} 

=head2 write_Feature

 Title   : write_Feature
 Usage   : $obj->write_Feature($contig,@features)
 Function: Writes a feature on the genomic sequence of a contig into the database
 Example :
 Returns : nothing
 Args    : Bio::EnsEMBL::SeqFeatureI


=cut

sub write_Feature {
    my ($self,$contig,@features) = @_;

    $self->warn("Obj->write_Feature is a deprecated method! 
Calling Feature_Obj->write_Feature instead!");
    
    $self->get_Feature_Obj->write($contig,@features);
} 

=head2 write_supporting_evidence

 Title   : write_supporting_evidence
 Usage   : $obj->write_supporting_evidence
 Function: Writes supporting evidence features to the database
 Example :
 Returns : nothing
 Args    : None


=cut

sub write_supporting_evidence {
    my ($self,$exon) = @_;

    $self->warn("Obj->write_supporting_evidence is a deprecated method!
Calling gene_Obj->write_supporting_evidence instead!");

    return $self->gene_Obj->write_supporting_evidence($exon);
}

=head2 get_supporting_evidence

 Title   : get_supporting_evidence
 Usage   : $obj->get_supporting_evidence
 Function: Writes supporting evidence features to the database
 Example :
 Returns : nothing
 Args    : array of exon objects, needed to know which exon to attach the evidence to


=cut

sub get_supporting_evidence {
    my ($self,@exons) = @_;

    $self->warn("Obj->get_supporting_evidence is a deprecated method! 
Calling gene_Obj->get_supporting_evidence instead!");

   return $self->gene_Obj->get_supporting_evidence(@exons);
}

=head2 write_Analysis

 Title   : write_Analysis
 Usage   : $obj->write_Analysis($anal)
 Function: Writes analysis details to the database
           Checks first whether this analysis entry already exists
 Example :
 Returns : int
 Args    : Bio::EnsEMBL::AnalysisI


=cut

sub write_Analysis {
    my ($self,$analysis) = @_;

    $self->warn("Obj->write_Analysis is a deprecated method! 
Calling Feature_Obj->write_Analysis instead!");
    
    $self->get_Feature_Obj->write_Analysis($analysis);
} 
    
=head2 exists_Homol_Feature

 Title   : exists_Homol_Feature
 Usage   : $obj->exists_Homol_Feature($feature)
 Function: Tests whether this feature already exists in the database
 Example :
 Returns : nothing
 Args    : Bio::SeqFeature::Homol


=cut

sub exists_Homol_Feature {
    my ($self,$feature,$analysisid,$contig) = @_;

    $self->warn("Obj->exists_Homol_Feature is a deprecated method! 
Calling Feature_Obj->exists_Homol_Feature instead!");
    
    $self->get_Feature_Obj->exists($feature,$analysisid,$contig);
} 
    
=head2 get_Analysis

 Title   : get_Analysis
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Analysis {
    my ($self,$id) = @_;

    $self->warn("Obj->get_Analysis is a deprecated method! 
Calling Feature_Obj->get_Analysis instead!");
    
    $self->get_Feature_Obj->get_Analysis($id);
} 

=head2 exists_Analysis

 Title   : exists_Analysis
 Usage   : $obj->exists_Analysis($anal)
 Function: Tests whether this feature already exists in the database
 Example :
 Returns : Analysis id if the entry exists
 Args    : Bio::EnsEMBL::Analysis


=cut

sub exists_Analysis {
    my ($self,$analysis) = @_;
    
    $self->warn("Obj->exists_Analysis is a deprecated method! 
Calling Feature_Obj->exists_Analysis instead!");
    
    $self->get_Feature_Obj->exists_Analysis($analysis);
} 
 
    
=head2 write_Transcript

 Title   : write_Transcript
 Usage   : $obj->write_Transcript($trans,$gene)
 Function: writes a particular transcript *but not the exons* into
           the database
 Example :
 Returns : 
 Args    :


=cut

sub write_Transcript{
   my ($self,$trans,$gene) = @_;

   $self->warn("Obj->write_Transcript is a deprecated method! 
Calling gene_Obj->write_Transcript instead!");

   return $self->gene_Obj->write_Transcript($trans,$gene);
}

=head2 write_Translation

 Title   : write_Translation
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub write_Translation{
    my ($self,$translation) = @_;

    $self->warn("Obj->write_Translation is a deprecated method
Calling gene_Obj->write_Translation instead!");

    return $self->gene_Obj->write_Translation($translation);
}


=head2 write_Exon

 Title   : write_Exon
 Usage   : $obj->write_Exon($exon)
 Function: writes a particular exon into the database
 Example :
 Returns : 
 Args    :


=cut

sub write_Exon {
   my ($self,$exon) = @_;

   $self->warn("Obj->write_Exon is a deprecated method! 
Calling gene_Obj->write_Exon instead!");

   return $self->gene_Obj->write_Exon($exon);
}




=head2 get_PredictionFeature_as_Transcript

 Title   : get_PredictionFeature_as_Transcript
 Usage   :$obj->get_PredictionFeature_as_Transcript($genscan_id)
 Function:Call get_PredictionFeature_as_Transcript in Feature_obj object,see documentation for this method
 Example :
 Returns : 
 Args    :


=cut

sub get_PredictionFeature_as_Transcript{
   my ($self,$genscan_id) = @_;

   $self->warn("Deprecated method : use FeatureAdaptor->fetch_PredictionFeature_as_Transcript");
   return $self->get_FeatureAdaptor->fetch_PredictionFeature_as_Transcript($genscan_id);

}


=head2 write_Clone

 Title   : write_Clone
 Usage   : DEPRECATED
 Function: DEPRECATED
 Example : DEPRECATED
 Returns : DEPRECATED
 Args    : DEPRECATED


=cut

sub write_Clone {
    my ($self,$clone) = @_;


    $self->warn("this method is being depreciated please use Bio::EnsEMBL::CloneAdaptor->store()\n");

    my $clone_ad = $self->get_CloneAdaptor();

    $clone_ad->store($clone);
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

 Title   : write_Contig
 Usage   : $obj->write_Contig($contig,$clone)
 Function: Writes a contig and its dna into the database
 Example :
 Returns : 
 Args    :


=cut

sub write_Contig {
    my($self, $contig, $clone)  = @_;
       

    $self->warn("this method is depreciated please use Bio::EnsEMBL::DBSQL::RawContigAdaptor->store()\n");
    
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

 Title   : _insertSequence
 Usage   : $obj->_insertSequence
 Function: Insert the dna sequence and date into the dna table.
 Example :
 Returns : 
 Args    : $sequence, $date


=cut

sub _insertSequence {
    my ($self, $sequence, $date) = @_;
    
    $sequence =~ tr/atgcn/ATGCN/;
    
    $self->warn("this method is depreciated please use Bio::EnsEMBL::RawContigAdaptor->_insertSequence\n");
    
    my $rca = $self->get_RawContigAdaptor();

    $rca->_insertSequence($sequence, $date);

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

 Title   : write_Species
 Usage   : $obj->write_Species
 Function: writes a species object into the database
 Example :
 Returns : 
 Args    :


=cut

sub write_Species {
    my ($self,$species) = @_;

    $self->throw("there isn't a species table in the new schema\n");

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

 Title   : get_Update_Obj
 Usage   : DEPRECATED
 Function: DEPRECATED
 Example : DEPRECATED
 Returns : DEPRECATED
 Args    : DEPRECATED

=cut

sub get_Update_Obj {
    my ($self) = @_;

    $self->warn("get_Update_Obj is deprecated. There is no longer an " .
		"Update_Obj object");

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

 Title   : get_all_chr_ids
 Usage   : DEPRECATED    
           Use $dba->get_ChromosomeAdaptor()->fetch_all() instead.
 Function: DEPRECATED returns all the valid FPC contigs from given golden path
 Example : DEPRECATED
 Returns : DEPRECATED
 Args    : DEPRECATED static golden path type (typically, 'UCSC')

=cut

sub get_all_chr_ids {
   my ($self, $type) = @_;

   $self->warn("DBAdaptor->get_all_chr_ids is deprecated\n" .
	       "Use \$dba->get_ChromosomeAdaptor()->fetch_all() instead");

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

 Title   : get_all_fpcctg_ids
 Usage   : DEPRECATED
           Use $dba->get_StaticGoldenPathAdaptor()->get_all_fpc_ids() instead
 Function: DEPRECATED returns all the valid FPC contigs from given golden path
 Example : DEPRECATED
 Returns : DEPRECATED
 Args    : DEPRECATED static golden path type (typically, 'UCSC')


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

 Title   : get_object_by_wildcard
 Usage   : DEPRECATED
 Function: DEPRECATED
 Example : DEPRECATED
 Returns : DEPRECATED
 Args    : DEPRECATED


=cut

sub get_object_by_wildcard{
   my ($self,$type,$string) = @_;

   $self->warn("DBAdaptor->get_object_by_wildcard is deprecated. Use " . 
	       "specific adaptor accessor functions instead");

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

 Title   : write_Chromosome
 Usage   : DEPRECATED
 Function: DEPRECATED writes a chromosome into the database
 Example : DEPRECATED
 Returns : DEPRECATED
 Args    : DEPRECATED


=cut

sub write_Chromosome {
  my ($self,$chromosome,$length, $known_genes, $unknown_genes, $snps) = @_;

  $self->throw("DBAdaptor->write_Chromosome() is deprecated.  " .
	       "No replacement has been written yet. If a replacement " .
	       "is created it will be within ChromosomeAdaptor");

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

 Title   : _analysis_cache
 Usage   : DEPRECATED $obj->_analysis_cache()
 Function: DEPRECATED
 Returns : DEPRECATED reference to a hash
 Args    : DEPRECATED newvalue (optional)


=cut

sub _analysis_cache{
   my $self = shift;

   $self->throw("DBAdaptor->_analysis_cache is deprecated\n");

#   if( @_ ) {
#      my $value = shift;
#      $obj->{'_analysis_cache'} = $value;
#    }
#    return $obj->{'_analysis_cache'};
}

=head2 _contig_seq_cache

 Title   : _contig_seq_cache
 Usage   : DEPRECATED
 Function: DEPRECATED
 Example : DEPRECATED
 Returns : DEPRECATED
 Args    : DEPRECATED


=cut

sub _contig_seq_cache{
   my ($self,$id,$seq) = @_;

   $self->throw("DBAdaptor->_contig_seq_cache is deprecated\n");


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

 Title   : _flush_seq_cache
 Usage   : DEPRECATED
 Function: DEPRECATED
 Example : DEPRECATED
 Returns : DEPRECATED
 Args    : DEPRECATED

=cut

sub _flush_seq_cache{
   my ($self,@args) = @_;

   $self->throw("DBAdaptor->_flush_seq_cache is deprecated\n");


   $self->{'_contig_seq_cache'} = {};
}



=head2 _lock_tables

 Title   : _lock_tables
 Usage   : DEPRECATED
 Function: DEPRECATED
 Example : DEPRECATED
 Returns : DEPRECATED
 Args    : DEPRECATED


=cut

sub _lock_tables{
   my ($self,@tables) = @_;
   
   $self->warn("DBAdaptor->_lock_tables is deprecated\n");

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

 Title   : _unlock_tables
 Usage   : DEPRECATED
 Function: DEPRECATED
 Example : DEPRECATED
 Returns : DEPRECATED
 Args    : DEPRECATED


=cut

sub _unlock_tables{
   my ($self,@tables) = @_;

   $self->warn("DBAdaptor->_unlock_tables is deprecated\n");


   my $sth = $self->prepare("unlock tables");
   my $rv  = $sth->execute();
   $self->throw("Failed to unlock tables") unless $rv;
   %{$self->{'_lock_table_hash'}} = ();
}


=head2 get_Clone

 Title   : get_Clone
 Usage   : DEPRECATED 
           Use $db->get_CloneAdaptor()->fetch_by_accession_version($acc,$ver) 
           or $db->get_CloneAdaptor()->fetch_by_accession($acc) instead
 Function: DEPRECATED retrieve latest version of a clone from the database
 Example : DEPRECATED
 Returns : DEPRECATED
 Args    : DEPRECATED

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

 Title   : list_embl_version_by_Clone
 Usage   : DEPRECATED 
           use $db->get_CloneAdaptor()->list_embl_version_by_accession instead
 Function: DEPRECATED 
           retrieve list of embl_versions of a clone from the database
 Example : DEPRECATED
           @versions = $dbobj->list_embl_versions_by_Clone('AB000381');
 Returns : DEPRECATED @versions
 Args    : DEPRECATED $accession

=cut

sub list_embl_version_by_Clone { 
    my( $self, $accession ) = @_;

    $self->warn('DBAdaptor->list_embl_version_by_Clone() is deprecated.\n' .
	'Use \$db->get_CloneAdaptor()->list_embl_version_by_accession(\$acc)');

    my $ca = $self->get_CloneAdaptor;

    return $ca->list_embl_version_by_accession($accession);
}

=head2 get_Clone_by_version

 Title   : get_Clone_by_version
 Usage   : DEPRECATED
 Function: DEPRECATED retrieve specific version of a clone from the database
 Example : DEPRECATED
 Returns : DEPRECATED
 Args    : DEPRECATED


=cut

sub get_Clone_by_version { 
    my ($self,$accession,$ver) = @_;

    $self->warn('DBAdaptor->get_CLone_by_version() is deprecated.\n' .
      'Use \$db->get_CloneAdaptor()->fetch_by_accession_version(\$acc,\$ver)');

    my $ca = $self->get_CloneAdaptor;

    return $ca->fetch_by_accession_version($accession,$ver);
}
  
=head2 get_Contig

 Title   : get_Contig
 Usage   : DEPRECATED
 Function: DEPRECATED
 Example : DEPRECATED
 Returns : DEPRECATED
 Args    : DEPRECATED

=cut

sub get_Contig {
    my ($self,$id) = @_;

    $self->warn("DBAdaptor->get_Contig() is a deprecated.\n" .
	       "Use \$db->get_RawContigAdaptor()->fetch_by_name(\$id)"); 

    return $self->get_RawContigAdaptor->fetch_by_name($id);
}


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

  
 
sub get_Contig_by_international_id{
   my ($self,$int_id) = @_;
   $self->throw("DBADaptor->get_Contig_by_international_id() is deprecated\n" .
		"No replacement has been implemented\n");

   return undef;

#   my $sth=$self->prepare("select id from contig where international_id = '$int_id'");
#   $sth->execute;
#   my $row = $sth->fetchrow_hashref;
#   my $id  = $row->{'id'};

#   return $self->get_Contig($id);
}

=head2 get_Contigs_by_Chromosome

 Title   : get_Contig_by_Chromosome
 Usage   : DEPRECATED @contigs = $dbobj->get_Contig_by_Chromosome( $chrObj );
 Function: DEPRECATED
           retrieve contigs belonging to a certain chromosome from the
           database 
 Example : DEPRECATED
 Returns : DEPRECATED A list of Contig objects. Probably an empty list.
 Args    : DEPRECATED


=cut

sub get_Contigs_by_Chromosome {
    my ($self,$chromosome ) = @_;
    
    $self->throw("Obj->get_Contigs_by_Chromosome is a deprecated method! 
Call Contig->get_by_Chromosome instead!");
}

=head2 get_all_Clone_id

 Title   : get_all_Clone_id
 Usage   : DEPRECATED @cloneid = $obj->get_all_Clone_id
 Function: DEPRECATED returns all the valid (live) Clone ids in the database
 Example : DEPRECATED
 Returns : DEPRECATED
 Args    : DEPRECATED


=cut

sub get_all_Clone_id{
   my ($self) = @_;
   my @out;

   my $sth = $self->prepare("select id from clone");
   my $res = $sth->execute;

   while( my $rowhash = $sth->fetchrow_hashref) {
       push(@out,$rowhash->{'id'});
   }

   return @out;
}

=head2 get_all_Contig_id

 Title   : get_all_Contig_id
 Usage   : DEPRECATED @Contigid = $obj->get_all_Contig_id
           call $db->get_RawContigAdaptor->fetch_all() and iterate through
           returned raw contig objects instead,
 Function: DEPRECATED returns all the valid (live) Contig ids in the database
 Example : DEPRECATED
 Returns : DEPRECATED
 Args    : DEPRECATED


=cut

sub get_all_Contig_id{
   my ($self) = @_;
   my @out;

   $self->warn('DBAdaptor->get_all_Contig_id is deprecated. Use:
   \$rca = \$self->get_RawContigAdaptor();
   foreach \$contig (\$rca->fetch_all()) {
     push \@out, \$contig->name();
   }
   ');
   
   my $rca = $self->get_RawContigAdaptor();

   foreach my $contig ($rca->fetch_all()) {
     push @out, $contig->name();
   }

   return @out;
}


=head2 perl_only_sequences

 Title   : perl_only_sequences
 Usage   : DEPRECATED $obj->perl_only_sequences($newval)
 Function: DEPRECATED
 Returns : DEPRECATED value of perl_only_sequences
 Args    : DEPRECATED newvalue (optional)


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

 Title   : perl_only_contigs
 Usage   : DEPRECATED $obj->perl_only_contigs($newval)
 Function: DEPRECATED
 Returns : DEPRECATED value of perl_only_contigs
 Args    : DEPRECATED newvalue (optional)


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

 Title   : _crossdb
 Usage   : DEPRECATED $obj->_crossdb($newval)
 Function: DEPRECATED
 Returns : DEPRECATED value of _crossdb
 Args    : DEPRECATED newvalue (optional)

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


sub create_tables { 
  my $self = shift;

  $self->throw('DBAdaptor->create_tables is deprecated.\n' .
	       'No replacement has been written');

  # get all adaptors ...
  # call create_tables on them

  # create some tables without adaptors
  # (which should disappear once)
}

#=head2 get_FamilyAdaptor

# Title   : get_FamilyAdaptor
# Usage   :
# Function:
# Example :
# Returns : 
# Args    :


#=cut

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
    
    my( $fa );
    unless ($fa = $self->{'_externaldata_family_familyadaptor'}) {
        eval{
            require Bio::EnsEMBL::ExternalData::Family::FamilyAdaptor;
        };
        if ($@) {
            $self->throw(
                "Unable to load 'Bio::EnsEMBL::ExternalData::Family::FamilyAdaptor'\n"
                . "It is not part of the core Ensembl distribution.\n"
                . "Have you installed it?");
        }
        $fa = Bio::EnsEMBL::ExternalData::Family::FamilyAdaptor->new($self);
        $self->{'_externaldata_family_familyadaptor'} = $fa;
    }
    return $fa;
}



=head2 find_GenomeHits
    
 Title   : find_GenomeHits
 Usage   : DEPRECATED my @features = $self->find_GenomeHits($hid)
 Function: DEPRECATED Finds all features in the db that
           are hits to a sequence with id $hid
 Example : DEPRECATED
 Returns : DEPRECATED @ Bio::EnsEMBL::FeaturePair
 Args    : DEPRECATED string

=cut
 
sub find_GenomeHits {
    my ($self,$arg) = @_;
    $self->throw("Bio::EnsEMBL::DBSQL::DBAdaptor::find_GenomeHits is deprecated");
    return ();
    #    return $self->feature_Obj->find_GenomeHits($arg);
}


=head2 release_number

 Title   : DEPRECATED release_number
 Usage   : DEPRECATED 
 Function: DEPRECATED 
  #######SNEAKY METHOD FOR RELEASE NUMBER, VERY TEMPORAR%Y!!!!
 Example : DEPRECATED
 Returns : DEPRECATED
 Args    : DEPRECATED


=cut

sub release_number{
   my ($self,@args) = @_;

   $self->throw("DBAdaptor::release_number is deprecated\n"); 

   return 110;
}

1;
