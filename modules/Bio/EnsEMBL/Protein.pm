
#
# BioPerl module for Protein.pm
#
# Cared for by Emmanuel Mongin <mongin@ebi.ac.uk>
#
# Copyright Emmanuel Mongin
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Protein.pm - DESCRIPTION of Object

=head1 SYNOPSIS

my $protein = Bio::EnsEMBL::Protein->new ( -seq =>$sequence,
					      -accession_number  => $id,
					      -display_id => $id,
					      -primary_id => $id,
					      -id => $id,
					      -desc => $desc,
					      -moltype => $moltype,
					      );


=head1 DESCRIPTION

This object inherit from Bio::SeqI for the implementation and from PrimarySeq
A protein object hold the basic PrimarySeq requirements, hold annotation, DBlink, protein_features object.
It will also soon hold family and SNPs objects. The DB adaptor for this object is called Protein_Adaptor 
Most of the calls should be done from the protein object, the protein object will then deal with calling the rigth DBadaptor (ProteinAdaptor, ProteinFeatureAdaptor, ...)

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Protein;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Root;

use Bio::SeqI;
use Bio::EnsEMBL::Transcript;
use Bio::SeqIO;
use Bio::Tools::SeqStats;

# Object preamble - inheriets from Bio::Root::Object

@ISA = qw( Bio::PrimarySeq Bio::SeqI);

=head2 SeqI interface features
=cut



=head2 adaptor

 Title   : adaptor
 Usage   : $self->adaptor($adaptor)
 Function: Getter/Setter for the ProteinAdaptor for this object
 Example : my $db = $prot->adaptor()->db();
 Returns : a ProteinAdaptor object
 Args    : (optional) a new ProteinAdaptor

=cut

sub adaptor {
  my ($self, $adaptor) = @_;

  if(defined $adaptor) {
    $self->{_adaptor} = $adaptor;
  }
  return $self->{_adaptor};
}



=head2 species

 Title   : species
 Usage   : $obj->species($newval)
 Function: 
 Returns : value of species
 Args    : newvalue (optional)


=cut

sub species{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'species'} = $value;
    }
    return $obj->{'species'};

}

=head2 primary_seq

 Title   : primary_seq
 Usage   : $obj->primary_seq($newval)
 Function: 
 Returns : value of primary_seq
 Args    : newvalue (optional)


=cut

sub primary_seq{
    my ($self) = @_;   
    return $self;
}



=head2 gene

 Title   : gene
 Usage   : my $gene = $prot->gene()
 Function: Getter/Setter for the gene associated with this protein
 Returns : Bio::EnsEMBL::Gene
 Args    : newvalue (optional)

=cut

sub gene {
  my ($self, $gene) = @_;

  if(defined $gene) {
    unless($gene->isa('Bio::EnsEMBL::Gene')) {
      $self->throw("$gene is not a valid Bio::EnsEMBL::Gene object");
    }
    $self->{'_gene'} = $gene;
  }

  return $self->{'_gene'};
}


=head2 transcript

 Title   : transcript
 Usage   : my $transcript = $prot->transcript()
 Function: Getter/Setter for the gene associated with this protein
 Returns : Bio::EnsEMBL::Transcript
 Args    : newvalue (optional)

=cut

sub transcript {
  my ($self, $transcript) = @_;

  if(defined $transcript) {
    unless($transcript->isa('Bio::EnsEMBL::Transcript')) {
      $self->throw("$transcript is not a valid" .
		   "Bio::EnsEMBL::Transcript object");
    }
    $self->{'_transcript'} = $transcript;
  }

  return $self->{'_transcript'};
}


=head2 molecule

 Title   : molecule
 Usage   : $obj->molecule($newval)
 Function: 
 Returns : value of molecule
 Args    : newvalue (optional)


=cut

sub molecule{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'molecule'} = $value;
    }
    return $obj->{'molecule'};

}

=head2 moltype

 Title   : moltype
 Usage   : $obj->moltype($newval)
 Function: 
 Returns : value of moltype
 Args    : newvalue (optional)


=cut

sub moltype{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'moltype'} = $value;
    }
    return $obj->{'moltype'};

}

=head2 top_SeqFeatures
    
 Title   : top_SeqFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut
    
sub top_SeqFeatures {
    my ($self,@args) = @_;

    return @{$self->get_all_ProteinFeatures};
}




sub get_all_ProteinFeatures {
  my $self = shift;

  my @f = ();

  push(@f, @{$self->get_all_DomainFeatures()});
    
  push(@f, @{$self->get_all_blastpFeatures()});
    
  push(@f, @{$self->get_all_SigpFeatures()});
  
  push(@f, @{$self->get_all_TransmembraneFeatures()});
  
  push(@f, @{$self->get_all_CoilsFeatures()});
  
  push(@f, @{$self->get_all_LowcomplFeatures()});

  return \@f;
}




=head2 get_Family

 Title   : get_Family
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_Families{
  my ($self) = @_;

  unless(defined ($self->{'_families'})) {
    $self->{'_families'} = [];
    my $proteinid = $self->id();
    my $fa = $self->adaptor()->db()->get_FamilyAdaptor();

    #should this return multiple families?
    my $family = $fa->get_Family_of_Ensembl_pep_id($proteinid);
    $self->add_Family($family);
   }

  return $self->{'_families'};
}



=head2 add_Family

 Title   : add_Family
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_Family{
  my ($self,$value) = @_;
  
  unless(defined $value && 
	 $value->isa('Bio::EnsEMBL::ExternalData::Family::Family')) {
    $self->throw("[$value] is not a Family object");
  }

  push(@{$self->{'_family'}},$value);  
}

=head2 Protein Specific Features
=cut

=head2 get_all_DomainFeatures

 Title   : get_all_DomainFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_DomainFeatures{
 my ($self) = @_;

 if (defined ($self->{'_domains'})) {
   return $self->{'_domains'};
 }

 my @f = ();
 push(@f,@{$self->get_all_PrintsFeatures()});
 push(@f,@{$self->get_all_PfamFeatures()});
 push(@f,@{$self->get_all_PrositeFeatures()});
 push(@f,@{$self->get_all_SuperfamilyFeatures()});  
 push(@f,@{$self->get_all_ProfileFeatures()});

 return \@f;
}


=head2 get_all_ProfileFeatures

 Title   : get_all_ProfileFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_ProfileFeatures{
  my ($self) = @_;
  
  unless(defined ($self->{'_profile'})) {
    $self->{'_profile'} = [];
    my $proteinid = $self->id();
    my $pfa = $self->adaptor->db()->get_ProteinFeatureAdaptor;
    my $array_features = 
      $pfa->fetch_all_by_feature_and_dbID('PROFILE',$proteinid);
    foreach my $in (@$array_features) {
      $self->add_Profile($in);
    }
  }
  
  return ( $self->{'_profile'} || [] );    
}

=head2 add_Profile

 Title   : add_Profile
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_Profile{
  my ($self,$value) = @_;
  
  if ((!defined $value) || (!$value->isa('Bio::EnsEMBL::ProteinFeature'))) {
    $self->throw("The Protein Feature added is not defined or is not a "
		 ."protein feature object");
  }

  push(@{$self->{'_profile'}},$value); 
  
}

=head2 get_all_blastpFeatures 

 Title   : get_all_blastpFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_blastpFeatures{
  my ($self) = @_;
  
  unless(defined ($self->{'_blastp'})) {
    $self->{'_blastp'} = [];
    my $proteinid = $self->id();
    my $pfa = $self->adaptor->db()->get_ProteinFeatureAdaptor();
    my $array_features = 
      $pfa->fetch_all_by_feature_and_dbID('blastp',$proteinid);
    foreach my $in (@$array_features) {
      $self->add_blastp($in);
    }
  }

  return ( $self->{'_blastp'} || [] );
}

=head2 add_blastp

 Title   : add_blastp
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_blastp{
  my ($self,$value) = @_;
        
  if ((!defined $value) || (!$value->isa('Bio::EnsEMBL::ProteinFeature'))) {
    $self->throw("[$value] is not defined or is not a protein feature object");
  }

  push(@{$self->{'_blastp'}},$value); 
}

=head2 get_all_PrintsFeatures

 Title   : get_all_PrintsFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=cut

sub get_all_PrintsFeatures{
  my ($self) = @_;

  unless(defined ($self->{'_prints'})) {
    $self->{'_prints'} = [];	# init in case there are no PRINTS
    my $proteinid = $self->id();
    my $pfa = $self->adaptor()->db()->get_ProteinFeatureAdaptor();
    my $array_features = 
      $pfa->fetch_all_by_feature_and_dbID('PRINTS',$proteinid);
    foreach my $in (@$array_features) {
      $self->add_Prints($in);
    } 
  }
  return ( $self->{'_prints'} || [] );
}

=head2 add_Prints

 Title   : add_Prints
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_Prints{
    my ($self,$value) = @_;
        
    if ((!defined $value) || (!$value->isa('Bio::EnsEMBL::ProteinFeature'))) {
	$self->throw("[$value] is not a protein feature object");
    }

   push(@{$self->{'_prints'}},$value); 
}


=head2 get_all_PfamFeatures

 Title   : get_all_PfamFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_PfamFeatures{
  my ($self) = @_;
  
  unless($self->{'_pfam'}) {
    $self->{'_pfam'} = [];
    my $proteinid = $self->id();
    my $pfa = $self->adaptor()->db()->get_ProteinFeatureAdaptor();
    my $array_features = 
      $pfa->fetch_all_by_feature_and_dbID('Pfam',$proteinid);
    foreach my $in (@$array_features) {
      $self->add_Pfam($in);
    }
  }
  return ( $self->{'_pfam'} || [] );
}

=head2 add_Pfam

 Title   : add_Pfam
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=cut

sub add_Pfam{
 my ($self,$value) = @_;
        
 if ((!defined $value) || (!$value->isa('Bio::EnsEMBL::ProteinFeature'))) {
   $self->throw("[$value] is not a protein feature object");
 }

 push(@{$self->{'_pfam'}},$value); 
}


=head2 get_all_PrositeFeatures

 Title   : get_all_PrositeFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=cut

sub get_all_PrositeFeatures{
  my ($self) = @_;
  
  unless(defined ($self->{'_prosite'})) {
    $self->{'_prosite'} = [];
    my $proteinid = $self->id();
    my $pfa = $self->adaptor()->db()->get_ProteinFeatureAdaptor();
    my $array_features = 
      $pfa->fetch_all_by_feature_and_dbID('PROSITE',$proteinid);
    foreach my $in (@$array_features) {
      $self->add_Prosite($in);
    }
  }
  return ( $self->{'_prosite'} || [] );
}

=head2 add_Prosite

 Title   : add_Prosite
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=cut

sub add_Prosite{
  my ($self,$value) = @_;
    
  if ((!defined $value) || (!$value->isa('Bio::EnsEMBL::ProteinFeature'))) {
    $self->throw("[$value] is not a protein feature object");
  }

  push(@{$self->{'_prosite'}},$value); 
}

=head2 get_all_SigpFeatures

 Title   : get_all_SigpFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_SigpFeatures{
  my ($self) = @_;

  unless(defined ($self->{'_sigp'})) {
    $self->{'_sigp'} = [];
    my $proteinid = $self->id();
    my $pfa = $self->adaptor()->db()->get_ProteinFeatureAdaptor();
    my $array_features = 
      $pfa->fetch_all_by_feature_and_dbID('Signalp',$proteinid);
    foreach my $in (@$array_features) {
      $self->add_Sigp($in);
    }
  }
  return ( $self->{'_sigp'} || [] );
}

=head2 add_Sigp

 Title   : add_Sigp
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_Sigp{
  my ($self,$value) = @_;
  
  if ((!defined $value) || (!$value->isa('Bio::EnsEMBL::ProteinFeature'))) {
    $self->throw("[$value] is not defined or is not a protein feature object");
  }
  push(@{$self->{'_sigp'}},$value); 
}


=head2 get_all_TransmembraneFeatures

 Title   : get_all_TransmembraneFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_TransmembraneFeatures{
 my ($self) = @_;

 unless(defined ($self->{'_transmembrane'})) {
   $self->{'_transmembrane'} = []; 
   my $proteinid = $self->id();
   my $pfa = $self->adaptor()->db()->get_ProteinFeatureAdaptor();
   my $array_features = 
     $pfa->fetch_all_by_feature_and_dbID('Tmhmm',$proteinid);
   foreach my $in (@$array_features) {
     $self->add_Transmembrane($in);
   }
 }

 return ( $self->{'_transmembrane'} || [] );    

}

=head2 add_Transmembrane

 Title   : add_Transmembrane
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_Transmembrane{
  my ($self,$value) = @_;
  
  if ((!defined $value) || (!$value->isa('Bio::EnsEMBL::ProteinFeature'))) {
    $self->throw("The Protein Feature added is not defined" .
		 "or is not a protein feature object");
    }
  push(@{$self->{'_transmembrane'}},$value); 
}

=head2 get_all_CoilsFeatures

 Title   : get_all_CoilsFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_CoilsFeatures{
 my ($self) = @_;

 unless(defined ($self->{'_coils'})) {
   $self->{'_coils'} = [];
   my $proteinid = $self->id();
   my $pfa = $self->adaptor()->db()->get_ProteinFeatureAdaptor();
   my $array_features = 
     $pfa->fetch_all_by_feature_and_dbID('ncoils',$proteinid);
   foreach my $in (@$array_features) {
     $self->add_Coils($in);
   }
 }
 return ( $self->{'_coils'} || [] );
}

=head2 add_Coils

 Title   : add_Coils
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_Coils{
  my ($self,$value) = @_;
  
  if ((!defined $value) || (!$value->isa('Bio::EnsEMBL::ProteinFeature'))) {
    $self->throw("[$value] is not a protein feature object");
    }
  push(@{$self->{'_coils'}},$value); 
}

=head2 get_all_LowcomplFeatures

 Title   : get_all_LowcomplFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=cut

sub get_all_LowcomplFeatures{
 my ($self) = @_;
 
 unless (defined ($self->{'_lowcompl'})) {
   my $proteinid = $self->id();
   my $pfa = $self->adaptor()->db()->get_ProteinFeatureAdaptor();
   my $array_features = 
     $pfa->fetch_all_by_feature_and_dbID('Seg',$proteinid);
   foreach my $in (@$array_features) {
     $self->add_Lowcompl($in);
   }
 }
 return ( $self->{'_lowcompl'} || [] );
}

=head2 add_LowCompl

 Title   : add_LowCompl
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_Lowcompl{
  my ($self,$value) = @_;

  if ((!defined $value) || (!$value->isa('Bio::EnsEMBL::ProteinFeature'))) {
    $self->throw("[$value] is not a protein feature object");
  }

  push(@{$self->{'_lowcompl'}},$value); 
}

=head2 get_all_SuperfamilyFeatures

 Title   : get_all_SuperfamilyFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_SuperfamilyFeatures{
  my ($self) = @_;
  
  unless(defined ($self->{'_superfamily'})) {
    $self->{'_superfamily'} = [];
    my $proteinid = $self->id();
    my $pfa = $self->adaptor()->db->get_ProteinFeatureAdaptor();
    my $array_features = 
     $pfa->fetch_all_by_feature_and_dbID('superfamily',$proteinid);
    foreach my $in (@$array_features) {
      $self->add_Superfamily($in);
    }
  }
  return ( $self->{'_superfamily'} || [] );
}

=head2 add_Superfamily

 Title   : add_Superfamily
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_Superfamily{
  my ($self,$value) = @_;
  
  if ((!defined $value) || (!$value->isa('Bio::EnsEMBL::ProteinFeature'))) {
    $self->throw("[$value] is not a protein feature object");
  } 

  push(@{$self->{'_superfamily'}},$value); 
}




=head2 length

 Title   : length
 Usage   :my $lenght = $protein->length
 Function:
 Example :
 Returns :length of the aa sequence 
 Args    :none


=cut

sub length{
   my ($self) = @_;
   my $length = length($self->seq);

   return $length;

}

=head2 get_all_DBlinks

 Title   : get_all_DBlinks
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_DBLinks{
 my ($self) = @_;

 unless(defined $self->{'_dblinks'}) {
   $self->{'_dblinks'} = [];
   my $transcript = $self->transcript();

   unless(defined $transcript) {
     throw("Cannot fetch Protein DBLinks without a valid Transcript.\n");
   }

   my $dbea = $self->adaptor()->db()->get_DBEntryAdaptor();
   $self->{'_dblinks'} = 
     $dbea->fetch_all_by_Translation($transcript->translation);
 }
 
 return $self->{'_dblinks'};
}


=head2 get_all_DASFactories

  Arg [1]   : none
  Function  : Retrieves a listref of registered DAS objects
              TODO: Abstract to a DBLinkContainer obj
  Returntype: [ DAS_objects ]
  Exceptions:
  Caller    :
  Example   : $dasref = $prot->get_all_DASFactories

=cut

sub get_all_DASFactories {
   my $self = shift;
   return [ $self->adaptor()->db()->_each_DASFeatureFactory ];
}

=head2 get_all_DASFeatures

  Arg [1]    : none
  Example    : $features = $prot->get_all_DASFeatures;
  Description: Retreives a hash reference to a hash of DAS feature
               sets, keyed by the DNS, NOTE the values of this hash
               are an anonymous array containing:
                (1) a pointer to an array of features;
                (2) a pointer to the DAS stylesheet
              TODO: Abstract to a DBLinkContainer obj
  Returntype : hashref of Bio::SeqFeatures
  Exceptions : ?
  Caller     : webcode

=cut

sub get_all_DASFeatures{
  my ($self,@args) = @_;
  $self->{_das_features} ||= {}; # Cache
  my %das_features;
  foreach my $dasfact( @{$self->get_all_DASFactories} ){
    my $dsn  = $dasfact->adaptor->dsn;
    my $name = $dasfact->adaptor->name;
    $name ||= $dasfact->adaptor->url .'/'. $dsn;
    if( $self->{_das_features}->{$name} ){ # Use cached
      $das_features{$name} = $self->{_das_features}->{$name};
      next;
    }
    else{ # Get fresh data
      my @featref = $dasfact->fetch_all_by_DBLink_Container( $self );
      $self->{_das_features}->{$name} = [@featref];
      $das_features{$name} = [@featref];
    }
  }
  return \%das_features;
}

1;










