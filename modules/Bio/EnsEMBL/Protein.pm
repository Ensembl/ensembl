
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

Protein.pm - DEPRECATED - Use Bio::EnsEMBL::Translation instead

=head1 SYNOPSIS

=head1 DESCRIPTION

Do not use this object.  Use Bio::EnsEMBL::Translation instead.

=head1 METHODS

=cut

package Bio::EnsEMBL::Protein;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Utils::Exception qw(throw deprecate);
use Bio::EnsEMBL::Translation;

@ISA = qw( Bio::EnsEMBL::Translation );


sub new {
  deprecate("The Protein class is deprecated. Use the Translation and " .
            "Transcript classes instead.\n");
  return bless {}, 'Bio::EnsEMBL::Protein';
}


=head2 species

  Description : DEPRECATED  no replacement exists.

=cut

sub species{
   my $obj = shift;

   deprecate("You should not have to use this method.");
   if( @_ ) {
      my $value = shift;
      $obj->{'species'} = $value;
    }
    return $obj->{'species'};

}

=head2 primary_seq

  Description: DEPRECATE do not use

=cut

sub primary_seq{
    my ($self) = @_;   
    return $self;
}



=head2 gene

  Description : DEPRECATED 
     use GeneAdaptor::fetch_by_translation_stable_id instead

=cut

sub gene {
  my $self = shift;

  deprecate("Use GeneAdaptor::fetch_by_translation_stable_id");

  my $adaptor = $self->adaptor();
  if(!$adaptor) {
    warning("Cannot get gene without attached adaptor");
    return undef;
  }

  my $stable_id = $self->stable_id();
  if(!$stable_id) {
    warning("Cannot get gene without translation stable_id");
    return undef;
  }
  my $g_adaptor = $adaptor->db->get_GeneAdaptor();
  return $g_adaptor->fetch_by_translation_stable_id($stable_id);
}


=head2 transcript

  Description : DEPRECATED
   use Bio::EnsEMBL::TranscriptAdaptor::fetch_by_translation_stable_id instead

=cut

sub transcript {
  my $self = shift;
  deprecate("Use TranscriptAdaptor::fetch_by_translation_stable_id");

  my $adaptor = $self->adaptor();
  if(!$adaptor) {
    warning("Cannot get transcript without attached adaptor");
    return undef;
  }

  my $stable_id = $self->stable_id();
  if(!$stable_id) {
    warning("Cannot get transcript without translation stable_id");
    return undef;
  }

  my $tr_adaptor = $adaptor->db->get_TranscriptAdaptor();
  return $tr_adaptor->fetch_by_translation_stable_id($stable_id);
}


=head2 molecule

  Description: DEPRECATED - No replacement.

=cut

sub molecule{
   my $obj = shift;
   deprecate("You should not have to use this method.");
   if( @_ ) {
      my $value = shift;
      $obj->{'molecule'} = $value;
    }
    return $obj->{'molecule'};

}

=head2 moltype
 
  Description: DEPRECATED - No replacement.

=cut

sub moltype{
   my $obj = shift;
   deprecate("You should not have to use this method.");
   if( @_ ) {
      my $value = shift;
      $obj->{'moltype'} = $value;
    }
    return $obj->{'moltype'};

}

=head2 top_SeqFeatures
    
  Description: DEPRECATED - use TranslationAdaptor::get_all_ProteinFeatures

=cut
    
sub top_SeqFeatures {
  my $self = shift;
  deprecate("Use TranslationAdaptor::get_all_ProteinFeatures");
  return @{$self->get_all_ProteinFeatures()};
}



=head2 get_all_Families

  Description: DEPRECATED

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

  Description: DEPRECATED

=cut

sub add_Family{
  my ($self,$value) = @_;
  
  unless(defined $value && 
	 $value->isa('Bio::EnsEMBL::ExternalData::Family::Family')) {
    throw("[$value] is not a Family object");
  }

  push(@{$self->{'_family'}},$value);  
}



=head2 get_all_ProfileFeatures

  Description : DEPRECATED - use Translation::get_all_ProteinFeatures('pfscan')

=cut

sub get_all_ProfileFeatures{
  my ($self) = @_;
  deprecate("Use Translation::get_all_ProteinFeatures('pfscan') instead");
  return $self->get_all_ProteinFeatures('pfscan');
}

=head2 add_Profile

  Description : DEPRECATED 

=cut

sub add_Profile{
  my ($self,$value) = @_;
  
  if ((!defined $value) || (!$value->isa('Bio::EnsEMBL::ProteinFeature'))) {
    throw("The Protein Feature added is not defined or is not a "
		 ."protein feature object");
  }

  push(@{$self->{'_profile'}},$value); 
  
}

=head2 get_all_blastpFeatures 

  Description : DEPRECATED use Translation::get_all_ProteinFeatures('blastp')

=cut

sub get_all_blastpFeatures{
  my $self = shift;
  deprecate("Use Translation::get_all_ProteinFeatures('blastp') instead");
  return $self->get_all_ProteinFeatures('blastp');
}

=head2 add_blastp

  Description: DEPRECATED

=cut

sub add_blastp{
  my ($self,$value) = @_;
        
  if ((!defined $value) || (!$value->isa('Bio::EnsEMBL::ProteinFeature'))) {
    throw("[$value] is not defined or is not a protein feature object");
  }

  push(@{$self->{'_blastp'}},$value); 
}

=head2 get_all_PrintsFeatures

  Description: DEPRECATED - use Translation::get_all_ProteinFeatures('prints');

=cut

sub get_all_PrintsFeatures{
  my ($self) = @_;
  deprecate('Use Translation::get_all_ProteinFeatures("prints") instead.');
  return $self->get_all_ProteinFeatures('prints');
}

=head2 add_Prints

  Description: DEPRECATED

=cut

sub add_Prints{
  my ($self,$value) = @_;
  if ((!defined $value) || (!$value->isa('Bio::EnsEMBL::ProteinFeature'))) {
    throw("[$value] is not a protein feature object");
  }
  
  push(@{$self->{'_prints'}},$value); 
}


=head2 get_all_PfamFeatures

  Description:  DEPRECATED use Translation::get_all_ProteinFeatures('pfam');

=cut

sub get_all_PfamFeatures{
  my ($self) = @_;
  deprecate('Use Translation::get_all_ProteinFeatures("pfam") instead');
  return $self->get_all_ProteinFeatures('pfam');
}


=head2 add_Pfam

  Description:  DEPRECATED

=cut

sub add_Pfam{
 my ($self,$value) = @_;
        
 if ((!defined $value) || (!$value->isa('Bio::EnsEMBL::ProteinFeature'))) {
   throw("[$value] is not a protein feature object");
 }

 push(@{$self->{'_pfam'}},$value); 
}


=head2 get_all_PrositeFeatures

  Description: DEPRECATED use Translation::get_all_PrositeFeatures('scanprosite');

=cut

sub get_all_PrositeFeatures{
  my ($self) = @_;
  deprecate('Use Translation::get_all_ProteinFeatures("scanprosite")');
  return $self->get_all_ProteinFeatures("scanprosite");
}

=head2 add_Prosite

  Description:  DEPRECATED

=cut

sub add_Prosite{
  my ($self,$value) = @_;
    
  if ((!defined $value) || (!$value->isa('Bio::EnsEMBL::ProteinFeature'))) {
    throw("[$value] is not a protein feature object");
  }

  push(@{$self->{'_prosite'}},$value); 
}

=head2 get_all_SigpFeatures

  Description: DEPRECATED - use Translation::get_all_ProteinFeatures instead

=cut

sub get_all_SigpFeatures{
  my ($self) = @_;
  deprecate("Use Translation::get_all_ProteinFeatures('signalp') instead.");
  return $self->get_all_ProteinFeatures('signalp');

}

=head2 add_Sigp

  Description: DEPRECATED

=cut

sub add_Sigp{
  my ($self,$value) = @_;
  
  if ((!defined $value) || (!$value->isa('Bio::EnsEMBL::ProteinFeature'))) {
    throw("[$value] is not defined or is not a protein feature object");
  }
  push(@{$self->{'_sigp'}},$value); 
}


=head2 get_all_TransmembraneFeatures

  Description: DEPRECATED - use Translation::get_all_ProteinFeatures

=cut

sub get_all_TransmembraneFeatures{
  my $self = shift;
  deprecate("Use Translation::get_all_ProteinFeatures('tmhmm') instead");
  return $self->get_all_ProteinFeatures('tmhmm');
}

=head2 add_Transmembrane

  Description: DEPRECATED

=cut

sub add_Transmembrane{
  my ($self,$value) = @_;
  
  if ((!defined $value) || (!$value->isa('Bio::EnsEMBL::ProteinFeature'))) {
    throw("The Protein Feature added is not defined" .
		 "or is not a protein feature object");
    }
  push(@{$self->{'_transmembrane'}},$value); 
}

=head2 get_all_CoilsFeatures

  Description: DEPRECATED use Translation::get_all_ProteinFeatures instead

=cut

sub get_all_CoilsFeatures{
  my $self = shift;
  deprecate("Use Translation::get_all_ProteinFeatures('ncoils') instead.");
  return $self->get_all_ProteinFeatures('ncoils');
}

=head2 add_Coils

  Description: DEPRECATED

=cut

sub add_Coils{
  my ($self,$value) = @_;
  
  if ((!defined $value) || (!$value->isa('Bio::EnsEMBL::ProteinFeature'))) {
    throw("[$value] is not a protein feature object");
    }
  push(@{$self->{'_coils'}},$value); 
}


=head2 get_all_LowcomplFeatures

  Description: DEPRECATED use Translation::get_all_ProteinFeatures('seg');

=cut

sub get_all_LowcomplFeatures{
  my $self = shift;
  deprecate("Use Translation::get_all_ProteinFeatures('seg')");
  return $self->get_all_ProteinFeatures('seg');
}

=head2 add_LowCompl

  Description: DEPRECATED

=cut

sub add_Lowcompl{
  my ($self,$value) = @_;

  if ((!defined $value) || (!$value->isa('Bio::EnsEMBL::ProteinFeature'))) {
    throw("[$value] is not a protein feature object");
  }

  push(@{$self->{'_lowcompl'}},$value); 
}

=head2 get_all_SuperfamilyFeatures

Description: DEPRECATED - use Translation::get_all_ProteinFeatures instead

=cut

sub get_all_SuperfamilyFeatures{
  my $self = shift;
  deprecate("Use Translation::get_all_ProteinFeatures('superfamily') instead");
  return $self->get_all_ProteinFeatures('superfamily');
}

=head2 add_Superfamily

  Description: DEPRECATED

=cut

sub add_Superfamily{
  my ($self,$value) = @_;
  
  if ((!defined $value) || (!$value->isa('Bio::EnsEMBL::ProteinFeature'))) {
    throw("[$value] is not a protein feature object");
  } 

  push(@{$self->{'_superfamily'}},$value); 
}


=head2 length

  DEPRECATED : use $transcript->translate->length instead

=cut

sub length{
  my $self = shift;
  deprecate("use transcript->translate->length instead");
  return $self->transcript->translate->length();
}


sub primary_id {
  my $self = shift;
  deprecate("Use Translation::dbID or Translation::stable_id instead.");
  return $self->stable_id || $self->dbID;
}


sub id {
  my $self = shift;
  deprecate("Use Translation::dbID or Translation::stable_id instead.");
  return $self->stable_id || $self->dbID;
}


sub seq {
  my $self = shift;
  deprecate("Use Translation::translate instead.");
  return $self->transcript->translate()->seq();
}


1;










