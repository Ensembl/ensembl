
#
# BioPerl module for Bio::EnsEMBL::FeatureFactory
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::FeatureFactory - A static class for providing features

=head1 SYNOPSIS

  # get features,featurepairs and analysis objects
  # These could be very different implementations, eg
  # C structs from the ensembl-cext module

  $feature = Bio::EnsEMBL::FeatureFactory->new_feature();
  $fp      = Bio::EnsEMBL::FeatureFactory->new_feature_pair();
  $ana     = Bio::EnsEMBL::FeatureFactory->new_analysis();

  # fill in feature attributes normally
  # (you cannot set them at the new function stage)
  $feature->start(10);

  # feature pair has a one-shot setting mechanism

  $fp->set_all_fields($start,$end,$strand,$score,$source,$primary,$seqname,$hstart,$hend,$hstrand,$hscore,$hsource,$hprimary,$hseqname);

  # you have to use specialised analysis objects from the 
  # feature factory as they have to come from the same implementation
  # of the feature factory

  $f->analysis($ana);


=head1 DESCRIPTION

FeatureFactory is a way of abstracting out the ability to get 
features, featurepairs (two features linked together in a database search
result) and analysis without needing to know the implementation.

The idea is that in this module we take full advantage of Perl's run
time object loading to figure out the best implementation to use
in this situation.

=head1 AUTHOR - Ewan Birney

Email birney@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...



package Bio::EnsEMBL::FeatureFactory;
use vars qw($USE_PERL_ONLY $ENSEMBL_EXT_LOADED);
use strict;

use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Analysis;


$USE_PERL_ONLY = 0;

BEGIN {
    
    eval { 
	require EnsemblExt;
    };
    if( $@ ) {
	$ENSEMBL_EXT_LOADED = 0;
    } else {
	$ENSEMBL_EXT_LOADED = 1;
    }
}


=head2 new_feature

 Title   : new_feature
 Usage   : $f = Bio::EnsEMBL::FeatureFactory->new_feature()
 Function: Provides a new feature for ensembl but without any
           knowledge of its implementation
 Returns : A Bio::EnsEMBL::SeqFeatureI compliant object
 Args    : none (you have to fill in the fields after construction)


=cut

sub new_feature{
   my ($self,@args) = @_;

  if( $ENSEMBL_EXT_LOADED == 1 && $USE_PERL_ONLY == 0 ) {
      # catch for @args being passed in.
      $self = Bio::EnsEMBL::Ext::SeqFeature->new();
      return $self;
  }

   return Bio::EnsEMBL::SeqFeature->new();

}

=head2 new_feature_pair

 Title   : new_feature_pair
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub new_feature_pair{
    my $self;
    if( $ENSEMBL_EXT_LOADED == 1 && $USE_PERL_ONLY == 1 ) {
	# catch for @args being passed in.
	$self = Bio::EnsEMBL::Ext::FeaturePair->new();
	return $self;
    }

    $self = Bio::EnsEMBL::FeaturePair->new();
    my $one = Bio::EnsEMBL::SeqFeature->new();
    my $two = Bio::EnsEMBL::SeqFeature->new();

    $self->feature1($one);
    $self->feature2($two);

    return $self;

}

=head2 new_analysis

 Title   : new_analysis
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub new_analysis{

    if( $ENSEMBL_EXT_LOADED == 1 && $USE_PERL_ONLY == 1 ) {
	# catch for @args being passed in.
	my $self = Bio::EnsEMBL::Ext::Analysis->new();
	return $self;
    }

    return Bio::EnsEMBL::Analysis->new();
}

1;


