# EnsEMBL Gene reading writing adaptor for mySQL
#
# Copyright EMBL-EBI 2001
#
# Author: Graham McVicker
# 
# Date : 22.07.2002
#

=head1 NAME

Bio::EnsEMBL::DBSQL::ProxyDnaAlignFeatureAdaptor

=head1 SYNOPSIS

Designed as an abstraction over the database specific DnaAlignFeatureAdaptors. The proxy gene adaptor normally behaves just as a normal core 
DnaAlignFeatureAdaptor, however, for certain requests it may decide to 
instead forward the request to another database (such as the EST database 
if it is available).

=head1 CONTACT

  Graham McVicker : mcvicker@ebi.ac.uk
  Arne Stabenau: stabenau@ebi.ac.uk
  Ewan Birney  : birney@ebi.ac.uk

=head1 APPENDIX

=cut

use strict;

use Bio::EnsEMBL::DBSQL::ProxyAdaptor;
use Bio::EnsEMBL::DBSQL::DnaAlignFeatureAdaptor;

package Bio::EnsEMBL::DBSQL::ProxyDnaAlignFeatureAdaptor;

use vars qw(@ISA);

#for most methods we want to go through the DnaAlignFeatureAdaptor 
#implementation. We want to make a decision at the lower level however,
#and since most calls to the adaptor chain through other methods without
#going back through the proxy, we have to extend the DnaAlignFeatureAdaptor
#for this proxy to correctly function.  We don't want the proxy decision
#to be made too early anotherwards...
@ISA = qw(Bio::EnsEMBL::DBSQL::ProxyAdaptor 
	  Bio::EnsEMBL::DBSQL::DnaAlignFeatureAdaptor);

#override generic_fetch_method to call appropriate db

=head2 generic_fetch

  Arg [1]    : string $constraint
  Arg [2]    : string $logic_name
  Example    : none
  Description: Overrides the Bio::EnsEMBL::DNAAlignFeatureAdaptor generic_fetch
               method.  Normally requests will still be forwarded to the 
               core (primary) database with the exception of requests that
               have a logic_name 'ex_e2g_feat' which will use the external
               EST database if it is available.
  Returntype : Bio::EnsEMBL::DBSQL::DnaAlignFeatureAdaptor
  Exceptions : none
  Caller     : general

=cut

sub generic_fetch {
  my($self, $constraint, $logic_name, $mapper, $slice) = @_;

  if((defined $logic_name)&&($logic_name eq 'ex_e2g_feat')) {
    my $est_db = $self->db()->get_db_adaptor('est');
    
    if(defined $est_db) {
      #forward request to the EST database
      my $est_adaptor = $est_db->get_DnaAlignFeatureAdaptor();
      return $est_adaptor->generic_fetch($constraint, "", $mapper, $slice);
    }
  }
  
  #use the core adaptor
  return $self->{'_primary_adaptor'}->generic_fetch($constraint, $logic_name, $mapper, $slice);
}


=head2 store

  Arg [1]    : Bio::EnsEMBL::DnaDnaAlignFeature $feature
  Example    : none
  Description: overrides the store method to ensure that the est database
               is used to store ex_e2g_features
  Returntype : see DnaAlignFeatureAdaptor::store
  Exceptions : none
  Caller     : general

=cut

sub store {
  my ($self, @sfs) = @_;

  my $est_features = [];
  my $core_features = [];

  my $est_db = $self->db()->get_db_adaptor('est');
  if(defined $est_db) {
    foreach my $f (@sfs) {
      if($f->analysis() && $f->analysis()->logic_name() eq 'ex_e2g_feat') {
	push @$est_features, $f;
      } else {
	push @$core_features, $f;
      }
      #forward request to the EST db
      my $est_adaptor = $est_db->get_DnaAlignFeatureAdaptor();
      if(scalar @$est_features) {
	$est_adaptor->store(@$est_features);
      }
      if(scalar @$core_features) {
	$self->{'_primary_adaptor'}->store(@$core_features);
      }
    }
  }

  #use the core adaptor
  return $self->{'_primary_adaptor'}->store(@sfs);
}


=head2 remove

  Arg [1]    : Bio::EnsEMBL::DnaDnaAlignFeature $feature
  Example    : none
  Description: overrides the remove method to ensure that the est database
               is used to remove ex_e2g_features
  Returntype : none
  Exceptions : none
  Caller     : general

=cut

sub remove {
  my ($self, $feature) = @_;

  if($feature && $feature->analysis && 
     $feature->analysis()->logic_name() eq 'ex_e2g_feat') {
    my $est_db = $self->db()->get_db_adaptor('est');

    if(defined $est_db) {
      #forward request to the EST db
      my $est_adaptor = $est_db->get_DnaAlignFeatureAdaptor();
      return $est_adaptor->remove($feature);
    }
  }

  #use the core adaptor
  return $self->{'_primary_adaptor'}->remove($feature);
}

1;
