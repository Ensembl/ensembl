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

Designed as an abstraction over the database specific DnaAlignFeatureAdaptors. The proxy gene adaptor normally behaves just as a normal core GeneAdaptor, 
however, for certain requests it may decide to instead forward
the request to another database (such as the EST database if it is available).

=head1 CONTACT

  Arne Stabenau: stabenau@ebi.ac.uk
  Ewan Birney  : birney@ebi.ac.uk
  Graham McVicker : mcvicker@ebi.ac.uk

=head1 APPENDIX

=cut

use strict;

package Bio::EnsEMBL::DBSQL::ProxyDnaAlignFeatureAdaptor;

use Bio::EnsEMBL::DBSQL::DnaAlignFeatureAdaptor;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::DBSQL::DnaAlignFeatureAdaptor);


#override super class constructor
sub new {
  my($class, $db, $core_adaptor) = @_;

  #call superclass constructor
  my $self = $class->SUPER::new($db);

  $self->{'_core_adaptor'} = $core_adaptor;
  return $self;
}

#override generic_fetch_method to call appropriate db
sub generic_fetch {
  my($self, $constraint, $logic_name) = @_;

  my $est_db = $self->db()->est_DBAdaptor();


  if(defined $est_db && $logic_name eq 'ex_e2g_feat') {
    #forward request to the EST database
    my $est_adaptor = $est_db->get_DnaAlignFeatureAdaptor();
    return $est_adaptor->generic_fetch($constraint, $logic_name);
  }
  
  #use the core adaptor
  return $self->{'_core_adaptor'}->generic_fetch($constraint, $logic_name);
}










