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

package Bio::EnsEMBL::DBSQL::ProxyDnaAlignFeatureAdaptor;

use vars qw(@ISA);

#for most methods we want to go through the DnaAlignFeatureAdaptor 
#implementation. We want to make a decision at the lower level however,
#and since most calls to the adaptor chain through other methods without
#going back through the proxy, we have to extend the DnaAlignFeatureAdaptor
#for this proxy to correctly function.  We don't want the proxy decision
#to be made too early anotherwards...
@ISA = qw(Bio::EnsEMBL::DBSQL::ProxyAdaptor 
	  Bio::EnsEMBL::DBSQL::DNAAlignFeatureAdaptor);

#override generic_fetch_method to call appropriate db
sub generic_fetch {
  my($self, $constraint, $logic_name) = @_;



  if($logic_name eq 'ex_e2g_feat') {
    my $est_db = $self->db()->get_db_adaptor('est');
    
    if(defined $est_db) {
      #forward request to the EST database
      my $est_adaptor = $est_db->get_DnaAlignFeatureAdaptor();
      return $est_adaptor->generic_fetch($constraint, $logic_name);
    }
  }
    
  #use the core adaptor
  return $self->{'_primary_adaptor'}->generic_fetch($constraint, $logic_name);
}










