# EnsEMBL Gene reading writing adaptor for mySQL
#
# Copyright EMBL-EBI 2002
#
# Author: Graham McVicker
# 
# Date : 05.08.2002
#

=head1 NAME

Bio::EnsEMBL::DBSQL::ProxySNPAdaptor

=head1 SYNOPSIS

Designed as an abstraction over the database specific SNPAdaptors. This is 
written right now to serve as a replacement for a core SNPadaptor which 
doesn''t even exist yet and probably never will since SNPs are taken from
external databases. In the future some sort of DBRegistry may remove the
need for this class.  

=head1 CONTACT

  Arne Stabenau: stabenau@ebi.ac.uk
  Ewan Birney  : birney@ebi.ac.uk
  Graham McVicker : mcvicker@ebi.ac.uk

=head1 APPENDIX

=cut

use strict;

package Bio::EnsEMBL::DBSQL::ProxySNPAdaptor;

use Bio::EnsEMBL::DBSQL::ProxyAdaptor;

use vars '@ISA';

@ISA = qw(Bio::EnsEMBL::DBSQL::ProxyAdaptor);


=head2 fetch_by_Slice

  Arg [1]    : list of arbitrary args @args
  Example    : none
  Description: Forwards request to the Lite database if it is available
               If the lite database is not available the request is 
               forwarded to the SNP database. 
  Returntype : listref of Bio::EnsEMBL::ExternalData::Variation
  Exceptions : thrown if neither the SNP nor Lite databases are available
  Caller     : snpview

=cut

sub fetch_all_by_Slice {
  my ($self, @args) = @_;

  my $lite_db = $self->db()->get_db_adaptor('lite');
  my $snp_db = $self->db()->get_db_adaptor('SNP');

  if(defined $lite_db) {
    #use the Lite database if it is available
    return $lite_db->get_SNPAdaptor()->fetch_all_by_Slice(@args);
  } elsif(defined $snp_db) {
    #use the snp database if it is available
    return $snp_db->get_SNPAdaptor()->fetch_all_by_Slice(@args);
  }

  #There is no core SNPAdaptor so throw an exception if lite and SNP
  #databases are unavailable
  $self->throw("Lite and SNP databases are unavailable. " .
	       "Unable to create SNP adaptor");

  return undef;
}


=head2 fetch_by_SNP_id

  Arg [1]    : list of arbitrary args @args
  Example    : none
  Description: Forwards request to the Lite database if it is available
               If the lite database is not available the request is 
               forwarded to the SNP database. 
  Returntype : Bio::EnsEMBL::ExternalData::Variation
  Exceptions : thrown if neither the SNP nor Lite databases are available
  Caller     : snpview

=cut

sub fetch_by_SNP_id {
  my ($self, @args) = @_;

  my $snp_db = $self->db()->get_db_adaptor('SNP');

  if(defined $snp_db) {
    return $snp_db->get_SNPAdaptor()->fetch_by_SNP_id(@args);
  }

  $self->throw("SNP database unavailable. Unable to create SNP adaptor");
}



1;

__END__
