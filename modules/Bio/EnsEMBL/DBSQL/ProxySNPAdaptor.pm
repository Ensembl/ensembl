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

It would be nicer if this class and the other SNPAdaptors implemented
a common interface. 

=head1 CONTACT

  Arne Stabenau: stabenau@ebi.ac.uk
  Ewan Birney  : birney@ebi.ac.uk
  Graham McVicker : mcvicker@ebi.ac.uk

=head1 APPENDIX

=cut

use strict;

package Bio::EnsEMBL::DBSQL::ProxySNPAdaptor;

use Bio::EnsEMBL::DBSQL::SNPAdaptorI;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;

use vars '@ISA';

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor Bio::EnsEMBL::SNPAdaptorI);

#implement the SNPAdaptorI interface 
use implements qw(Bio::EnsEMBL::DBSQL::SNPAdaptorI);



sub fetch_by_Slice {
  my ($self, @args) = @_;

  my $lite_db = $self->db()->lite_DBAdaptor();
  my $snp_db = $self->db()->SNP_DBAdaptor();

  if(defined $lite_db) {
    #use the Lite database if it is available
    return $lite_db->get_SNPAdaptor()->fetch_by_Slice(@args);
  } elsif(defined $snp_db) {
    #use the snp database if it is available
    return $snp_db->get_SNPAdaptor()->fetch_by_Slice(@args);
  }

  #There is no core SNPAdaptor so throw an exception if lite and SNP
  #databases are unavailable
  $self->throw("Lite database unavailable. Unable to create SNP adaptor");

  return undef;
}

sub fetch_by_SNP_id {
  my ($self, @args) = @_;

  my $snp_db = $self->db()->SNP_DBAdaptor();

  if(defined $snp_db) {
    return $snp_db->get_SNPAdaptor()->fetch_by_SNP_id(@args);
  }

  $self->throw("SNP database unavailable. Unable to create SNP adaptor");
}



1;

__END__
