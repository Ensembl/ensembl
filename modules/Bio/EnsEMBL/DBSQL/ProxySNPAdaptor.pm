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

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

#implement the SNPAdaptorI interface 
use implements qw(Bio::EnsEMBL::DBSQL::SNPAdaptorI);



sub fetch_by_Slice {
  my ($self, @args) = @_;

  my $lite_db = $self->db()->lite_DBAdaptor();
 
  if(defined $lite_db) {
    #use the Lite database if it is available
    return $lite_db->get_SNPAdaptor()->fetch_by_Slice(@args);
  }

  #currently the only type of SNP adaptor is the lite one.  If lite is
  #unavailable throw an exception
  $self->throw("Lite database unavailable. Unable to create SNP adaptor");

  return undef;
}



1;

__END__
