# EnsEMBL Gene reading writing adaptor for mySQL
#
# Copyright EMBL-EBI 2002
#
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

  Post questions to the Ensembl development list <ensembl-dev@ebi.ac.uk>

=head1 METHODS

=cut

use strict;

package Bio::EnsEMBL::DBSQL::ProxySNPAdaptor;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;


use vars ('@ISA', '$AUTOLOAD');

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);



=head2 AUTOLOAD

  Arg [1]    : list of arbitrary values @args
               a list of arguments to pass to the request method
  Example    : -
  Description: AUTOLOAD method should not be called directly.  It is called
               implicitly when a method requested from this class cannot be
               found. This method first tries to execute the requested method
               in the primary adaptor.  If the method cannot be found then 
               it searches the other attached databases for equivalent adaptors
               and tries then one at a time.
  Returntype : arbitrary
  Exceptions : thrown if the requested method cannot be found on the primary
               adaptor or on any of the attached databases.
  Caller     : called implicitly by perl

=cut

sub AUTOLOAD {
  my ($self, @args) =  @_;
  
  #determine the method which was called
  my $method = $AUTOLOAD;
  
  #strip out fully qualified method name
  $method =~ s/.*:://;


  my $lite_db = Bio::EnsEMBL::Registry->get_db($self->db(),'lite');
  my $snp_db = Bio::EnsEMBL::Registry->get_db($self->db(),'SNP');
  
  #
  # First try the primary adaptor
  #
  if( defined $lite_db ) {
      my $snp_adaptor = $lite_db->get_SNPAdaptor();
      if($snp_adaptor->can($method)) {
	  return $snp_adaptor->$method(@args);
      } 
  }
  
  my $snp_adaptor = $snp_db->get_SNPAdaptor();
  if($snp_adaptor->can($method)) {
      return $snp_adaptor->$method(@args);
  }


  throw("The requested method $method could not be found in lite or snp" );
}


1;

__END__
