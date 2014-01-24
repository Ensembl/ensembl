=head1 LICENSE

Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::DBSQL::ProxySNPAdaptor

=head1 SYNOPSIS

Designed as an abstraction over the database specific SNPAdaptors. This
is written right now to serve as a replacement for a core SNPadaptor
which doesn''t even exist yet and probably never will since SNPs are
taken from external databases. In the future some sort of DBRegistry may
remove the need for this class.

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::ProxySNPAdaptor;

use strict;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw);

use vars ('@ISA', '$AUTOLOAD');

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

=head2 fetch_attributes_only

  Arg [1]    : int refsnp_id
  Arg [2]    : (optional) string source
  Example    : none
  Description: Retrieves a snp objcet from the SNP database but does not
               populate the location information.  This is necessary given 
               the current state of the snp database because location 
               information has to be retrieved differently for different 
               species!
  Returntype : Bio::EnsEMBL::SNP
  Exceptions : none
  Caller     : snpview

=cut



sub fetch_attributes_only{
    my ( $self, @args ) = @_;

  my $lite_db = Bio::EnsEMBL::Registry->get_db($self->db(),'lite');
  my $snp_db = Bio::EnsEMBL::Registry->get_db($self->db(),'SNP');

  if( defined $snp_db ) {
      my $snp_adaptor = $snp_db->get_SNPAdaptor();
      return $snp_adaptor->fetch_attributes_only( @args );
  }

  if( defined $lite_db ) {
      my $snp_adaptor = $lite_db->get_SNPAdaptor();
      return $snp_adaptor->fetch_attributes_only( @args );
  }

}




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

  if( defined $lite_db ) {
      my $snp_adaptor = $lite_db->get_SNPAdaptor();
      if($snp_adaptor->can($method)) {
	  return $snp_adaptor->$method(@args);
      } 
  }

  if( defined $snp_db ) {
      my $snp_adaptor = $snp_db->get_SNPAdaptor();
      if($snp_adaptor->can($method)) {
	  return $snp_adaptor->$method(@args);
      }
  }



  throw("The requested method $method could not be found in lite or snp" );
}

sub DESTROY {
}

1;

__END__
