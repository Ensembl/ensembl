# EnsEMBL Gene reading writing adaptor for mySQL
#
# Copyright EMBL-EBI 2001
#
# Author: Arne Stabenau
# 
# Date : 15.07.2002
#

=head1 NAME

Bio::EnsEMBL::DBSQL::GeneLiteAdaptor - MySQL Database queries to retrieve genes quickly from denormalized tables.

=head1 SYNOPSIS

=head1 CONTACT

  Arne Stabenau: stabenau@ebi.ac.uk
  Ewan Birney  : birney@ebi.ac.uk

=head1 APPENDIX

=cut

use strict;

package Bio::EnsEMBL::Lite::SNPAdaptor;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::SNP;

use vars '@ISA';

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);


=head2 fetch_by_Slice

  Arg  1    : Bio::EnsEMBL::Slice $slice
              The slice we want SNPs on
  Function  : retrieve all the SNPs on this slice. 
              uses Lite databases transcript to get info
  Returntype: list of Bio::EnsEMBL::ExternalData::Variation
  Exceptions: none
  Caller    : Bio::EnsEMBL::Slice

=cut

sub fetch_by_Slice {
  my ($self, $slice ) = @_;

  my $slice_start = $slice->chr_start();
  my $slice_end   = $slice->chr_end();
  
  my $sth = $self->prepare(
      "SELECT chr_start, chr_strand, 
              refsnpid, tscid, hgbaseid, type
       FROM   snp
       WHERE  chr_name = ? AND chr_start >= ? AND chr_start <= ?");

  $sth->execute($slice->chr_name(), $slice_start, $slice_end);
  
  my($chr_start, $chr_strand, $refsnpid, $tscid, $hgbaseid, $type); 
  
  $sth->bind_columns(\$chr_start, \$chr_strand, \$refsnpid, 
		     \$tscid, \$hgbaseid, \$type);

  my @snps;  

  my %link_hash;
  my $link;

  while($sth->fetch()) {
    # this snp has not been seen before, create a new snp object
    
    my $snp = Bio::EnsEMBL::SNP->new_fast(
		  { '_gsf_start'  => $chr_start - $slice_start + 1,
		    '_gsf_end'    => $chr_start - $slice_start + 1,
		    '_snp_strand' => $chr_strand,
		    '_gsf_score'  => 1,
		    '_type'       => $type });
    

    #Add db links to the snp variation object
    unless($link = $link_hash{"dbSNP:$refsnpid"}) {
      $link = new Bio::Annotation::DBLink;
      $link->database('dbSNP');
      $link->primary_id($refsnpid);
      $link_hash{"dbSNP:$refsnpid"} = $link;
    }
    $snp->add_DBLink($link);
   
    if ($hgbaseid) {
      unless($link = $link_hash{"HGBASE:$hgbaseid"}) {
	$link = new Bio::Annotation::DBLink;
	$link->database('HGBASE');
	$link->primary_id($hgbaseid);
	$link_hash{"HGBASE:$hgbaseid"} = $link;
      }
      $snp->add_DBLink($link);
    }
    if ($tscid) {
      unless($link = $link_hash{"TSC-CSHL:$tscid"}) {
	$link = new Bio::Annotation::DBLink;
	$link->database('TSC-CSHL');
	$link->primary_id($tscid);
	$link_hash{"TSC-CSHL:$tscid"} = $link;
      }
      $snp->add_DBLink($link);
    }

    push @snps, $snp;
  }
	
  return @snps;
}
    
    
    
1;

__END__
