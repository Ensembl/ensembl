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
use Bio::EnsEMBL::DBSQL::SNPAdaptorI;

use vars '@ISA';

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor Bio::EnsEMBL::SNPAdaptorI);

#implement the SNPAdaptorI interface
use implements qw(Bio::EnsEMBL::DBSQL::SNPAdaptorI);


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

  
  my $sql = 
      "SELECT snp_chrom_start, chrom_strand, 
              refsnpid, anosnpid, tscid, hgbaseid, type, clone
       FROM snp
       WHERE chr_name = '"  . $slice->chr_name() . "'
       AND   snp_chrom_start > " . $slice->chr_start() . "
       AND   snp_chrom_start < " . $slice->chr_end();

  my $sth = $self->prepare($sql);
  
  eval {
    $sth->execute();
  };

  return () if $@;

  my @snps;
  
  while(my $arr = $sth->fetchrow_arrayref()) {
    my ($snp_start, $strand, $refsnp_id, 
	$tscid, $hgbaseid, $type, $acc) = @{$arr};



    # this snp has not been seen before, create a new snp object
    my $snp = Bio::EnsEMBL::External::Variation->new(
			  -start      => $snp_start - $slice->chr_start() + 1
			  -end        => $snp_start - $slice->chr_start() + 1
			  -strand     => $strand
                          -score      => 1,
			  -source_tag => 'dbSNP' );

    #Add db links to the snp variation object
    my $link = new Bio::Annotation::DBLink;
    $link->database('dbSNP');
    $link->primary_id($refsnp_id);
    $snp->add_DBLink($link);

    if ($hgbaseid) {
      my $link2 = new Bio::Annotation::DBLink;
      $link2->database('HGBASE');
      $link2->primary_id($hgbaseid);
      $snp->add_DBLink($link2);
    }
    if ($tscid) {
      my $link3 = new Bio::Annotation::DBLink;
      $link3->database('TSC-CSHL');
      $link3->primary_id($tscid);
      $snp->add_DBLink($link3);
    }

    push @snps, $snp;
  }
	
  return @snps;
}
    
    
    
1;

__END__
