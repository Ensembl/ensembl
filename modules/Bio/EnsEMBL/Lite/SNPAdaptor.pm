# EnsEMBL Gene reading writing adaptor for mySQL
#
# Copyright EMBL-EBI 2001
#
# Author: Arne Stabenau
# 
# Date : 15.07.2002
#

=head1 NAME

Bio::EnsEMBL::DBSQL::GeneLiteAdaptor - 
MySQL Database queries to retrieve genes quickly from denormalized tables.

=head1 SYNOPSIS

=head1 CONTACT

  Arne Stabenau: stabenau@ebi.ac.uk
  Ewan Birney  : birney@ebi.ac.uk

=head1 APPENDIX

=cut

use strict;

package Bio::EnsEMBL::Lite::SNPAdaptor;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::SNP;

use vars '@ISA';

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);


=head2 fetch_all_by_Slice

  Arg  1    : Bio::EnsEMBL::Slice $slice
              The slice we want SNPs on
  Function  : retrieve all the SNPs on this slice. 
              uses Lite databases transcript to get info
  Returntype: listreference of Bio::EnsEMBL::ExternalData::Variation
  Exceptions: none
  Caller    : Bio::EnsEMBL::Slice

=cut

sub fetch_all_by_Slice {
  my ($self, $slice ) = @_;

  my $slice_start = $slice->chr_start();
  my $slice_end   = $slice->chr_end();
  
  my %SNPS = qw( 12 dbSNP 13 WI 14 HGBASE 15 TSC-CSHL 16 ANO );
  my $QUERY = "select internal_id, chr_start, chr_end, chr_strand, type, range_type,
		      validated, alleles, snpclass, mapweight, ambiguity, source,
		      id_refsnp, id_wi, id_hgbase, id_tsc, id_ano 
                 FROM snp
                WHERE chr_name = ? AND chr_start >= ? and chr_start <= ? AND chr_end >= ?";

  my $sth = $self->prepare( $QUERY ); 
  eval {
  $sth->execute($slice->chr_name(), $slice_start - 500 , $slice_end, $slice_start);
  };
  return [] if $@;
  my @snps = ();  

  my %link_hash;
  my $link;

  while(my $arrayref = $sth->fetchrow_arrayref()) {
    my @links = ();
    foreach( sort keys %SNPS ) {
       my $V = $arrayref->[ $_ ];
       if( $V && $V ne '' ) {
         unless($link = $link_hash{"$SNPS{$_}:$V"}) {
           $link_hash{"$SNPS{$_}:$V"} = $link = Bio::EnsEMBL::DBEntry->new_fast( {'_dbname'     => $SNPS{$_}, '_primary_id' => $V });
         }
         push @links, $link;
       }
    }

    #create a snp object through a fast (hacky) constructor
    my $snp = Bio::EnsEMBL::SNP->new_fast(
		  { 'dbID'       => $arrayref->[0], 
		   '_gsf_start'  => $arrayref->[1] - $slice_start + 1,
		    '_gsf_end'    => $arrayref->[2] - $slice_start + 1,
		    '_snp_strand' => $arrayref->[3],
		    '_gsf_score'  => 1,
		    '_type'       => $arrayref->[4],
                    '_range_type' => $arrayref->[5],
                    '_validated'  => $arrayref->[6],
                    'alleles'    => $arrayref->[7],
                    '_ambiguity_code' => $arrayref->[10],
                    '_snpclass'   => $arrayref->[8],
                    '_mapweight'  => $arrayref->[9],
		    '_source_tag' => $arrayref->[11],
		    'link'        => \@links });
    push @snps, $snp;
  }
	
  return \@snps;
}

1;

