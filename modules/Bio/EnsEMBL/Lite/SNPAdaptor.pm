=head1 LICENSE

  Copyright (c) 1999-2009 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::DBSQL::GeneLiteAdaptor - 
MySQL Database queries to retrieve genes quickly from denormalized tables.

=head1 SYNOPSIS

=head1 METHODS

=cut


package Bio::EnsEMBL::Lite::SNPAdaptor;

use strict;
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

  my @snps;

  my $from_strand = $slice->strand();

  # warn "LITE snps....";
  # wherever this slice is, it needs to be converted to
  # a toplevel slices since all snps in the lite database 
  # are stored on toplevel seqregions

  my @projection = @{$slice->project('toplevel')};

  my $link_col_idx = 12;
  my @link_dbs = ('dbSNP', 'WI', 'HGBASE', 'TSC-CSHL', 'ANO');

  my %link_cache = ();

  foreach my $segment (@projection) {
    my $from_start = $segment->from_start();
    my $from_end   = $segment->from_end();
    my $top_slice  = $segment->to_Slice();
    my $top_slice_start  = $top_slice->start();
    my $top_slice_end    = $top_slice->end();
    my $top_slice_strand = $top_slice->strand();

    my $sth = $self->prepare
      ("SELECT internal_id, chr_start, chr_end, chr_strand, type, " .
       "       range_type, validated, alleles, snpclass, mapweight, ".
       "       ambiguity, source, id_refsnp, id_wi, id_hgbase, id_tsc, " .
       "       id_ano " .
       "FROM   snp " .
       "WHERE  chr_name = ? " .
       "AND    chr_start >= ? " .
       "AND    chr_start <= ? " .
       "AND    chr_end >= ?");

    $sth->execute($top_slice->seq_region_name(),
                  $top_slice_start - 500,
                  $top_slice_end,
                  $top_slice_start);

    while(my $arrayref = $sth->fetchrow_arrayref()) {
      my @links = ();

      # loop over the last columns of the row to retrieve a single external
      # database link for each column

      for(my $i = 0; $i < scalar(@link_dbs); $i++) {
        my $link_id = $arrayref->[$link_col_idx + $i];

        next if(!$link_id);

        my $link_db = $link_dbs[$i];
        my $link = $link_cache{"$link_db:$link_id"};

        if(!$link) {
          $link = Bio::EnsEMBL::DBEntry->new_fast
            ({'dbname' => $link_db,
              'primary_id'   => $link_id,
              'display_name' => $link_id});
          $link_cache{"$link_db:$link_id"} = $link;
        }

        push @links, $link;
      }

      #create a snp object through a fast (hacky) constructor
      my $status = $arrayref->[6];
      $status =~ s/-/ /;
      if($status && $status ne 'no info') {
        $status = "proven $status";
      } else {
        $status = 'suspected';
      }

      # coordinates must be adjusted so that they are first 
      # relative to the start of the top level slice (rather than absolute)
      # and then adjusted so they are relative to the start of the
      # original requested slice (w/ from_start)
      my($start,$end,$strand);
      if($top_slice_strand == 1) {
        $start = $arrayref->[1] - $top_slice_start + $from_start;
        $end   = $arrayref->[2] - $top_slice_start + $from_start;
        $strand = $arrayref->[3];
      } else {
        $start = $top_slice_end - $arrayref->[2] + $from_start;
        $end   = $top_slice_end - $arrayref->[1] + $from_start;
        $strand = $arrayref->[3] * -1;
      }

      push @snps, Bio::EnsEMBL::SNP->new_fast({
        'dbID'       => $arrayref->[0],
        '_gsf_start'  => $start,
        '_gsf_end'    => $end,
        '_snp_strand' => $strand,
        '_gsf_score'  => 1,
        '_type'       => $arrayref->[4],
        '_range_type' => $arrayref->[5],
        '_validated'  => $arrayref->[6],
        'status'     => $status,
        'alleles'    => $arrayref->[7],
        '_ambiguity_code' => $arrayref->[10],
        '_snpclass'   => $arrayref->[8],
        '_mapweight'  => $arrayref->[9],
        '_source' => $arrayref->[11],
        '_source_tag' => $arrayref->[11],
        'link'        => \@links,
        '_unique_id'  => "$arrayref->[0]:$arrayref->[1]"
      });
    }
  }
	
  return \@snps;
}

sub fetch_all_by_Slice_transcript_ids {
  my $self = shift;
  my $slice = shift;
  my $transcript_ids = shift;
  my $DB             = shift || 'core';
  my $snps = $self->fetch_all_by_Slice( $slice );
  my %SNPS = ();
  foreach my $transid ( @{$transcript_ids||[]} ) {
 # warn "TRANSCRIPT: $transid";
    my $sth = $self->prepare(qq(select gs.snp_id, gs.type, gs.aminoacid_start,
                gs.aminoacid_offset, gs.wildtype_aminoacid,
                gs.aminoacids, s.internal_id, s.chr_start
           from gene_snp as gs, snp as s
          where gs.transcript_id = ? and gs.db = "$DB" and
                gs.snp_id = s.snp_id
    ));
    $sth->execute( $transid );
    while(my $a = $sth->fetchrow_arrayref()) {
      $SNPS{"$a->[6]:$a->[7]"}{$transid} = [ @$a ];
    }
  }
  foreach my $snp ( @$snps ) {
    $snp->{'_transcripts'} = {};
    my $snptype = '99:';
    if( $SNPS{$snp->{'_unique_id'}} ) {
      #warn ">>> $snp->{'_unique_id'}";
      foreach my $transid ( keys %{$SNPS{$snp->{'_unique_id'}}} ) {
        my $a = $SNPS{$snp->{'_unique_id'}}{$transid};
        $snp->{'_transcripts'}{$transid} = $a;
        $snptype = $a->[1] if $a->[1] lt $snptype;
      }
    }
    $snp->{'_local_type'} = $snptype;
  }
  return $snps;
}

sub fetch_attributes_only_lite{
  my $self = shift;

  my $refsnp_id = shift;
  my $source = shift || 'dbSNP';

  my $WHERE = $source eq 'dbSNP' ? "id_refsnp = ? and source='dbSNP'" : "id_ano=? and source='non-dbSNP'";
  my %SNPS = qw( 12 dbSNP 13 WI 14 HGBASE 15 TSC-CSHL 16 ANO );
  my $QUERY = "select internal_id, chr_start, chr_end, chr_strand, type, range_type,
                      validated, alleles, snpclass, mapweight, ambiguity, source,
                      id_refsnp, id_wi, id_hgbase, id_tsc, id_ano, chr_name
                 FROM snp
                WHERE $WHERE";

  my $sth = $self->prepare( $QUERY );
  eval { $sth->execute($refsnp_id);};
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
    my $STATUS = $arrayref->[6];
    $STATUS =~s/-/ /;
    $STATUS = ( $STATUS && $STATUS ne 'no info' ) ? "proven $STATUS" : 'suspected';
    my $snp = Bio::EnsEMBL::SNP->new_fast(
                  { 'dbID'       => $arrayref->[0],
                    '_snp_strand' => $arrayref->[3],
                    '_gsf_score'  => 1,
                    '_type'       => $arrayref->[4],
                    '_range_type' => $arrayref->[5],
                    '_validated'  => $arrayref->[6],
                    'status'     => $STATUS,
                    'alleles'    => $arrayref->[7],
                    '_ambiguity_code' => $arrayref->[10],
                    '_snpclass'   => $arrayref->[8],
                    '_mapweight'  => $arrayref->[9],
                    '_source'     => $arrayref->[11],
                    '_source_tag' => $arrayref->[11],
                    'link'        => \@links });
    return $snp;
  }
  return undef;
}

1;

