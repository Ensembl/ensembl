# This module encapsulates the SQL needed for the idmapping
# Author: Arne Stabenau 3.10.2001

# ( why not use EnsEMBL here? Currently to slow?? This piece of software is 
#  going to be replaced very soon by a java version )

use DBI;

package SQL;



# takes a database handle and produces all old exon information needed

sub orig_exon_information {
  my $dbh  = shift;
  my %res;

  my $sth = $dbh->prepare( "
    SELECT exon_stable, transcript_stable, 
           gene_stable, exon_id, gene_id, transcript_id, 
           clone_id, clone_version, 
           contig_offset, contig_length,
           exon_start, exon_end, exon_strand,
           chr_name, chr_start, chr_end, 
           raw_start, raw_end, raw_ori
     FROM  exon_temp 
    ORDER by clone_id
" );
#    LIMIT 10000" );
  $sth->execute();
  return _store_in_arrayref( $sth );
}

sub target_exon_information {
  my $dbh = shift;
  
  my $sth = $dbh->prepare( "
    SELECT exon_id, transcript_id, gene_id,
           clone_id, clone_version, 
           contig_offset, contig_length,
           exon_start, exon_end, exon_strand,
           chr_name, chr_start, chr_end, 
           raw_start, raw_end, raw_ori
     FROM exon_temp
    ORDER by clone_id
 " );
#   LIMIT 10000");
  $sth->execute();
  return _store_in_arrayref( $sth );
}



sub old_orig_exon_information {
  my $dbh = shift;
  my %res;
  $dbh->do( "drop table exon_temp" );

  my $sth = $dbh->prepare( "
    CREATE TABLE exon_temp
    SELECT esi.stable_id as exon_stable, tsi.stable_id as transcript_stable, 
           gsi.stable_id as gene_stable, e.exon_id as exon_id, t.gene_id as gene_id,
           t.transcript_id as transcript_id,
           cl.embl_id as clone_id, cl.embl_version as clone_version, 
           c.offset as contig_offset, c.length as contig_length,
           e.seq_start as exon_start, e.seq_end as exon_end, e.strand as exon_strand,
           sgp.chr_name as chr_name, sgp.chr_start as chr_start, sgp.chr_end as chr_end, 
           sgp.raw_start as raw_start, sgp.raw_end as raw_end, sgp.raw_ori as raw_ori
      FROM exon e, exon_transcript et, exon_stable_id esi, contig c, clone cl,
           static_golden_path sgp, transcript t, gene_stable_id gsi, 
           transcript_stable_id tsi
     WHERE e.exon_id = et.exon_id
       AND esi.exon_id = e.exon_id
       AND et.transcript_id = t.transcript_id
       AND tsi.transcript_id = et.transcript_id
       AND t.gene_id = gsi.gene_id
       AND c.internal_id = e.contig_id
       AND cl.internal_id = c.clone
       AND sgp.raw_id = c.internal_id
     ORDER by clone_id, clone_version,contig_offset, exon_strand
" );
  $sth->execute();
#     ORDER by gene_stable, transcript_stable, exon_stable, e.sticky_rank

#  return _store_in_arrayref( $sth );
}

sub old_target_exon_information {
  my $dbh = shift;
  $dbh->do( "drop table exon_temp" );
  my $sth = $dbh->prepare( "
    CREATE TABLE exon_temp
    SELECT e.exon_id, t.transcript_id, t.gene_id,
           cl.embl_id as clone_id, cl.embl_version as clone_version, 
           c.offset as contig_offset, c.length as contig_length,
           e.seq_start as exon_start, e.seq_end as exon_end, e.strand as exon_strand,
           sgp.chr_name as chr_name, sgp.chr_start as chr_start, sgp.chr_end as chr_end, 
           sgp.raw_start as raw_start, sgp.raw_end as raw_end, sgp.raw_ori as raw_ori
      FROM exon e, exon_transcript et, contig c, clone cl,
           static_golden_path sgp, transcript t
     WHERE e.exon_id = et.exon_id
       AND et.transcript_id = t.transcript_id
       AND c.internal_id = e.contig_id
       AND cl.internal_id = c.clone
       AND sgp.raw_id = c.internal_id
     ORDER by clone_id, clone_version,contig_offset, exon_strand
" );
  
  $sth->execute();

#  return _store_in_arrayref( $sth );
}


sub _store_in_arrayref {
  my $sth = shift;
  my $res = [];

  my $count = 0;

  while( my $hr = $sth->fetchrow_hashref() ) {
    my %tempHash = %$hr;
    push( @$res, \%tempHash );
  }
  return $res;
}

sub exon_sequence {
  my $dbh = shift;
  my $exon_id = shift;
  my $exonSeq;

  my $sth = $dbh->prepare( "
    SELECT e.exon_id, substring( d.sequence, e.seq_start, e.seq_end-e.seq_start+1 ), 
           e.seq_start as start, e.seq_end as end,
           e.strand as strand, e.sticky_rank as rank
      FROM exon e, contig c, dna d
     WHERE e.contig_id = c.internal_id
       AND c.dna = d.id
       AND e.exon_id = ?
     ORDER by rank
  " );
  $sth->execute( $exon_id );
  
  while( my $arrRef = $sth->fetchrow_arrayref ) {
    my $part  = $arrRef->[1];
    if( $arrRef->[4] != 1 ) {
      $part = reverse( $part );
      $part =~ tr/ACGTacgt/TGCATGCA/;
    }
    $exonSeq .= $part;
  }

  return $exonSeq;
}


1;
