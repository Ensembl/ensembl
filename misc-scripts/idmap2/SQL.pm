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
    SELECT esi.stable_id as exon_stable, tsi.stable_id as transcript_stable, 
           gsi.stable_id as gene_stable, e.exon_id as exon_id, 
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
     ORDER by gene_stable, transcript_stable, exon_stable, e.sticky_rank
" );
  $sth->execute();

  return _store_in_arrayref( $sth );
}

sub target_exon_information {
  my $dbh = shift;
  
  my $sth = $dbh->prepare( "
    SELECT e.exon_id, t.transcript_id, t.gene_id,
           cl.embl_id as clone_id, cl.embl_version as clone_version, 
           c.offset as contig_offset, c.length as contig_length, e.exon_id as exon_id,
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
" );

  $sth->execute ();
  return _store_in_arrayref( $sth );
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

1;
