package XrefParser::RFAMCoreParser;

use strict;
use warnings;
use Carp;

use base qw( XrefParser::RFAMParser );
use Bio::EnsEMBL::Registry;

sub get_rfam_transcript_stable_ids {
  my ( $self, $mapper ) = @_;
  if(!defined $mapper) {
    croak "Need to connect to a core database; please use a mapper based configuration";
  }
  my $dba      = $mapper->core()->dba();
  my $rfam_sql = <<'SQL';
select distinct t.stable_id, hit_name 
from analysis a 
join transcript t on (a.analysis_id = t.analysis_id and a.logic_name = ? and t.biotype != ?) 
join exon_transcript et on (t.transcript_id = et.transcript_id) 
join supporting_feature sf on (et.exon_id = sf.exon_id and sf.feature_type = ?) 
join dna_align_feature df on (sf.feature_id = df.dna_align_feature_id)
order by hit_name
SQL

  my $sth = $dba->dbc->prepare($rfam_sql);
  $sth->execute( 'ncRNA', 'miRNA', 'dna_align_feature' );

#hash keyed on RFAM accessions, value is an array of ensembl transcript stable_ids
  my %rfam_transcript_stable_ids;

  while ( my ( $stable_id, $hit_name ) = $sth->fetchrow_array ) {
    my $rfam_id;
    if ( $hit_name =~ /^(RF\d+)/ ) {
      $rfam_id = $1;
    }
    if ($rfam_id) {
      push @{ $rfam_transcript_stable_ids{$rfam_id} }, $stable_id;
    }
  }
  $sth->finish;
  return \%rfam_transcript_stable_ids;
}

1;
