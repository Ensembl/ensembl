use strict;
use warnings;

use Bio::EnsEMBL::DBSQL::DBAdaptor;

my $user = 'ensro';
my $pass = undef;
my $dbname = 'mcvicker_chimp_agp';
my $port = 3350;
my $host = 'ecs4';

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new
  (-user => $user,
   -host => $host,
   -pass => $pass,
   -dbname => $dbname,
   -port => $port);

my $slice_adaptor = $db->get_SliceAdaptor();
my $csa = $db->get_CoordSystemAdaptor();
my $from_cs = $csa->fetch_by_name("scaffold", "BROAD1");
my $to_cs   = $csa->fetch_by_name("chromosome", "BROAD1");

my $mapper = $db->get_AssemblyMapperAdaptor()->fetch_by_CoordSystems($from_cs,
                                                                     $to_cs);

while(<>) {
  my($asm_sr_id, $cmp_sr_id,$asm_start, $asm_end,
     $cmp_start, $cmp_end, $ori) = split;

  my $slice = $slice_adaptor->fetch_by_seq_region_id($cmp_sr_id);

  if(!$slice) {
    die("unknown seq_region $asm_sr_id");
  }

  my $cmp_name = $slice->seq_region_name();

  my($new_name, $new_start, $new_end, $new_strand) =
    $mapper->fastmap($cmp_name, $cmp_start, $cmp_end, 1, $from_cs);

  if(!$new_name) {
    die("Mapping failed");
  }

  $slice = $slice_adaptor->fetch_by_region('chromosome', $new_name,
                                           undef,undef,undef, "BROAD1");

  my $new_sr_id = $slice->get_seq_region_id();

  $ori *= $new_strand;

  print join "\t", $asm_sr_id, $new_sr_id, $asm_start, $asm_end, $new_start,
    $new_end, $ori, "\n";
}
