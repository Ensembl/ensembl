use strict;
use warnings;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::MiscFeature;
use Bio::EnsEMBL::Attribute;
use Bio::EnsEMBL::MiscSet;

use DBI;

use Getopt::Long;


my ( $host, $user, $pass, $port, $dbname,
     $mhost, $muser, $mpass, $mport, $mdbname);


GetOptions( "host=s", \$host,
            "user=s", \$user,
            "pass=s", \$pass,
            "port=i", \$port,
            "dbname=s", \$dbname,
            "mhost=s", \$mhost,
            "muser=s", \$muser,
            "mpass=s", \$mpass,
            "mport=i", \$mport,
            "mdbname=s", \$mdbname
          );



($mdbname && $mhost && $muser && $user && $mhost && $mdbname) || usage();

my $mart = DBI->connect("DBI:mysql:dbname=$mdbname;port=$mport;host=$mhost",
                        $muser, $mpass, {'RaiseError' => 1});


my $core = Bio::EnsEMBL::DBSQL::DBAdaptor->new
  (-host => $host,
   -user => $user,
   -pass => $pass,
   -port => $port,
   -dbname => $dbname);


my $set = Bio::EnsEMBL::MiscSet->new
  (-longest_feature => 2e6,
   -description     => 'ENCODE Regions',
   -code            => 'encode_regions',
   -name            => 'ENCODE Regions');


my $slice_adaptor = $core->get_SliceAdaptor();
my $misc_feat_adaptor = $core->get_MiscFeatureAdaptor();

my $rows = $mart->selectall_arrayref
  (qq{SELECT glook_encode_region_name, filt_chr_name, filt_chrom_start,
             filt_chrom_end, silent_type, encode_description
      FROM hsapiens__encode__look});


foreach my $row (@$rows) {
  my $slice = $slice_adaptor->fetch_by_region('chromosome', $row->[1]);

  my $name_attrib = Bio::EnsEMBL::Attribute->new
    (-CODE => 'name',
     -NAME => 'Name',
    -VALUE => $row->[0]);

  my $desc_attrib = Bio::EnsEMBL::Attribute->new
    (-CODE => 'type',
     -NAME => 'Type of Feature',
     -VALUE => $row->[4]);

  my $misc_feature = Bio::EnsEMBL::MiscFeature->new
    (-START  => $row->[2],
     -END    => $row->[3],
     -STRAND => 1,
     -SLICE  => $slice);

  $misc_feature->add_Attribute($name_attrib);
  $misc_feature->add_Attribute($desc_attrib);
  $misc_feature->add_MiscSet($set);

  $misc_feat_adaptor->store($misc_feature);
}



sub usage {
  print STDERR
    qq{
usage: \n
 perl import_encode_regions -mhost <mart_host> -muser <mart_user>
        -mdbname <mart_dbname> [-mpass <mart_pass> -mport <mart_port>]
        -host <ensembl_host> -dbname <ensembl_dbname> -user <ensembl_user>
        [-port <ensembl_port> -pass <ensembl_pass>]
};
  exit;
}
