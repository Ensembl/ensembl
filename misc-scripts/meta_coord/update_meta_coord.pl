#!/usr/local/ensembl/bin/perl -w

use strict;
use Bio::EnsEMBL::DBSQL::DBConnection;
use Getopt::Long;

my $help = 0;
my ($host, $port, $user, $pass, $dbname);

my $usage = "\n$0 --host ecs2 --port 3364 --user ensadmin --pass ensembl --dbname homo_sapiens_core_27_35a

[--help] displays this menu.

This script will dump the current meta_coord table in homo_sapiens_core_27_35a.meta_coord file.
Then it will update the meta_coord table for all the following table names one by one
                     assembly_exception
                     gene
                     exon
                     density_feature
                     ditag
                     ditag_feature
                     dna_align_feature
                     karyotype
                     marker_feature
                     misc_feature
                     qtl_feature
                     oligo_feature
                     prediction_exon
                     prediction_transcript
                     protein_align_feature
                     regulatory_feature
                     regulatory_search_region
                     repeat_feature
                     simple_feature
                     transcript\n";

GetOptions('help'     => \$help,
           'host=s'   => \$host,
           'port=i'   => \$port,
           'user=s'   => \$user,
           'pass=s'   => \$pass,
           'dbname=s' => \$dbname);

#print "help: $help argv:" . scalar(@ARGV) . "$host $port $user $pass $dbname\n";

#if ($help || scalar @ARGV == 0) {
#  print $usage,"\n";
#  exit 0;
#}

my $dbc = new Bio::EnsEMBL::DBSQL::DBConnection(-host   => $host,
                                                -port   => $port,
                                                -user   => $user,
                                                -pass   => $pass,
                                                -dbname => $dbname);

my @table_names = qw(
                     assembly_exception
                     gene
                     exon
                     density_feature
                     ditag
                     ditag_feature
                     dna_align_feature
                     karyotype
                     marker_feature
                     misc_feature
                     qtl_feature
                     oligo_feature
                     prediction_exon
                     prediction_transcript
                     protein_align_feature
                     regulatory_feature
                     regulatory_search_region
                     repeat_feature
                     simple_feature
                     transcript
                    );

unless (system("mysql -h$host -P$port -u$user -p$pass -N -e 'select * from meta_coord' $dbname > $dbname.meta_coord.backup") ==0) {
  print STDERR "Can't dump the original meta_coord for back up\n";
  exit 1;
} else {
  print STDERR "original meta_coord table backed up in $dbname.meta_coord.backup\n";
}

foreach my $table_name (@table_names) {
  print STDERR "Updating $table_name table entries...";
  my $sql = "delete from meta_coord where table_name = ? ";
  my $sth = $dbc->prepare($sql);
  $sth->execute($table_name);
  $sth->finish;

  $sql = "insert into meta_coord ".
         "select '$table_name',s.coord_system_id, max(t.seq_region_end-t.seq_region_start+1) ".
         "from $table_name t, seq_region s ".
         "where t.seq_region_id=s.seq_region_id group by s.coord_system_id";
  $sth = $dbc->prepare($sql);
  $sth->execute;
  $sth->finish;
  print STDERR "Done\n";
}

