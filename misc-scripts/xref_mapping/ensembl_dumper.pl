# Dump Ensembl sequences to fasta file

use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;

my $host   = 'ecs2';
my $port   = 3364;
my $user   = 'ensro';
my $dbname = 'homo_sapiens_core_25_34e';

my $file = "ensembl_transcripts.fasta";

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $host,
					    -port   => $port,
                                            -user   => $user,
                                            -dbname => $dbname);

my $slice_adaptor = $db->get_SliceAdaptor();
my $slice = $slice_adaptor->fetch_by_region('chromosome', '21');

my $t = 0;

open(FILE, ">" . $file);

foreach my $gene (@{$slice->get_all_Genes()}) {

  foreach my $trans (@{$gene->get_all_Transcripts()}) {

    print FILE ">" . $trans->dbID() . "\n" . $trans->spliced_seq() . "\n";
    $t++;

  }
}

close(FILE);

print "Wrote $t transcripts to $file\n";
