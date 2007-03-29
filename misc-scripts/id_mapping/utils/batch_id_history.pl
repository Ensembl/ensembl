#!/usr/local/bin/perl

=head1 NAME

batch_id_history.pl - find stable IDs in the archive

=head1 SYNOPSIS

find_stable_ids_in_archive.pl [arguments]

Required arguments:

  --host=hOST                 database host HOST
  --port=PORT                 database port PORT
  --user=USER                 database username USER
  --dbname=NAME               database name NAME
  --stable_id_file=FILE       read stable ID list from FILE

Optional arguments:

  --pass=PASS                 database passwort PASS
  --outfile=FILE              write output to FILE
  --pep_seq                   print peptide sequence


=head1 DESCRIPTION

This script reads a list of stable IDs from a file and sees if it can find them
in the stable ID archive. It will print the ID history for each of them and
optinally the peptide sequence found there as well. Note that this will not 
print the full history network, but rather branch out from your focus stable ID
only. If you are interested in the full network, have a look at
Bio::EnsEMBL::StableIdHistoryTree and related modules.

=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http://www.ensembl.org/info/about/code_licence.html for details.

=head1 AUTHOR

Patrick Meidl <meidl@ebi.ac.uk>, Ensembl core API team

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<ensembl-dev@ebi.ac.uk>

=cut

use strict;
use warnings;
no warnings 'uninitialized';

use Getopt::Long;
use Bio::EnsEMBL::SimpleFeature;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

$| = 1;

my ($host, $port, $user, $pass, $dbname, $stable_id_file, $outfile, $pep_seq);

GetOptions(
    "host=s",             \$host,
    "port=i",             \$port,
    "user=s",             \$user,
    "pass=s",             \$pass,
    "dbname=s",           \$dbname,
    "stable_id_file=s",   \$stable_id_file,
    "outfile=s",          \$outfile,
    "pep_seq",            \$pep_seq,
);

# check required params
unless ($host && $port && $user && $dbname && $stable_id_file) {
  die "ERROR: Unable to run script.\nNeed host, port, user, dbname and stable_id_file parameters.\n";
}

# connect to database and get adaptors
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -HOST     => $host,
    -PORT     => $port,
    -USER     => $user,
    -PASS     => $pass,
    -DBNAME   => $dbname,
);

my $aa = $db->get_ArchiveStableIdAdaptor;

# read list of stable IDs from file
my $infh;
open($infh, "<", $stable_id_file) or
  die("Can't open $stable_id_file for reading: $!");

# get output filehandle
my $outfh;
if ($outfile) {
  open($outfh, ">", $outfile) or die("Can't open $outfile for writing: $!");
} else {
  $outfh = \*STDOUT;
}

while (my $sid = <$infh>) {
  
  # skip comments
  next if (/^#/);

  chomp($sid);
  print $outfh "\n$sid\n";

  # fetch initial archive id
  my $archiveStableId = $aa->fetch_by_stable_id($sid);

  if ($archiveStableId) {
  
    # get archive id history
    my $history = $aa->fetch_archive_id_history($archiveStableId);

    # group archive ids by release
    my $release;
    foreach my $arch_id (@$history) {
      push @{ $release->{$arch_id->release} }, $arch_id;
    }

    # loop over releases and print results
    foreach my $r (sort { $a <=> $b } keys %$release) {
      
      print $outfh "\n  Release $r (".$release->{$r}->[0]->assembly.", ".
        $release->{$r}->[0]->db_name.")\n";

      # loop over archive ids
      foreach my $a (@{ $release->{$r} }) {
        print $outfh "    ".$a->stable_id.".".$a->version."\n";
        
        # print peptide sequence
        if ($pep_seq) {
          foreach my $pep (@{ $a->get_all_translation_archive_ids }) {
            print $outfh "      Peptide ".$pep->stable_id."> ".
              $pep->get_peptide."\n";
          }
        }

      }
    }
    
  } else {
    print $outfh "  Not found in archive.\n\n";
  }

}

close($infh);
close($outfh);


