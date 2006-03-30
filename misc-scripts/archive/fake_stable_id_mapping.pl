#!/usr/local/bin/perl

=head1 NAME

fake_stable_id_mapping.pl - fix the stable ID archive for a database after some
genes were deleted

=head1 SYNOPSIS

fake_stable_id_mapping.pl [options]

Options:

    --conffile, --conf=FILE             read parameters from FILE
                                        (default: conf/Conversion.ini)

    --dbname, db_name=NAME              use new database NAME
    --host, --dbhost, --db_host=HOST    use new database host HOST
    --port, --dbport, --db_port=PORT    use new database port PORT
    --user, --dbuser, --db_user=USER    use new database username USER
    --pass, --dbpass, --db_pass=PASS    use new database passwort PASS
    --ensembldbname=NAME                use old database NAME
    --ensemblhost=HOST                  use old database host HOST
    --ensemblport=PORT                  use old database port PORT
    --ensembluser=USER                  use old database username USER
    --ensemblpass=PASS                  use old database passwort PASS

    --mapping_session_id=ID             latest mapping session
    --stable_id_file=FILE               the path of the file containing a list
                                        of gene stable Ids that were deleted

    --logfile, --log=FILE               log to FILE (default: *STDOUT)
    --logpath=PATH                      write logfile to PATH (default: .)
    --logappend, --log_append           append to logfile (default: truncate)
    -v, --verbose                       verbose logging (default: false)
    -i, --interactive=0|1               run script interactively (default: true)
    -n, --dry_run, --dry=0|1            don't write results to database
    -h, --help, -?                      print help (this message)

=head1 DESCRIPTION

This script fakes a stable ID archive run for a database where some genes were
deleted. A new mapping session is created and all stable IDs other than the
deleted ones are mapped to themselves. For the deleted genes, appropriate
entries in gene_archive and peptide archive are created. All this is done to
the new database, whereas stable Ids of deleted objects are looked up in the
old database (if you haven't deleted them yet, old and new can point to the
same db).

Please note that when using two different databases as input and one is from
the last release, you might have to use a cvs checkout of a matching branch.

=head1 LICENCE

This code is distributed under an Apache style licence:
Please see http://www.ensembl.org/code_licence.html for details

=head1 AUTHOR

Patrick Meidl <pm2@sanger.ac.uk>

=head1 CONTACT

Post questions to the EnsEMBL development list ensembl-dev@ebi.ac.uk

=cut

use strict;
use warnings;
no warnings 'uninitialized';

use FindBin qw($Bin);
use vars qw($SERVERROOT);

BEGIN {
    $SERVERROOT = "$Bin/../../..";
    unshift(@INC, "$SERVERROOT/ensembl/modules");
    unshift(@INC, "$SERVERROOT/bioperl-live");
}

use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Utils::ConversionSupport;
use Digest::MD5 qw(md5_hex);

$| = 1;

my $support = new Bio::EnsEMBL::Utils::ConversionSupport($SERVERROOT);

# parse options
$support->parse_common_options(@_);
$support->parse_extra_options(
  'ensemblhost=s',
  'ensemblport=n',
  'ensembluser=s',
  'ensemblpass=s',
  'ensembldbname=s',
  'mapping_session_id=n',
  'stable_id_file=s',
);
$support->allowed_params(
  $support->get_common_params,
  'ensemblhost',
  'ensemblport',
  'ensembluser',
  'ensemblpass',
  'ensembldbname',
  'mapping_session_id',
  'stable_id_file',
);

if ($support->param('help') or $support->error) {
    warn $support->error if $support->error;
    pod2usage(1);
}

# ask user to confirm parameters to proceed
$support->confirm_params;

# get log filehandle and print heading and parameters to logfile
$support->init_log;

$support->check_required_params(
  'ensemblhost',
  'ensemblport',
  'ensembluser',
  'ensemblpass',
  'ensembldbname',
  'mapping_session_id',
  'stable_id_file',
);

# connect to database and get adaptors
my $dba_new = $support->get_database('core');
my $dba_old = $support->get_database('ensembl', 'ensembl');
my $dbh_new = $dba_new->dbc->db_handle;
my $ta = $dba_old->get_TranscriptAdaptor;
my $ga = $dba_old->get_GeneAdaptor;

my $sql;
my $c;
my $sth;

# read list of deleted gene_stable_ids from file
$support->log_stamped("Reading list of deleted gene_stable_ids from file, and fetching associated transcript and translation stable IDs from the db...\n");
my @gsi;
my @tsi;
my @tlsi;
my @genes;
my $infh = $support->filehandle('<', $support->param('stable_id_file'));

while (my $g = <$infh>) {
  chomp $g;
  my $gene = $ga->fetch_by_stable_id($g);

  # skip non-protein-coding genes
  unless ($gene->biotype eq 'protein_coding') {
    $support->log_warning("Gene ".$gene->stable_id." is non-protein_coding, skipping.\n", 1);
    next;
  }
  
  push @gsi, $g;
  push @genes, $gene;

  # fetch associated transcript and translation stable IDs from the 37 db
  foreach my $transcript (@{ $gene->get_all_Transcripts }) {
    push @tsi, $transcript->stable_id;
    push @tlsi, $transcript->translation->stable_id;
  }
}

my $gsi_string = "'".join("', '", @gsi)."'";
my $tsi_string = "'".join("', '", @tsi)."'";
my $tlsi_string = "'".join("', '", @tlsi)."'";
$support->log_stamped("Done reading ".scalar(@gsi)." gene and fetching ".scalar(@tsi)." transcript and ".scalar(@tlsi)." translation stable IDs.\n\n");

# exit now if doing a dry run
if ($support->param('dry_run')) {
  $support->log("Nothing else to do for a dry run. Exiting.\n\n");
  $support->finish_log;
  exit;
}

# backup archive tables in case you screw up
$support->log_stamped("Creating backup of stable_id_event, gene_archive and peptide_archive...\n");
$sql = qq(
  CREATE TABLE stable_id_event_bak
  SELECT * FROM stable_id_event
);
$c = $dbh_new->do($sql);
$sql = qq(
  CREATE TABLE gene_archive_bak
  SELECT * FROM gene_archive
);
$c = $dbh_new->do($sql);
$sql = qq(
  CREATE TABLE peptide_archive_bak
  SELECT * FROM peptide_archive
);
$c = $dbh_new->do($sql);
$support->log_stamped("Done.\n\n");

# create a new mapping session
$support->log("Creating new mapping session...\n");
my $old_db_name = $support->param('ensembldbname');
my $new_db_name = $support->param('dbname');
$sql = qq(
  INSERT INTO mapping_session (old_db_name, new_db_name, created)
  VALUES ('$old_db_name', '$new_db_name', NOW())
);
$c = $dbh_new->do($sql);
my $mapping_session_id = $dbh_new->{'mysql_insertid'};
$support->log("Done.\n\n");

# create stable_id_event entries for all objects, mapping to themselves
$support->log_stamped("Creating stable_id_event entries for all objects, mapping to themselves...\n");
my $msi = $support->param('mapping_session_id');
$sql = qq(
  SELECT new_stable_id, new_version, mapping_session_id, type
  FROM stable_id_event
  WHERE mapping_session_id = $msi
  AND new_stable_id IS NOT NULL
);
$sth = $dbh_new->prepare($sql);
$sth->execute;

$c = 0;
my %unique_sie;

while (my $r = $sth->fetchrow_hashref) {
  $unique_sie{$r->{'new_stable_id'}.":".$r->{'new_version'}} = $r;
  $c++;
}

$sth->finish;

$sql = qq(
  INSERT INTO stable_id_event 
  VALUES (?, ?, ?, ?, ?, ?)
);
$sth = $dbh_new->prepare($sql);

foreach my $k (keys %unique_sie) {
  $sth->execute(
    $unique_sie{$k}->{'new_stable_id'},
    $unique_sie{$k}->{'new_version'},
    $unique_sie{$k}->{'new_stable_id'},
    $unique_sie{$k}->{'new_version'},
    $mapping_session_id,
    $unique_sie{$k}->{'type'},
  );
}

$sth->finish;

my $u = scalar(keys(%unique_sie));
$support->log_stamped("Done inserting $u entries (ignoring ".($c-$u)." duplicates).\n\n");

# set stable_id_event.new_stable_id to NULL for deleted objects
$support->log_stamped("Setting new_stable_id to NULL for deleted objects...\n");

$support->log("Genes... ", 1);
$sql = qq(
  UPDATE stable_id_event
  SET new_stable_id = NULL, new_version = 0
  WHERE new_stable_id IN ($gsi_string)
  AND mapping_session_id = $mapping_session_id
);
$c = $dbh_new->do($sql);
$support->log("[$c]\n");

$support->log("Transcripts... ", 1);
$sql = qq(
  UPDATE stable_id_event
  SET new_stable_id = NULL, new_version = 0
  WHERE new_stable_id IN ($tsi_string)
  AND mapping_session_id = $mapping_session_id
);
$c = $dbh_new->do($sql);
$support->log("[$c]\n");

$support->log("Translations... ", 1);
$sql = qq(
  UPDATE stable_id_event
  SET new_stable_id = NULL, new_version = 0
  WHERE new_stable_id IN ($tlsi_string)
  AND mapping_session_id = $mapping_session_id
);
$c = $dbh_new->do($sql);
$support->log("[$c]\n");

$support->log_stamped("Done.\n\n");

# populate gene_archive and peptide_archive
$support->log_stamped("Populating gene_archive and peptide_archive...\n");

my $sth_gene = $dbh_new->prepare(qq(
  INSERT INTO gene_archive
  VALUES (?, ?, ?, ?, ?, ?, ?, ?)
));
my $sth_pep = $dbh_new->prepare(qq(
  INSERT INTO peptide_archive (md5_checksum, peptide_seq)
  VALUES (?, ?)
));
my $sth_sie = $dbh_new->prepare(qq(
  DELETE sie FROM stable_id_event sie, mapping_session ms
  WHERE sie.new_stable_id = ?
  AND sie.new_version = ?
  AND sie.mapping_session_id = ms.mapping_session_id
  AND ms.old_db_name = 'ALL'
));

$c = 0;

foreach my $gene (@genes) {
  # delete ALL mapping session entries from stable_id_event where gene was
  # deleted
  $sth_sie->execute($gene->stable_id, $gene->version);

  foreach my $trans (@{ $gene->get_all_Transcripts }) {
    my $tl = $trans->translation;

    # delete ALL mapping session entries from stable_id_event where object was
    # deleted
    $sth_sie->execute($trans->stable_id, $trans->version);
    $sth_sie->execute($tl->stable_id, $tl->version);
  
    # add peptide_archive entry
    $sth_pep->execute(md5_hex($trans->translate->seq), $trans->translate->seq);
    my $pid = $dbh_new->{'mysql_insertid'};

    # add gene_archive entry
    $sth_gene->execute(
        $gene->stable_id,
        $gene->version,
        $trans->stable_id,
        $trans->version,
        $tl->stable_id,
        $tl->version,
        $pid,
        $mapping_session_id
    );

    $c++;
  }
}
$support->log_stamped("Done adding $c entries.\n\n");

# delete backup tables? (or remind to do so)
if ($support->user_proceed("Would you like to drop the temporary backup table?")) {
  $dbh_new->do(qq(DROP TABLE stable_id_event_bak));
  $dbh_new->do(qq(DROP TABLE gene_archive_bak));
  $dbh_new->do(qq(DROP TABLE peptide_archive_bak));
  $support->log("Done.\n\n");
}

# finish logfile
$support->finish_log;

