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
    --altdbname=NAME                    use old database NAME
    --althost=HOST                      use old database host HOST
    --altport=PORT                      use old database port PORT
    --altuser=USER                      use old database username USER
    --altpass=PASS                      use old database passwort PASS

    --mapping_session_id|msi=ID         latest mapping session
    --gene_stable_id_file|gsi_file|gsi=FILE
                                        the path of the file containing a list
                                        of gene stable Ids that were deleted
    --transcript_stable_id_file|tsi_file|tsi=FILE    
                                        (optional) the path of the file
                                        containing a list of transcript stable
                                        Ids that were deleted
    --skip_ncrna|skip_ncRNA|skip_nc=0|1 (optionally) skip ncRNAs

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
entries in gene_archive and peptide_archive are created. All this is done to
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
  'althost=s',
  'altport=n',
  'altuser=s',
  'altpass=s',
  'altdbname=s',
  'mapping_session_id|msi=n',
  'gene_stable_id_file|gsi_file|gsi=s',
  'transcript_stable_id_file|tsi_file|tsi=s',
  'skip_ncrna|skip_ncRNA|skip_nc=s'
);
$support->allowed_params(
  $support->get_common_params,
  'althost',
  'altport',
  'altuser',
  'altpass',
  'altdbname',
  'mapping_session_id',
  'gene_stable_id_file',
  'transcript_stable_id_file',
  'skip_ncrna',
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
  'althost',
  'altport',
  'altuser',
  'altdbname',
  'mapping_session_id',
  'gene_stable_id_file',
);

# connect to database and get adaptors
my $dba_new = $support->get_database('core');
my $dba_old = $support->get_database('ensembl', 'alt');
my $dbh_new = $dba_new->dbc->db_handle;
my $ta = $dba_old->get_TranscriptAdaptor;
my $ga = $dba_old->get_GeneAdaptor;

my $sql;
my $c;
my $sth;

#
# read list of deleted gene_stable_ids from file
#
$support->log_stamped("Reading list of deleted gene_stable_ids from file, and fetching associated transcript and translation stable IDs from the db...\n");
my %gsi;
my %tsi;
my %tlsi;
my %genes;
my $gfh = $support->filehandle('<', $support->param('gene_stable_id_file'));

while (my $g = <$gfh>) {
  chomp $g;
  my $gene = $ga->fetch_by_stable_id($g);

  # skip non-protein-coding genes
  unless ($gene->biotype eq 'protein_coding') {
    $support->log_warning("Gene ".$gene->stable_id." is non-protein_coding, skipping.\n", 1);
    next;
  }
  
  $genes{$g} = $gene;
  $gsi{$g} = 1;

  # fetch associated transcript and translation stable IDs from the 37 db
  foreach my $transcript (@{ $gene->get_all_Transcripts }) {
    $tsi{$transcript->stable_id} = 1;
    $tlsi{$transcript->translation->stable_id} = 1;
  }
}

#
# read list of deleted transcript_stable_ids from file
#
if ($support->param('transcript_stable_id_file')) {

  $support->log_stamped("Reading list of deleted transcript_stable_ids from file, and fetching associated translation stable IDs from the db...\n");

  my $tfh = $support->filehandle('<', $support->param('transcript_stable_id_file'));

  while (my $t = <$tfh>) {
    chomp $t;
    my $transcript = $ta->fetch_by_stable_id($t);

    # skip non-protein-coding genes
    unless ($transcript->biotype eq 'protein_coding') {
      $support->log_warning("Transcript ".$transcript->stable_id." is non-protein_coding, skipping.\n", 1);
      next;
    }
    
    $tsi{$transcript->stable_id} = 1;
    $tlsi{$transcript->translation->stable_id} = 1;
    
    my $gene = $ga->fetch_by_transcript_id($transcript->dbID);
    $genes{$gene->stable_id} = $gene;
  }
}

my $gsi_string = "'".join("', '", keys(%gsi))."'";
my $tsi_string = "'".join("', '", keys(%tsi))."'";
my $tlsi_string = "'".join("', '", keys(%tlsi))."'";

$support->log_stamped("Done loading ".scalar(keys(%gsi))." gene, ".scalar(keys(%tsi))." transcript and ".scalar(keys(%tlsi))." translation stable IDs.\n\n");

# exit now if doing a dry run
if ($support->param('dry_run')) {
  $support->log("Nothing else to do for a dry run. Exiting.\n\n");
  $support->finish_log;
  exit;
}

# create a new mapping session
$support->log("Creating new mapping session...\n");
my $old_db_name = $support->param('altdbname');
my $new_db_name = $support->param('dbname');
$sql = qq(
  INSERT INTO mapping_session (old_db_name, new_db_name, created)
  VALUES ('$old_db_name', '$new_db_name', NOW())
);
$c = $dbh_new->do($sql);
my $mapping_session_id = $dbh_new->{'mysql_insertid'};
$support->log("Done.\n\n");

#
# create stable_id_event entries for all objects, mapping to themselves
#
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
  VALUES (?, ?, ?, ?, ?, ?, ?)
);
$sth = $dbh_new->prepare($sql);

# optionally skip ncRNAs
my %nc_genes = ();
if ($support->param('skip_ncrna')) {
  foreach my $biotype (qw(miRNA misc_RNA Mt-tRNA Mt-rRNA rRNA snoRNA snRNA)) {
    map { $nc_genes{$_->stable_id} = 1 }
      @{ $ga->fetch_all_by_biotype($biotype) };
  }
}

foreach my $k (keys %unique_sie) {
  # optionally skip ncRNAs
  next if ($nc_genes{$unique_sie{$k}->{'new_stable_id'}});

  $sth->execute(
    $unique_sie{$k}->{'new_stable_id'},
    $unique_sie{$k}->{'new_version'},
    $unique_sie{$k}->{'new_stable_id'},
    $unique_sie{$k}->{'new_version'},
    $mapping_session_id,
    $unique_sie{$k}->{'type'},
    1
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

$c = 0;

foreach my $gsi (keys(%genes)) {
  my $gene = $genes{$gsi};

  foreach my $trans (@{ $gene->get_all_Transcripts }) {
  
    # skip transcripts that were not deleted (since %genes may contain genes
    # were only some but not all transcripts were deleted)
    next unless ($tsi{$trans->stable_id});
  
    my $tl = $trans->translation;

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

# finish logfile
$support->finish_log;

