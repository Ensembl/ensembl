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

    --gene_stable_id_file|gsi_file|gsi=FILE
                                        the path of the file containing a list
                                        of gene stable Ids that were deleted
    --transcript_stable_id_file|tsi_file|tsi=FILE    
                                        (optional) the path of the file
                                        containing a list of transcript stable
                                        Ids that were deleted
    --skip_ncrna|skip_ncRNA|skip_nc=0|1 (optionally) skip ncRNAs
    --skip_biotypes=LIST                (optionally) skip LISTed biotypes

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
  'gene_stable_id_file|gsi_file|gsi=s',
  'transcript_stable_id_file|tsi_file|tsi=s',
  'skip_ncrna|skip_ncRNA|skip_nc=s',
  'skip_biotypes=s@'
);
$support->allowed_params(
  $support->get_common_params,
  'althost',
  'altport',
  'altuser',
  'altpass',
  'altdbname',
  'gene_stable_id_file',
  'transcript_stable_id_file',
  'skip_ncrna',
  'skip_biotypes',
);

if ($support->param('help') or $support->error) {
    warn $support->error if $support->error;
    pod2usage(1);
}

$support->comma_to_list('skip_biotypes');

# ask user to confirm parameters to proceed
$support->confirm_params;

# get log filehandle and print heading and parameters to logfile
$support->init_log;

$support->check_required_params(
  'althost',
  'altport',
  'altuser',
  'altdbname',
  'gene_stable_id_file',
);

# connect to database and get adaptors
my $dba_new = $support->get_database('core');
my $dba_old = $support->get_database('ensembl', 'alt');
my $dbh_new = $dba_new->dbc->db_handle;

# define some globally used variables
my %genes;
my %gsi;
my %tsi;
my %tlsi;
my $gsi_string;
my $tsi_string;
my $tlsi_string;

# read list of deleted gene and transcript stable IDs from file(s)
&parse_deleted_files;

# create new mapping session
my $mapping_session_id = &create_mapping_session;

# create stable_id_event entries for all objects, mapping to themselves
&create_stable_id_events;

# set stable_id_event.new_stable_id to NULL for deleted objects
&mark_deleted;

# populate gene_archive and peptide_archive
&populated_archive;

# finish logfile
$support->finish_log;


### END main ###

sub parse_deleted_files {

  my $ta = $dba_old->get_TranscriptAdaptor;
  my $ga = $dba_old->get_GeneAdaptor;

  #
  # read list of deleted gene_stable_ids from file
  #
  $support->log_stamped("Reading list of deleted gene_stable_ids from file, and fetching associated transcript and translation stable IDs from the db...\n");
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

    # fetch associated transcript and translation stable IDs from the old db
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

  $gsi_string = "'".join("', '", keys(%gsi))."'";
  $tsi_string = "'".join("', '", keys(%tsi))."'";
  $tlsi_string = "'".join("', '", keys(%tlsi))."'";

  $support->log_stamped("Done loading ".scalar(keys(%gsi))." gene, ".scalar(keys(%tsi))." transcript and ".scalar(keys(%tlsi))." translation stable IDs.\n\n");

}


sub create_mapping_session {

  # create a new mapping session
  $support->log("Creating new mapping session...\n");
  my $old_db_name = $support->param('altdbname');
  my $new_db_name = $support->param('dbname');
  my $sql = qq(
    INSERT INTO mapping_session (old_db_name, new_db_name, created)
    VALUES ('$old_db_name', '$new_db_name', NOW())
  );
  my $c = $dbh_new->do($sql) unless ($support->param('dry_run'));
  my $mapping_session_id = $dbh_new->{'mysql_insertid'};
  $support->log("Done.\n\n");

  return $mapping_session_id;
}


sub create_stable_id_events {

  #
  # create stable_id_event entries for all objects, mapping to themselves
  #
  $support->log_stamped("Creating stable_id_event entries for all objects, mapping to themselves...\n");

  my $ga = $dba_old->get_GeneAdaptor;
  my @genes = @{ $ga->fetch_all };

  my $sql = qq(
    INSERT INTO stable_id_event 
    VALUES (?, ?, ?, ?, ?, ?, ?)
  );
  my $sth = $dbh_new->prepare($sql);

  # optionally skip ncRNAs
  #
  # this is the complete list of ncRNA biotype; you might need to update it (the
  # code below will try to help you with this)
  my @nc_biotypes = qw(
    miRNA
    miRNA_pseudogene
    misc_RNA
    misc_RNA_pseudogene
    Mt_rRNA
    Mt_tRNA
    Mt_tRNA_pseudogene
    rRNA
    rRNA_pseudogene
    scRNA
    scRNA_pseudogene
    snoRNA
    snoRNA_pseudogene
    snRNA
    snRNA_pseudogene
    tRNA_pseudogene
  );

  my %skip_biotypes = ();

  if ($support->param('skip_ncrna')) {

    %skip_biotypes = map { $_ => 1 } @nc_biotypes;

    # make sure we have a complete list of ncRNA biotypes
    $sql = qq(SELECT DISTINCT biotype from gene);
    my $sth1 = $dbh_new->prepare($sql);
    $sth1->execute;
    my @biotypes_db;
    
    while ((my $biotype) = $sth1->fetchrow_array) {
      push @biotypes_db, $biotype unless ($skip_biotypes{$biotype});
    }

    $sth1->finish;

    if (@biotypes_db) {
      print "These are the non-ncRNA biotypes found in the db:\n";
      map { print "  $_\n" } @biotypes_db;
      print "\nPlease check that the list of ncRNA biotypes is still complete, otherwise adapt the script.\n";
      exit unless $support->user_proceed("Continue?");
    }
  }

  # optionally skip other biotypes
  if ($support->param('skip_biotypes')) {
    %skip_biotypes = map { $_ => 1 } $support->param('skip_biotypes');
  }

  my %stats = map { $_ => 0 } qw(g tr tl g_tot tr_tot);
  my $num_genes = scalar(@genes);
  my $i;

  while (my $gene = shift(@genes)) {

    $support->log_progress($num_genes, ++$i);
  
    $stats{g_tot}++;
    
    next if ($skip_biotypes{$gene->biotype});

    unless ($support->param('dry_run')) {
      $sth->execute(
        $gene->stable_id,
        $gene->version,
        $gene->stable_id,
        $gene->version,
        $mapping_session_id,
        'gene',
        1
      );
    }

    $stats{g}++;

    # transcripts
    my @transcripts = @{ $gene->get_all_Transcripts };
    while (my $tr = shift(@transcripts)) {
      
      $stats{tr_tot}++;
    
      next if ($skip_biotypes{$tr->biotype});

      unless ($support->param('dry_run')) {
        $sth->execute(
          $tr->stable_id,
          $tr->version,
          $tr->stable_id,
          $tr->version,
          $mapping_session_id,
          'transcript',
          1
        );
      }

      $stats{tr}++;

      # translations
      if (my $tl = $tr->translation) {
        
        unless ($support->param('dry_run')) {
          $sth->execute(
            $tl->stable_id,
            $tl->version,
            $tl->stable_id,
            $tl->version,
            $mapping_session_id,
            'translation',
            1
          );
        }

        $stats{tl}++;
      }

    }

  }

  $sth->finish;

  $support->log_stamped("Done inserting entries for $stats{g} (of $stats{g_tot}) genes, $stats{tr} (of $stats{tr_tot}) transcripts, $stats{tl} translations.\n\n");

}

sub mark_deleted {

  # set stable_id_event.new_stable_id to NULL for deleted objects
  $support->log_stamped("Setting new_stable_id to NULL for deleted objects...\n");

  $support->log("Genes... ", 1);
  my $sql = qq(
    UPDATE stable_id_event
    SET new_stable_id = NULL, new_version = 0
    WHERE new_stable_id IN ($gsi_string)
    AND mapping_session_id = $mapping_session_id
  );
  my $c = $dbh_new->do($sql) unless ($support->param('dry_run'));
  $support->log("[$c]\n");

  $support->log("Transcripts... ", 1);
  $sql = qq(
    UPDATE stable_id_event
    SET new_stable_id = NULL, new_version = 0
    WHERE new_stable_id IN ($tsi_string)
    AND mapping_session_id = $mapping_session_id
  );
  $c = $dbh_new->do($sql) unless ($support->param('dry_run'));
  $support->log("[$c]\n");

  $support->log("Translations... ", 1);
  $sql = qq(
    UPDATE stable_id_event
    SET new_stable_id = NULL, new_version = 0
    WHERE new_stable_id IN ($tlsi_string)
    AND mapping_session_id = $mapping_session_id
  );
  $c = $dbh_new->do($sql) unless ($support->param('dry_run'));
  $support->log("[$c]\n");

  $support->log_stamped("Done.\n\n");

}


sub populated_archive {

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

  my $c = 0;

  foreach my $gsi (keys(%genes)) {
    my $gene = $genes{$gsi};

    foreach my $trans (@{ $gene->get_all_Transcripts }) {
    
      # skip transcripts that were not deleted (since %genes may contain genes
      # where only some but not all transcripts were deleted)
      next unless ($tsi{$trans->stable_id});
    
      my $tl = $trans->translation;

      # add peptide_archive entry
      unless ($support->param('dry_run')) {
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
      }

      $c++;
    }
  }
  $support->log_stamped("Done adding $c entries.\n\n");

}

