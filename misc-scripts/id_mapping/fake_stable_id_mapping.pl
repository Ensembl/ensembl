#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Don't change the above line.
# Change the PATH in the myRun.ksh script if you want to use another perl.

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
                                        (deprectated) the path of the file
                                        containing a list of gene stable Ids
                                        that were deleted
    --transcript_stable_id_file|tsi_file|tsi=FILE    
                                        (deprectated, optional) the path of the
                                        file containing a list of transcript
                                        stable Ids that were deleted
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
deleted.

It assumes that the new database already has its *_stable_id tables populated.
A new mapping session is created and all stable IDs other than the deleted ones
are mapped to themselves. For the deleted genes, appropriate entries in
gene_archive and peptide_archive are created. All this is done to the new
database, whereas stable Ids of deleted objects are looked up in the old
database. The scripts also increments the stable ID versions of genes where
transcripts were deleted (but the gene is still there due to other retained
transcripts).

Please note that when using two different databases as input and one is from
the last release, you might have to patch it to the current schema.


=head1 AUTHOR

Patrick Meidl <pm2@sanger.ac.uk>

=head1 CONTACT

Post questions to the EnsEMBL development list http://lists.ensembl.org/mailman/listinfo/dev

=cut

use strict;
use warnings;
no warnings 'uninitialized';

use Getopt::Long;
use Pod::Usage;
use FindBin qw($Bin);
use Bio::EnsEMBL::Utils::ConversionSupport;
use Digest::MD5 qw(md5_hex);

$| = 1;

my $support = new Bio::EnsEMBL::Utils::ConversionSupport("$Bin/../../..");

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
);

# connect to database and get adaptors
my $dba_new = $support->get_database('core');
my $dba_old = $support->get_database('ensembl', 'alt');
my $dbh_new = $dba_new->dbc->db_handle;

# define some globally used variables
my %genes = ();
my %genes_mod = ();
my %gsi_del = ();
my %gsi_mod = ();
my %tsi_del = ();
my %tlsi_del = ();
my $gsi_string;
my $gsi_mod_string;
my $tsi_string;
my $tlsi_string;
my %skip_biotypes = ();

#
# find out which genes and transcripts were deleted
#
if ($support->param('gene_stable_id_file') or
    $support->param('transcript_stable_id_file')) {

  # read list of deleted gene and transcript stable IDs from file(s)
  # this is error-prone and therefore DEPRECATED!
  &parse_deleted_files;

} else {
  # infer deleted objects from dbs (more robust)
  &determine_deleted;
}

# create new mapping session
my $mapping_session_id = &create_mapping_session;

# create stable_id_event entries for all objects, mapping to themselves
&create_stable_id_events;

# increment gene version for all genes where transcripts were deleted
&increment_gene_versions;

# set stable_id_event.new_stable_id to NULL for deleted objects
&mark_deleted;

# populate gene_archive and peptide_archive
&populated_archive;

# finish logfile
$support->finish_log;


### END main ###


=head2 parse_deleted_files

  Example     : &parse_deleted_files;
  Description : DEPRECATED
                Read list of deleted gene and transcript stable IDs from file(s).
                Note that this method of determining which objects were deleted
                is now DEPRECATED (because it was error-prone when dealing with 
                both whole gene and individual transcript deletions).
  Return type : none
  Exceptions  : thrown on missing files
  Caller      : main()
  Status      : At Risk
              : under development

=cut

sub parse_deleted_files {

  $support->log_warning("DEPRECATED. Don't use stable ID files (this is error-prone), rather let the script determine which objects were deleted from the dbs.\n");

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

    $genes{$g} = $gene;
    $gsi_del{$g} = 1;

    # fetch associated transcript and translation stable IDs from the old db
    foreach my $transcript (@{ $gene->get_all_Transcripts }) {
      $tsi_del{$transcript->stable_id} = 1;
      my $tl = $transcript->translation;
      $tlsi_del{$tl->stable_id} = 1 if ($tl);
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
      
      $tsi_del{$transcript->stable_id} = 1;
      my $tl = $transcript->translation;
      $tlsi_del{$tl->stable_id} = 1 if ($tl);
      
      my $gene = $ga->fetch_by_transcript_id($transcript->dbID);
      $genes{$gene->stable_id} = $gene;
      $gsi_mod{$gene->stable_id} = 1;
    }
  }

  $gsi_string = "'".join("', '", keys(%gsi_del))."'";
  $gsi_mod_string = "'".join("', '", keys(%gsi_mod))."'";
  $tsi_string = "'".join("', '", keys(%tsi_del))."'";
  $tlsi_string = "'".join("', '", keys(%tlsi_del))."'";

  $support->log_stamped("Done loading ".scalar(keys(%gsi_del))." gene, ".scalar(keys(%tsi_del))." transcript and ".scalar(keys(%tlsi_del))." translation stable IDs.\n\n");

}


=head2 determine_deleted

  Example     : &determine_deleted;
  Description : Infer deleted genes/transcripts from dbs by comparing which
                objects are in the old and new db.
  Return type : none
  Exceptions  : none
  Caller      : main()
  Status      : At Risk
              : under development

=cut

sub determine_deleted {
  
  $support->log_stamped("Determining list of deleted gene, transcript and translation stable IDs by comparing dbs...\n");

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

  if ($support->param('skip_ncrna')) {

    %skip_biotypes = map { $_ => 1 } @nc_biotypes;

    # make sure we have a complete list of ncRNA biotypes
    my $sql = qq(SELECT DISTINCT biotype from gene);
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

  # get old and new genes and transcripts from db
  my $ga_old = $dba_old->get_GeneAdaptor;
  my @genes_old = @{ $ga_old->fetch_all };

  my $ga_new = $dba_new->get_GeneAdaptor;
  my %genes_new = map { $_->stable_id => $_ } @{ $ga_new->fetch_all };
  my $ta_new = $dba_new->get_TranscriptAdaptor;
  my %tsi_new = map { $_ => 1 } @{ $ta_new->list_stable_ids };

  while (my $g_old = shift(@genes_old)) {
    
    # skip biotypes
    next if ($skip_biotypes{$g_old->biotype});
  
    my $gsi = $g_old->stable_id;
  
    # mark gene as deleted
    unless ($genes_new{$gsi}) {
      $gsi_del{$gsi} = 1;
      $genes{$gsi} = $g_old;
    }

    # transcripts
    foreach my $tr (@{ $g_old->get_all_Transcripts }) {

      # skip biotypes
      next if ($skip_biotypes{$tr->biotype});
  
      my $tsi = $tr->stable_id;

      unless ($tsi_new{$tsi}) {
        # mark transcript and translation as deleted
        $tsi_del{$tsi} = 1;
        my $tl = $tr->translation;
        $tlsi_del{$tl->stable_id} = 1 if ($tl);

        # mark gene as modified
        $gsi_mod{$gsi} = 1;
        $genes_mod{$gsi} = $g_old;
        $genes{$gsi} = $g_old;
      }
    }
  }

  # create stable ID strings for use in mysql IN statements
  $gsi_string = "'".join("', '", keys(%gsi_del))."'";
  $gsi_mod_string = "'".join("', '", keys(%gsi_mod))."'";
  $tsi_string = "'".join("', '", keys(%tsi_del))."'";
  $tlsi_string = "'".join("', '", keys(%tlsi_del))."'";

  # stats
  my $fmt = "%-15s%6d\n";
  $support->log("Deleted objects found:\n", 1);
  $support->log(sprintf($fmt, "genes", scalar(keys(%gsi_del))), 2);
  $support->log(sprintf($fmt, "transcripts", scalar(keys(%tsi_del))), 2);
  $support->log(sprintf($fmt, "translations", scalar(keys(%tlsi_del))), 2);
  
  $support->log("Modified genes found: ".scalar(keys(%gsi_mod))."\n", 1);

  $support->log_stamped("Done.\n\n");
}


=head2 create_mapping_session

  Example     : my $msi = &create_mapping_session;
  Description : Creates a new mapping_session in the db.
  Return type : Int - mapping_session_id of the newly created entry
  Exceptions  : none
  Caller      : main()
  Status      : At Risk
              : under development

=cut

sub create_mapping_session {

  $support->log("Creating new mapping session...\n");
  
  my $old_db_name = $support->param('altdbname');
  my $new_db_name = $support->param('dbname');
  
  my $old_mca = $dba_old->get_MetaContainer;
  my ($old_release) = @{ $old_mca->list_value_by_key('schema_version') };
  my ($old_assembly) = @{ $old_mca->list_value_by_key('assembly.default') };

  my $new_mca = $dba_new->get_MetaContainer;
  my ($new_release) = @{ $new_mca->list_value_by_key('schema_version') };
  my ($new_assembly) = @{ $new_mca->list_value_by_key('assembly.default') };
  
  my $sql = qq(
    INSERT INTO mapping_session
    VALUES (NULL, '$old_db_name', '$new_db_name', '$old_release',
            '$new_release','$old_assembly', '$new_assembly', NOW())
  );
  my $c = $dbh_new->do($sql) unless ($support->param('dry_run'));
  my $mapping_session_id = $dbh_new->{'mysql_insertid'};
  
  my $fmt = "%-23s%-40s\n";
  $support->log(sprintf($fmt, 'mapping_session_id', $mapping_session_id), 1);
  $support->log(sprintf($fmt, 'old_db_name', $old_db_name), 1);
  $support->log(sprintf($fmt, 'new_db_name', $new_db_name), 1);
  $support->log(sprintf($fmt, 'old_release', $old_release), 1);
  $support->log(sprintf($fmt, 'new_release', $new_release), 1);
  $support->log(sprintf($fmt, 'old_assembly', $old_assembly), 1);
  $support->log(sprintf($fmt, 'new_assembly', $new_assembly), 1);

  $support->log("Done.\n\n");

  return $mapping_session_id;
}


=head2 create_stable_id_events 

  Example     : &create_stable_id_events
  Description : Creates stable_id_event entries for all objects found in the old
                db, mapping them to themselves. Optionally, some biotypes will
                be skipped (this is useful if a separate script is run to deal
                with ncRNAs).
  Return type : none
  Exceptions  : none
  Caller      : main()
  Status      : At Risk
              : under development

=cut

sub create_stable_id_events {

  $support->log_stamped("Creating stable_id_event entries for all objects, mapping to themselves...\n");

  my $ga = $dba_old->get_GeneAdaptor;
  my @genes = @{ $ga->fetch_all };

  my $sql = qq(
    INSERT INTO stable_id_event 
    VALUES (?, ?, ?, ?, ?, ?, ?)
  );
  my $sth = $dbh_new->prepare($sql);

  my %stats = map { $_ => 0 } qw(g tr tl g_tot tr_tot);
  my $num_genes = scalar(@genes);
  my $i;

  while (my $gene = shift(@genes)) {

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


=head2 increment_gene_versions 

  Example     : &increment_gene_versions;
  Description : Increment version of all genes where transcripts were deleted.
                Also checks that gene_stable_id.stable_id is correct (and adjusts
                if necessary).
  Return type : none
  Exceptions  : none
  Caller      : main()
  Status      : At Risk
              : under development

=cut

sub increment_gene_versions {
  
  $support->log_stamped("Incrementing gene versions for genes where transcripts were deleted...\n");

  # update stable_id_event
  my $sql = qq(
    UPDATE stable_id_event
    SET new_version = new_version + 1
    WHERE new_stable_id IN ($gsi_mod_string)
    AND mapping_session_id = $mapping_session_id
  );
  my $c = 0;
  $c = $dbh_new->do($sql) unless ($support->param('dry_run'));
  $support->log("stable_id_event [$c]\n", 1);

  # update gene_stable_id
  $sql = qq(
    UPDATE gene_stable_id
    SET version = ?
    WHERE stable_id = ?
    AND version < ?
  );
  my $sth = $dbh_new->prepare($sql);

  $c = 0;

  foreach my $g (values(%genes_mod)) {
    my $version = $g->version + 1;
    $c += $sth->execute($version, $g->stable_id, $version)
      unless ($support->param('dry_run'));
  }
  
  $support->log("gene_stable_id [$c]\n", 1);

  $support->log_stamped("Done.\n\n");
}


=head2 mark_deleted

  Example     : &mark_deleted;
  Description : Sets stable_id_event.new_stable_id to NULL for deleted objects.
  Return type : none
  Exceptions  : none
  Caller      : main()
  Status      : At Risk
              : under development

=cut

sub mark_deleted {

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


=head2 populate_archive

  Example     : &populate_archive;
  Description : Populates gene_archive and peptide_archive for all deleted
                transcripts.
  Return type : none
  Exceptions  : none
  Caller      : main()
  Status      : At Risk
              : under development

=cut

sub populated_archive {

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
      next unless ($tsi_del{$trans->stable_id});
    
      my $tl = $trans->translation;

      # add peptide_archive entry
      unless ($support->param('dry_run')) {

        my $tl_stable_id = "";
        my $tl_version = 0;
        my $pid = 0;

        if ($tl) {
          $tl_stable_id = $tl->stable_id;
          $tl_version = $tl->version;
          my $pep_seq = $trans->translate->seq;
          $sth_pep->execute(md5_hex($pep_seq), $pep_seq);
          $pid = $dbh_new->{'mysql_insertid'};
        }

        # add gene_archive entry
        $sth_gene->execute(
            $gene->stable_id,
            $gene->version,
            $trans->stable_id,
            $trans->version,
            $tl_stable_id,
            $tl_version,
            $pid,
            $mapping_session_id
        );
      }

      $c++;
    }
  }

  $support->log_stamped("Done adding $c entries.\n\n");
}

