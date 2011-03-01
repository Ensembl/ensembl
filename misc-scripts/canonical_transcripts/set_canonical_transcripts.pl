=pod

=head1 SYNOPSIS

  Script to set canonical transcripts for each gene

=head1 DESCRIPTION

  The rule is if the cluster contains translated transcripts take the
  one with the longest cds otherwise take the one with the longest
  cdna

=head1 OPTIONS

=head2 Required DB connection arguments

  -dbname       Databsse name

  -dbhost       Database host

  -dbport       Database port

  -dbuser       Database user

  -dbpass       Database password

=head2  Optional DB connection arguments

  -dnadbname    DNA Databsse name

  -dnadbhost    DNA Database host

  -dnadbuser    DNA Database user

  -ccds_dbname  CCDS database name

  -ccds_host    CCDS database host

  -ccds_user    CCDS database user


=head2 Other optional arguments

  -coord_system_name    Coordinate system to use

  -include_non_ref      Specify if the non_reference regions should be _excluded_. (default: include) 

  -seq_region_name      Chromosome name if running a single seq_region

  -write                Specify if results should be written to the database

  -verbose              Increase verbosity of output messages

=head1 EXAMPLE

  perl set_canonical_transcripts.pl -dbhost host -dbuser user -dbpass *** -dbname dbname -dbport 3306 -coord_system toplevel -write

=head1 EXTRA

  To check the script has run correctly you can run the
  CoreForeignKeys healthcheck:

  ./run-healthcheck.sh -d dbname -output problem CoreForeignKeys

=cut

#!/usr/local/ensembl/bin/perl

use strict;
use warnings;
use Carp;
use Data::Dumper;
use DBI qw(:sql_types);
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw);

my ($host, $port, $dbname, $user,$pass);
my ($dnahost, $dnadbname, $dnauser);
my ($ccds_host, $ccds_dbname, $ccds_user);

my $coord_system;
my $seq_region_name;
my $write = 0;
my $include_non_ref = 1;
my $verbose = 0;

GetOptions( 'dbhost:s'            => \$host,
            'dbport:n'            => \$port,
            'dbname:s'            => \$dbname,
            'dbuser:s'            => \$user,
            'dbpass:s'            => \$pass,
            'dnahost:s'           => \$dnahost,
            'dnadbname:s'         => \$dnadbname,
            'dnauser:s'           => \$dnauser,
            'ccds_host:s'         => \$ccds_host,
            'ccds_dbname:s'       => \$ccds_dbname,
            'ccds_user:s'         => \$ccds_user,
            'coord_system_name:s' => \$coord_system,
            'seq_region_name:s'   => \$seq_region_name,
            'write!'              => \$write,
            'include_non_ref!'    => \$include_non_ref,
            'verbose!'            => \$verbose, );

unless ($write) {
  print "You have not used the -write option "
      . "so results will not be written into the database\n";
}

my $db =
  new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $host,
                                      -user   => $user,
                                      -port   => $port,
                                      -dbname => $dbname,
                                      -pass   => $pass );
my $ccds_db;

if ($ccds_dbname) {
  if (!$ccds_user || !$ccds_host) {
    throw ("You must provide user, host and dbname details to connect to CCDS DB!");
  }
  $ccds_db = 
  new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $ccds_host,
                                      -user   => $ccds_user,
                                      -port   => $port,
                                      -dbname => $ccds_dbname );
}

my $geneDB_has_DNA = check_if_DB_contains_DNA($db);

my $dnadb;

if ($geneDB_has_DNA == 0 && !$dnadbname) {
  throw ("Your gene DB contains no DNA. You must provide a DNA_DB to connect to");
} elsif ($geneDB_has_DNA == 0 && $dnadbname) {
  $dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
                                                 '-host'   => $dnahost,
                                                 '-user'   => $dnauser,
                                                 '-dbname' => $dnadbname,
                                                 '-port'   => $port,
                                             );
  $db->dnadb($dnadb);
}

my $sa = $db->get_SliceAdaptor;
my $dea = $db->get_DBEntryAdaptor;
my $slices;

# Only update the canonical transcripts in this region.
if ($seq_region_name) {
  print "Only updating genes on chromosome $seq_region_name\n";
  my $slice =
    $sa->fetch_by_region( $coord_system, $seq_region_name, $include_non_ref );
  push( @$slices, $slice );
} else {
  $slices = $sa->fetch_all( $coord_system, '', $include_non_ref );
}

my $gene_update_sql = "update gene set canonical_transcript_id = ? where gene_id = ?";
my $gene_sth = $db->dbc->prepare($gene_update_sql);

SLICE: foreach my $slice (@$slices) {
  print "\nGetting genes for " . $slice->name . "\n" if ($verbose);
  my $genes = $slice->get_all_Genes( undef, undef, 1 );
  my %canonical;

  GENE: foreach my $gene (@$genes) {
    print "\nLooking at gene: (dbID " . $gene->dbID . ")";
    if ($gene->stable_id) {
      print " :: " . $gene->stable_id ;
    }
    print "\n";
    my $transcripts = $gene->get_all_Transcripts;
    if ( @$transcripts == 1 ) {
      if ($ccds_db) {
        check_Ens_trans_against_CCDS($transcripts->[0], $ccds_db, $verbose);
      }
      $canonical{ $gene->dbID } = $transcripts->[0]->dbID;
      print "Transcript ". $transcripts->[0]->display_id . " of biotype " . $transcripts->[0]->biotype .
            " is chosen as the canonical transcript for gene ". $gene->display_id .
            " of biotype ". $gene->biotype. "\n";
      next GENE;
    }

    # Keep track of the transcript biotypes
    my %trans;
    foreach my $transcript ( @{$transcripts} ) {
      push @{ $trans{ $transcript->biotype } }, $transcript;
    }

    # If there is a single protein_coding transcript, make it the
    # canonical transcript, if it has a functional translation.
    foreach my $key ( keys(%trans) ) {
      if ( ( $key eq 'protein_coding' ) && ( @{ $trans{$key} } == 1 ) ) {
        my $trans_id = $trans{$key}->[0]->dbID;
        if ( $trans{$key}->[0]->translation->seq !~ /\*/ ) {
          if ($ccds_db) {
            check_Ens_trans_against_CCDS($trans{$key}->[0], $ccds_db, $verbose);
          }
          $canonical{ $gene->dbID } = $trans_id;
          print "Transcript ". $trans{$key}->[0]->display_id . " of biotype " . $trans{$key}->[0]->biotype .
            " is chosen as the canonical transcript for gene ". $gene->display_id . 
            " of biotype ". $gene->biotype. "\n";
          next GENE;
        } else {
          carp(   "Transcript $trans_id has internal stop(s) "
                . "and shouldn't if it is protein coding\n" );
        }
      } elsif (    !exists( $trans{'protein_coding'} )
              && $key eq 'nonsense_mediated_decay'
              && ( @{ $trans{$key} } == 1 ) ) {
        my $trans_id = $trans{$key}->[0]->dbID;
        if (    ( $trans{$key}->[0]->translation )
             && ( $trans{$key}->[0]->translation->seq !~ /\*/ ) ) {
          print "The NMD transcript $trans_id will become "
              . "the canonical transcript because there are no "
              . "protein coding transcripts.\n";
          if ($ccds_db) {
            check_Ens_trans_against_CCDS($trans{$key}->[0], $ccds_db, $verbose);
          }
          $canonical{ $gene->dbID } = $trans_id;
          print "Transcript ". $trans{$key}->[0]->display_id . " of biotype nonsense_mediated_decay" .
            " is chosen as the canonical transcript for gene ". $gene->display_id .
            " of biotype ". $gene->biotype. "\n";
          next GENE;
        } else {
          carp(   "Transcript $trans_id is an NMD with either "
                . " no translation or internal stop(s), skipping it...\n" );
        }
      }
    }

    my $has_translation = 0;
    my $count           = 0;
    my @with_translation;
    my @no_translation;

    foreach my $transcript ( @{$transcripts} ) {

      if (    ( $gene->biotype       ne 'pseudogene' )
           && ( $gene->biotype       ne 'processed_transcript' )
           && ( $transcript->biotype ne 'polymorphic_pseudogene' )
           && ( $transcript->biotype ne 'processed_transcript' )
           && $transcript->translation
           && ( $transcript->translation->seq !~ /\*/ ) )

      {
        push( @with_translation, $transcript );
      } else {
        push( @no_translation, $transcript );
      }
    }

    my @sorted;

    if (@with_translation) {

      my ( @ccds_prot_coding_len_and_trans, @merge_prot_coding_len_and_trans, @prot_coding_len_and_trans);
      my ( @merge_translateable_len_and_trans, @translateable_len_and_trans );
      my ( @merge_len_and_trans, @len_and_trans );

      foreach my $trans (@with_translation) {

        print "Looking at Ens trans dbID " . $trans->dbID . " :: biotype " . $trans->biotype . " :: Transl length " .  
              $trans->translate->length . "\n" if ($verbose);

        my $h = { trans => $trans, len => $trans->translate->length };
        my $ensembl_trans_found_in_CCDS = 0;

        if ( $trans->biotype eq 'protein_coding') {
          if ($ccds_db) {
            $ensembl_trans_found_in_CCDS = check_Ens_trans_against_CCDS($trans, $ccds_db, $verbose);
          }
          if ($ensembl_trans_found_in_CCDS) {
            push @ccds_prot_coding_len_and_trans, $h;
          } elsif ($trans->analysis->logic_name eq 'ensembl_havana_transcript') {
            push @merge_prot_coding_len_and_trans, $h;
          } else {
            push @prot_coding_len_and_trans, $h;
          }
          next;
        } else { # transcript biotype is not protein_coding but the transcript has translation
          if ($trans->analysis->logic_name eq 'ensembl_havana_transcript') {
            push @merge_translateable_len_and_trans, $h;
          } else {
            push @translateable_len_and_trans, $h;
          }
        }
      }

      my @tmp_sorted;
      if (@ccds_prot_coding_len_and_trans) {
        @tmp_sorted = sort { $b->{len} <=> $a->{len} } @ccds_prot_coding_len_and_trans;  # sort all the Ens transcripts which matched CCDS by coding length
      } elsif (@merge_prot_coding_len_and_trans) {
        @tmp_sorted = sort { $b->{len} <=> $a->{len} } @merge_prot_coding_len_and_trans;
      } elsif (@prot_coding_len_and_trans) {
        @tmp_sorted = sort { $b->{len} <=> $a->{len} } @prot_coding_len_and_trans;
      } elsif (@merge_len_and_trans) {
        print "There are no 'protein_coding' labelled transcripts, ".
              "will take the longest of merged tranlslateable transcripts\n" if ($verbose);
        @tmp_sorted = sort { $b->{len} <=> $a->{len} } @merge_translateable_len_and_trans;
      } else {
        print "There are neither 'protein_coding' labelled transcripts nor merged ".
              "translateable transcripts. Will take the longest of any unmerged ".
              "translateable transcript.\n" if ($verbose);
        @tmp_sorted = sort { $b->{len} <=> $a->{len} } @translateable_len_and_trans;
      }

      foreach my $h (@tmp_sorted) {
        print "Adding to sorted " . $h->{trans}->dbID . " :: " . $h->{trans}->biotype . "\n" if ($verbose);
        push @sorted, $h->{trans};
      }
    } else {
      # we merge genes without a translation too
      my ( @merge_len_and_trans, @len_and_trans );
      foreach my $trans (@no_translation) {
        if ($trans->analysis->logic_name eq 'ensembl_havana_transcript') {
          push @merge_len_and_trans, $trans;
        } else {
          push @len_and_trans, $trans;
        }
      }
      if (@merge_len_and_trans) {
        @sorted = sort { $b->length <=> $a->length } @merge_len_and_trans;
      } else {
        @sorted = sort { $b->length <=> $a->length } @len_and_trans;
      }
    } # end of loop for no translation

    # # #
    # set canonical transcirpt
    # # #
    print "Transcript ". $sorted[0]->display_id . " of biotype " . $sorted[0]->biotype .
          " is chosen as the canonical transcript for gene ". $gene->display_id . 
          " of biotype ". $gene->biotype. "\n";
    $canonical{ $gene->dbID } = $sorted[0]->dbID;

  } ## end foreach my $gene (@$genes)

  foreach my $gene_id ( keys(%canonical) ) {
    my $transcript_id = $canonical{$gene_id};
    print "Updating gene $gene_id "
        . "with canonical transcript: $transcript_id.\n";
    $gene_sth->execute( $transcript_id, $gene_id ) if ($write);
  }
} ## end foreach my $slice (@$slices)


sub check_if_DB_contains_DNA {
  my ($db) = @_;
  my $sql_command = "select count(*) from dna";
  my $sth = $db->dbc->prepare($sql_command);
  $sth->execute();
  my @dna_array = $sth->fetchrow_array;
  if ($dna_array[0] > 0) {
    print "Your DB ". $db->dbc->dbname ." contains DNA sequences. No need to attach a ".
          "DNA_DB to it.\n" if ($verbose);
    return 1;
  } else {
    print "Your DB ". $db->dbc->dbname ." does not contain DNA sequences.\n"
          if ($verbose);
    return 0;
  }
}

sub check_Ens_trans_against_CCDS {
  my ($ensembl_trans, $ccds_db, $verbose) = @_;
 
  # check if the ensembl transcript is on 
  # a reference or non-reference slice
  # CCDS is only on reference slices
  if (!$ensembl_trans->slice->is_reference) {
    print "Transcript dbID ".$ensembl_trans->dbID." is on non-reference slice ".$ensembl_trans->slice->name.", therefore no CCDS\n" if ($verbose);
    return 0;
  } else {
    print "Transcript dbID ".$ensembl_trans->dbID." is on reference slice ".$ensembl_trans->slice->name."\n" if ($verbose);
  }

  my @ensembl_translateable_exons = @{ $ensembl_trans->get_all_translateable_Exons };

  my $ext_slice = $ccds_db->get_SliceAdaptor->fetch_by_region( 'toplevel',
                     $ensembl_trans->slice->seq_region_name,
                     $ensembl_trans->seq_region_start, $ensembl_trans->seq_region_end );

  EXT_GENE: foreach my $ext_gene ( @{ $ext_slice->get_all_Genes } ) {
    EXT_TRANS: foreach my $ext_trans ( @{ $ext_gene->get_all_Transcripts } ) {
      my @ext_exons = @{ $ext_trans->get_all_Exons };

      if ( scalar(@ensembl_translateable_exons) == scalar(@ext_exons) ) {
        for ( my $i = 0 ; $i < scalar(@ensembl_translateable_exons) ; $i++ ) {
          if ( $ensembl_translateable_exons[$i]->coding_region_start($ensembl_trans) != $ext_exons[$i]->seq_region_start
             || $ensembl_translateable_exons[$i]->strand != $ext_exons[$i]->strand
             || $ensembl_translateable_exons[$i]->coding_region_end($ensembl_trans) != $ext_exons[$i]->seq_region_end )
          {
             next EXT_TRANS;
          }
        }
        print "Ensembl transcript " . $ensembl_trans->display_id . " found match " . 
          $ext_gene->display_id . " in CCDS DB.\n";
        return 1;
      }
    } # end foreach EXT_TRANS
  } # end foreach EXT_GENE 
} ## end sub

