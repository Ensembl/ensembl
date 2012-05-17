
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

  -include_duplicates    Specify if the duplicate regions should be _included_. eg. Human PAR on Y (default: exclude) 

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

#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Data::Dumper;
use DBI qw(:sql_types);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Utils::CliHelper;

my $cli_helper = Bio::EnsEMBL::Utils::CliHelper->new();

# get the basic options for connecting to a database server
my $optsd = [ @{ $cli_helper->get_dba_opts() },
              @{ $cli_helper->get_dba_opts('dna') },
              @{ $cli_helper->get_dba_opts('ccds') } ];
# add the print option
push( @{$optsd}, "coord_system_name:s" );
push( @{$optsd}, "logic_name:s" );
push( @{$optsd}, "write!" );
push( @{$optsd}, "include_non_ref!" );
push( @{$optsd}, "include_duplicates" );
push( @{$optsd}, "verbose!" );
# process the command line with the supplied options plus a help subroutine
my $opts = $cli_helper->process_args( $optsd, \&usage );

$opts->{write}           ||= 0;
$opts->{include_non_ref} ||= 1;
$opts->{verbose}         ||= 0;

unless ( $opts->{write} ) {
  print "You have not used the -write option "
    . "so results will not be written into the database\n";
}

my $db_args = $cli_helper->get_dba_args_for_opts($opts);
my $dnadb_args = $cli_helper->get_dba_args_for_opts( $opts, 'dna' );
if ( defined $dnadb_args && scalar(@$dnadb_args) != scalar(@$db_args) ) {
  croak "Different number of DBAs found for DB and DNADB";
}
my $ccdsdb_args = $cli_helper->get_dba_args_for_opts( $opts, 'ccds' );
if ( defined $ccdsdb_args && scalar(@$ccdsdb_args) != scalar(@$db_args) ) {
  croak "Different number of DBAs found for DB and CCDSDB";
}

for ( my $i = 0 ; $i < scalar(@$db_args) ; $i++ ) {

  my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new( %{ $db_args->[$i] } );

  my $geneDB_has_DNA = check_if_DB_contains_DNA($db);

  my $dnadb;

  if ( $geneDB_has_DNA == 0 && !$dnadb_args ) {
    throw(
       "Your gene DB contains no DNA. You must provide a DNA_DB to connect to"
    );
  } elsif ( $geneDB_has_DNA == 0 && $dnadb_args ) {
    $dnadb = Bio::EnsEMBL::DBSQL::DBAdaptor->new( %{ $dnadb_args->[$i] } );
    $db->dnadb($dnadb);
  }

  my $ccds_db;
  if ( defined $ccdsdb_args ) {
    $ccds_db = Bio::EnsEMBL::DBSQL::DBAdaptor->new( %{ $ccdsdb_args->[$i] } );
  }

  my $sa  = $db->get_SliceAdaptor;
  my $dea = $db->get_DBEntryAdaptor;
  my $slices;

  # Only update the canonical transcripts in this region.
  if ( $opts->{seq_region_name} ) {
    print "Only updating genes on chromosome "
      . $opts->{seq_region_name} . "\n";
    my $slice =
      $sa->fetch_by_region( $opts->{coord_system_name},
                            $opts->{seq_region_name},
                            $opts->{include_non_ref},
                            $opts->{include_duplicates} );
    push( @$slices, $slice );
  } else {
    $slices = $sa->fetch_all( $opts->{coord_system_name},
                              '',
                              $opts->{include_non_ref},
                              $opts->{include_duplicates} );
  }

  my $gene_update_sql =
    "update gene set canonical_transcript_id = ? where gene_id = ?";
  my $gene_sth = $db->dbc->prepare($gene_update_sql);

SLICE: foreach my $slice (@$slices) {
    # check whether the slice is reference or not
    # this check is important for species that have a ccds database
    # and that have non-reference sequence ie. human
    my $slice_is_reference;
    if ( !$slice->is_reference ) {
      print "Doing non-reference slice " . $slice->name . "\n";
      $slice_is_reference = 0;
    } else {
      print "Doing reference slice " . $slice->name . "\n";
      $slice_is_reference = 1;
    }

    # Now fetch the genes
    print "\nGetting genes for " . $slice->name . "\n"
      if ( $opts->{verbose} );
    my $genes = $slice->get_all_Genes( $opts->{logic_name}, undef, 1 );
    my %canonical;

  GENE: foreach my $gene (@$genes) {
      print "\nLooking at gene: (dbID " . $gene->dbID . ")";
      if ( $gene->stable_id ) {
        print " :: " . $gene->stable_id;
      }
      print "\n";
      my $transcripts = $gene->get_all_Transcripts;
      if ( @$transcripts == 1 ) {
        if ( $ccds_db && $slice_is_reference ) {
          check_Ens_trans_against_CCDS( $transcripts->[0], $ccds_db,
                                        $opts->{verbose} );
        }
        $canonical{ $gene->dbID } = $transcripts->[0]->dbID;
        print "Transcript "
          . $transcripts->[0]->display_id
          . " of biotype "
          . $transcripts->[0]->biotype
          . " is chosen as the canonical transcript for gene "
          . $gene->display_id
          . " of biotype "
          . $gene->biotype . "\n";
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
        if (    ( $key eq 'protein_coding' )
             && ( @{ $trans{$key} } == 1 ) )
        {
          my $trans_id = $trans{$key}->[0]->dbID;
          if ( $trans{$key}->[0]->translation->seq !~ /\*/ ) {
            if ( $ccds_db && $slice_is_reference ) {
              check_Ens_trans_against_CCDS( $trans{$key}->[0],
                                            $ccds_db, $opts->{verbose} );
            }
            $canonical{ $gene->dbID } = $trans_id;
            print "Transcript "
              . $trans{$key}->[0]->display_id
              . " of biotype "
              . $trans{$key}->[0]->biotype
              . " is chosen as the canonical transcript for gene "
              . $gene->display_id
              . " of biotype "
              . $gene->biotype . "\n";
            next GENE;
          } else {
            carp(   "Transcript $trans_id has internal stop(s) "
                  . "and shouldn't if it is protein coding\n" );
          }
        } elsif (    !exists( $trans{'protein_coding'} )
                  && $key eq 'nonsense_mediated_decay'
                  && ( @{ $trans{$key} } == 1 ) )
        {
          my $trans_id = $trans{$key}->[0]->dbID;
          if (    ( $trans{$key}->[0]->translation )
               && ( $trans{$key}->[0]->translation->seq !~ /\*/ ) )
          {
            print "The NMD transcript $trans_id will become "
              . "the canonical transcript because there are no "
              . "protein coding transcripts.\n";
            if ( $ccds_db && $slice_is_reference ) {
              check_Ens_trans_against_CCDS( $trans{$key}->[0],
                                            $ccds_db, $opts->{verbose} );
            }
            $canonical{ $gene->dbID } = $trans_id;
            print "Transcript "
              . $trans{$key}->[0]->display_id
              . " of biotype nonsense_mediated_decay"
              . " is chosen as the canonical transcript for gene "
              . $gene->display_id
              . " of biotype "
              . $gene->biotype . "\n";
            next GENE;
          } else {
            carp(   "Transcript $trans_id is an NMD with either "
                  . " no translation or internal stop(s), skipping it...\n" );
          }
        } ## end elsif ( !exists( $trans{'protein_coding'...
      }    # foreach biotype

      my $has_translation = 0;
      my $count           = 0;
      my @with_translation;
      my @no_translation;

      foreach my $transcript ( @{$transcripts} ) {

        if (    ( $gene->biotype ne 'pseudogene' )
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

        my ( @ccds_prot_coding_len_and_trans,
             @merge_prot_coding_len_and_trans,
             @prot_coding_len_and_trans );
        my ( @merge_translateable_len_and_trans,
             @translateable_len_and_trans );
        my ( @merge_len_and_trans, @len_and_trans );

        foreach my $trans (@with_translation) {

          print "Looking at Ens trans dbID "
            . $trans->dbID    . " :: biotype "
            . $trans->biotype . " :: Transl length "
            . $trans->translate->length . "\n"
            if ( $opts->{verbose} );

          my $h = { trans => $trans, len => $trans->translate->length };
          my $ensembl_trans_found_in_CCDS = 0;

          if ( $trans->biotype eq 'protein_coding' ) {
            if ( $ccds_db && $slice_is_reference ) {
              $ensembl_trans_found_in_CCDS =
                check_Ens_trans_against_CCDS( $trans, $ccds_db,
                                              $opts->{verbose} );
            }
            if ($ensembl_trans_found_in_CCDS) {
              push @ccds_prot_coding_len_and_trans, $h;
            } elsif (
                 $trans->analysis->logic_name eq 'ensembl_havana_transcript' )
            {
              push @merge_prot_coding_len_and_trans, $h;
            } else {
              push @prot_coding_len_and_trans, $h;
            }
            next;
          } else { # transcript biotype is not protein_coding but the transcript has translation
            if ( $trans->analysis->logic_name eq 'ensembl_havana_transcript' )
            {
              push @merge_translateable_len_and_trans, $h;
            } else {
              push @translateable_len_and_trans, $h;
            }
          }
        } ## end foreach my $trans (@with_translation)

        my @tmp_sorted;
        if (@ccds_prot_coding_len_and_trans) {
          @tmp_sorted =
            sort { $b->{len} <=> $a->{len} }
            @ccds_prot_coding_len_and_trans
            ; # sort all the Ens transcripts which matched CCDS by coding length
        } elsif (@merge_prot_coding_len_and_trans) {
          @tmp_sorted =
            sort { $b->{len} <=> $a->{len} } @merge_prot_coding_len_and_trans;
        } elsif (@prot_coding_len_and_trans) {
          @tmp_sorted =
            sort { $b->{len} <=> $a->{len} } @prot_coding_len_and_trans;
        } elsif (@merge_len_and_trans) {
          print "There are no 'protein_coding' labelled transcripts, "
            . "will take the longest of merged tranlslateable transcripts\n"
            if ( $opts->{verbose} );
          @tmp_sorted =
            sort { $b->{len} <=> $a->{len} }
            @merge_translateable_len_and_trans;
        } else {
          print "There are neither 'protein_coding' labelled transcripts nor merged "
              . "translateable transcripts. Will take the longest of any unmerged "
              . "translateable transcript.\n"
            if ( $opts->{verbose} );
          @tmp_sorted =
            sort { $b->{len} <=> $a->{len} } @translateable_len_and_trans;
        }

        foreach my $h (@tmp_sorted) {
          print "Adding to sorted "
            . $h->{trans}->dbID . " :: "
            . $h->{trans}->biotype . "\n"
            if ( $opts->{verbose} );
          push @sorted, $h->{trans};
        }
      } else {
        # we merge genes without a translation too
        my ( @merge_len_and_trans, @len_and_trans );
        foreach my $trans (@no_translation) {
          if ( $trans->analysis->logic_name eq 'ensembl_havana_transcript' ) {
            push @merge_len_and_trans, $trans;
          } else {
            push @len_and_trans, $trans;
          }
        }
        if (@merge_len_and_trans) {
          @sorted =
            sort { $b->length <=> $a->length } @merge_len_and_trans;
        } else {
          @sorted = sort { $b->length <=> $a->length } @len_and_trans;
        }
      }    # end of loop for no translation

      # # #
      # set canonical transcirpt
      # # #
      print "Transcript "
        . $sorted[0]->display_id
        . " of biotype "
        . $sorted[0]->biotype
        . " is chosen as the canonical transcript for gene "
        . $gene->display_id
        . " of biotype "
        . $gene->biotype . "\n";
      $canonical{ $gene->dbID } = $sorted[0]->dbID;

      if ( !exists $canonical{ $gene->dbID } ) {
        throw(
             "No canonical transcript has been set for gene " . $gene->dbID );
      }
    } ## end foreach my $gene (@$genes)

    foreach my $gene_id ( keys(%canonical) ) {
      my $transcript_id = $canonical{$gene_id};
      print "Updating gene $gene_id "
        . "with canonical transcript: $transcript_id.\n";
      $gene_sth->execute( $transcript_id, $gene_id )
        if ( $opts->{write} );
    }
  } ## end foreach my $slice (@$slices)

} ## end for ( my $i = 0 ; $i < ...

sub check_if_DB_contains_DNA {
  my ($db)        = @_;
  my $sql_command = "select count(*) from dna";
  my $sth         = $db->dbc->prepare($sql_command);
  $sth->execute();
  my @dna_array = $sth->fetchrow_array;
  if ( $dna_array[0] > 0 ) {
    print "Your DB "
      . $db->dbc->dbname
      . " contains DNA sequences. No need to attach a "
      . "DNA_DB to it.\n"
      if ( $opts->{verbose} );
    return 1;
  } else {
    print "Your DB " . $db->dbc->dbname . " does not contain DNA sequences.\n"
      if ( $opts->{verbose} );
    return 0;
  }
}

sub check_Ens_trans_against_CCDS {
  my ( $ensembl_trans, $ccds_db, $verbose ) = @_;

  my @ensembl_translateable_exons =
    @{ $ensembl_trans->get_all_translateable_Exons };

  my $ext_slice =
    $ccds_db->get_SliceAdaptor->fetch_by_region(
                                       'toplevel',
                                       $ensembl_trans->slice->seq_region_name,
                                       $ensembl_trans->seq_region_start,
                                       $ensembl_trans->seq_region_end );

EXT_GENE: foreach my $ext_gene ( @{ $ext_slice->get_all_Genes } ) {
  EXT_TRANS:
    foreach my $ext_trans ( @{ $ext_gene->get_all_Transcripts } ) {
      my @ext_exons = @{ $ext_trans->get_all_Exons };

      if ( scalar(@ensembl_translateable_exons) == scalar(@ext_exons) ) {
        for ( my $i = 0 ; $i < scalar(@ensembl_translateable_exons) ; $i++ ) {
          if ( $ensembl_translateable_exons[$i]
               ->coding_region_start($ensembl_trans) !=
               $ext_exons[$i]->seq_region_start
               || $ensembl_translateable_exons[$i]->strand !=
               $ext_exons[$i]->strand
               || $ensembl_translateable_exons[$i]
               ->coding_region_end($ensembl_trans) !=
               $ext_exons[$i]->seq_region_end )
          {
            next EXT_TRANS;
          }
        }
        print "Ensembl transcript "
          . $ensembl_trans->display_id
          . " found match "
          . $ext_gene->display_id
          . " in CCDS DB.\n";
        return 1;
      }
    }    # end foreach EXT_TRANS
  }    # end foreach EXT_GENE
} ## end sub check_Ens_trans_against_CCDS

