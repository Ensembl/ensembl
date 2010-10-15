=pod

=head1 SYNOPSIS

  Script to set canonical transcripts for each gene

=head1 DESCRIPTION

  The rule is if the cluster contains translated transcripts take the
  one with the longest cds otherwise take the one with the longest
  cdna

=head1 OPTIONS

=head2 Required arguments

  -dbname       Databsse name

  -dbhost       Database host

  -dbport       Database port

  -dbuser       Database user

  -dbpass       Database password

=head2  Optional arguments

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

my $host;
my $port;
my $dbname;
my $user;
my $pass;
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
  print "Getting genes for " . $slice->name . "\n" if ($verbose);
  my $genes = $slice->get_all_Genes( undef, undef, 1 );
  my %canonical;

GENE: foreach my $gene (@$genes) {
    print "Looking at gene: (" . $gene->dbID . ") :: " . $gene->stable_id . "\n";
    my $transcripts = $gene->get_all_Transcripts;
    if ( @$transcripts == 1 ) {
      $canonical{ $gene->dbID } = $transcripts->[0]->dbID;
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
          $canonical{ $gene->dbID } = $trans_id;
          next GENE;
        } else {
          carp(   "Transcript $trans_id has internal stop(s) "
                . "and shouldn't if it is protein coding\n" );
        }
      }
      elsif (    !exists( $trans{'protein_coding'} )
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

          $canonical{ $gene->dbID } = $trans_id;
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
      my ( @ccds_prot_coding_len_and_trans, @merge_prot_coding_len_and_trans, @prot_coding_len_and_trans, @merge_len_and_trans, @len_and_trans );
      foreach my $trans (@with_translation) {

        if ($verbose) {
          print $trans->dbID . " :: " . $trans->biotype . " :: " .  $trans->translate->length . "\n";
        }

        if ( $trans->biotype eq 'protein_coding') {
          my $h = { trans => $trans, len => $trans->translate->length };
          if (scalar(@{$dea->fetch_all_by_Transcript($trans, 'CCDS')})) {
            push @ccds_prot_coding_len_and_trans, $h;
          } elsif ($trans->analysis->logic_name eq 'ensembl_havana_transcript') {
            push @merge_prot_coding_len_and_trans, $h;
          } else {
            push @prot_coding_len_and_trans, $h;
          }
          next;
        } else {
          my $h = { trans => $trans, len => $trans->translate->length };
          if ($trans->analysis->logic_name eq 'ensembl_havana_transcript') {
            push @merge_len_and_trans, $h;
          } else {
            push @len_and_trans, $h;
          }
        }
      }

      my @tmp_sorted;
      if (@ccds_prot_coding_len_and_trans) {
        @tmp_sorted = sort { $b->{len} <=> $a->{len} } @ccds_prot_coding_len_and_trans;
      } elsif (@merge_prot_coding_len_and_trans) {
        @tmp_sorted = sort { $b->{len} <=> $a->{len} } @merge_prot_coding_len_and_trans;
      } elsif (@prot_coding_len_and_trans) {
        @tmp_sorted = sort { $b->{len} <=> $a->{len} } @prot_coding_len_and_trans;
      } elsif (@merge_len_and_trans) {
        @tmp_sorted = sort { $b->{len} <=> $a->{len} } @merge_len_and_trans;
      } else {
        if ($verbose) {
          print "There are no 'protein_coding' labelled transcripts, "
              . "will take the longest of the other coding transcripts\n";
        }
        @tmp_sorted = sort { $b->{len} <=> $a->{len} } @len_and_trans;
      }

      foreach my $h (@tmp_sorted) {

        if ($verbose) {
          print "Adding to sorted " . $h->{trans}->dbID . " :: " . $h->{trans}->biotype . "\n";
        }
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
    $canonical{ $gene->dbID } = $sorted[0]->dbID;

  } ## end foreach my $gene (@$genes)

  foreach my $gene_id ( keys(%canonical) ) {
    my $transcript_id = $canonical{$gene_id};
    print "Updating gene $gene_id "
        . "with canonical transcript: $transcript_id.\n";
    $gene_sth->execute( $transcript_id, $gene_id ) if ($write);
  }
} ## end foreach my $slice (@$slices)

