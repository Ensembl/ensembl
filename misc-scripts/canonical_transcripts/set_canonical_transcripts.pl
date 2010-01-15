=pod

SYNOPSIS

  Script to set canonical transcripts for each gene

DESCRIPTION

  The rule is if the cluster contains translated transcripts take the
  one with the longest cds otherwise take the one with the longest
  cdna

EXAMPLE

  perl set_canonical_transcripts.pl -dbhost host -dbuser user -dbpass *** -dbname dbname -dbport 3306 -coord_system toplevel -write

EXTRA

  To check the script has run correctly you can run the
  CoreForeignKeys healthcheck:

  ./run-healthcheck.sh -d dbname -output problem CoreForeignKeys

=cut

#!/usr/local/ensembl/bin/perl

use strict;
use Carp;
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
    print "Looking at gene: ", $gene->dbID, "\n";
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
        unless ( $trans{$key}->[0]->translation->seq =~ /\*/ ) {
          $canonical{ $gene->dbID } = $trans_id;
          next GENE;
        } else {
          carp(   "Transcript $trans_id has internal stop(s) "
                . "and it shouldn't if it is protein coding\n" );
        }
      }
    }

    my $has_translation = 0;
    my $count           = 0;
    my @with_translation;
    my @no_translation;

    foreach my $transcript ( @{$transcripts} ) {

      if (    ( $gene->biotype ne 'pseudogene' )
           && ( $gene->biotype       ne 'processed_transcript' )
           && ( $transcript->biotype ne 'nonsense_mediated_decay' )
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
      my @len_and_trans;
      foreach my $trans (@with_translation) {
        my $h = { trans => $trans, len => $trans->translate->length };
        push @len_and_trans, $h;
      }
      my @tmp_sorted = sort { $b->{len} <=> $a->{len} } @len_and_trans;

      foreach my $h (@tmp_sorted) {
        #print "Adding to sorted " . $h->{trans}->dbID . "\n";
        push @sorted, $h->{trans};
      }
    } else {
      @sorted = sort { $b->length <=> $a->length } @no_translation;
    }
    $canonical{ $gene->dbID } = $sorted[0]->dbID;
  } ## end foreach my $gene (@$genes)

  foreach my $gene_id ( keys(%canonical) ) {
    my $transcript_id = $canonical{$gene_id};
    print "Updating gene $gene_id\n";
    $gene_sth->execute( $transcript_id, $gene_id ) if ($write);
  }
} ## end foreach my $slice (@$slices)

