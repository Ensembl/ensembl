=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

package XrefMapper::CoordinateMapper;

# $Id$

# This is a set of subroutines used for creating Xrefs based on
# coordinate overlaps.


use strict;
use warnings;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Mapper::RangeRegistry;

use DBI qw( :sql_types );

use Carp;
use IO::File;
use File::Spec::Functions;


use base qw( Exporter XrefMapper::BasicMapper);

our @EXPORT = qw( run_coordinatemapping );

our $coding_weight = 2;
our $ens_weight    = 3;

our $transcript_score_threshold = 0.75;

my $verbose = 0;

sub new {
  my($class, $mapper) = @_;

  my $self ={};
  bless $self,$class;
  $self->core($mapper->core);
  $self->xref($mapper->xref);
  $self->mapper($mapper);
  $verbose = $mapper->verbose();
  return $self;
}

sub mapper{
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_mapper} = $arg );
  return $self->{_mapper};
}


sub run_coordinatemapping {
  my ( $self, $do_upload, $species_id) = @_;


  my $sth_stat = $self->xref->dbc->prepare(
                           "INSERT INTO process_status (status, date) "
                         . "VALUES('coordinate_xrefs_started',now())" );
  $sth_stat->execute();
  $sth_stat->finish();

  my $xref_db = $self->xref();
  my $core_db = $self->core();

  my $species = $core_db->species();
#  my $species_id = $self->mapper->core->species;
  $species_id = XrefMapper::BasicMapper::get_species_id_from_species_name( $xref_db, $species ) unless defined $species_id;
  if (!defined $species_id) { return; }


  # We only do coordinate mapping for mouse and human for now.
  if ( !( $species eq 'mus_musculus' || $species eq 'homo_sapiens' ) ) {
    #  if ( !( $species eq 'mus_musculus') ) {
    my $sth_stat = $self->xref->dbc->prepare(
                          "INSERT INTO process_status (status, date) "
                        . "VALUES ('coordinate_xref_finished',now())" );
    $sth_stat->execute();
    $sth_stat->finish();
    return;
  }
  print "species = $species ($species_id))\n";

  my $output_dir = $core_db->dir();

  my $xref_filename = catfile( $output_dir, 'xref_coord.txt' );
  my $object_xref_filename =
    catfile( $output_dir, 'object_xref_coord.txt' );
  my $unmapped_reason_filename =
    catfile( $output_dir, 'unmapped_reason_coord.txt' );
  my $unmapped_object_filename =
    catfile( $output_dir, 'unmapped_object_coord.txt' );

  my $xref_dbh = $xref_db->dbc()->db_handle();
  my $core_dbh = $core_db->dbc()->db_handle();

  ######################################################################
  # Figure out the last used 'xref_id', 'object_xref_id',              #
  # 'unmapped_object_id', and 'unmapped_reason_id' from the Core       #
  # database.                                                          #
  ######################################################################

  my $xref_id =
    $core_dbh->selectall_arrayref('SELECT MAX(xref_id) FROM xref')
    ->[0][0];
  my $object_xref_id = $core_dbh->selectall_arrayref(
                 'SELECT MAX(object_xref_id) FROM object_xref')->[0][0];
  my $unmapped_object_id = $core_dbh->selectall_arrayref(
         'SELECT MAX(unmapped_object_id) FROM unmapped_object')->[0][0];
  my $unmapped_reason_id = $core_dbh->selectall_arrayref(
         'SELECT MAX(unmapped_reason_id) FROM unmapped_reason')->[0][0];

  log_progress( "Last used xref_id            is %d\n", $xref_id );
  log_progress( "Last used object_xref_id     is %d\n",
                $object_xref_id );
  log_progress( "Last used unmapped_object_id is %d\n",
                $unmapped_object_id );
  log_progress( "Last used unmapped_reason_id is %d\n",
                $unmapped_reason_id );

  ######################################################################
  # Get an 'analysis_id', or discover that we need to add our analysis #
  # to the 'analyis' table later.                                      #
  ######################################################################

  my $analysis_params =
    sprintf( "weights(coding,ensembl)="
               . "%.2f,%.2f;"
               . "transcript_score_threshold=" . "%.2f",
             $coding_weight, $ens_weight, $transcript_score_threshold );

  my $analysis_sql = qq(
    SELECT  analysis_id
    FROM    analysis
    WHERE   logic_name = 'xrefcoordinatemapping'
    AND     parameters = ?
  );

  my $analysis_sth = $core_dbh->prepare($analysis_sql);
  $analysis_sth->bind_param( 1, $analysis_params, SQL_VARCHAR );

  $analysis_sth->execute();

  my $analysis_id = $analysis_sth->fetchall_arrayref()->[0][0];
  if ( !defined($analysis_id) ) {
    $analysis_id =
      $core_dbh->selectall_arrayref( "SELECT analysis_id FROM analysis "
               . "WHERE logic_name = 'xrefcoordinatemapping'" )->[0][0];

    if ( defined($analysis_id) && $do_upload ) {
      log_progress(   "Will update 'analysis' table "
                    . "with new parameter settings\n" );

      #-----------------------------------------------------------------
      # Update an existing analysis.
      #-----------------------------------------------------------------

      my $sql = qq(
        UPDATE  analysis
        SET     created = now(), parameters = ?
        WHERE   analysis_id = ?
      );

      $core_dbh->do( $sql, undef, $analysis_params, $analysis_id );

    } else {
      log_progress("Can not find analysis ID for this analysis:\n");
      log_progress("  logic_name = 'xrefcoordinatemapping'\n");
      log_progress( "  parameters = '%s'\n", $analysis_params );

      if ($do_upload) {

        #---------------------------------------------------------------
        # Store a new analysis.
        #---------------------------------------------------------------

        log_progress("A new analysis will be added\n");

        $analysis_id = $core_dbh->selectall_arrayref(
                       'SELECT MAX(analysis_id) FROM analysis')->[0][0];
        log_progress( "Last used analysis_id is %d\n", $analysis_id );

        my $sql =
            'INSERT INTO analysis '
          . 'VALUES(?, now(), ?, \N, \N, \N, ?, '
          . '\N, \N, ?, ?, \N, \N, \N)';
        my $sth = $core_dbh->prepare($sql);

        $sth->bind_param( 1, ++$analysis_id,          SQL_INTEGER );
        $sth->bind_param( 2, 'xrefcoordinatemapping', SQL_VARCHAR );
        $sth->bind_param( 3, 'xref_mapper.pl',        SQL_VARCHAR );
        $sth->bind_param( 4, $analysis_params,        SQL_VARCHAR );
        $sth->bind_param( 5, 'CoordinateMapper.pm',   SQL_VARCHAR );

        $sth->execute();
      }
    } ## end else [ if ( defined($analysis_id...
  } ## end if ( !defined($analysis_id...

  if ( defined($analysis_id) ) {
    log_progress( "Analysis ID                  is %d\n",
                  $analysis_id );
  }

  ######################################################################
  # Read and store available Xrefs from the Xref database.             #
  ######################################################################

  my %unmapped;
  my %mapped;

  my $xref_sql = qq(
    SELECT  coord_xref_id, source_id, accession
    FROM    coordinate_xref
    WHERE   species_id = ?
  );

  my $xref_sth = $xref_dbh->prepare($xref_sql);
  $xref_sth->bind_param( 1, $species_id, SQL_INTEGER );

  $xref_sth->execute();

  my $external_db_id;

  while ( my $xref = $xref_sth->fetchrow_hashref() ) {
    $external_db_id ||=
      $XrefMapper::BasicMapper::source_to_external_db{ $xref->{
        'source_id'} };
    $external_db_id ||= 11000;    # FIXME (11000 is 'UCSC')

    $unmapped{ $xref->{'coord_xref_id'} } = {
                   'external_db_id' => $external_db_id,
                   'accession'      => $xref->{'accession'},
                   'reason'         => 'No overlap',
                   'reason_full' =>
                     'No coordinate overlap with any Ensembl transcript'
    };
  }
  $xref_sth->finish();

  if(!defined($external_db_id)){
    die "External_db_id is undefined  species_id = $species_id\n";
  }

  ######################################################################
  # Do coordinate matching.                                            #
  ######################################################################

  my $core_db_adaptor =
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
                                   -host => $core_db->dbc()->host(),
                                   -port => $core_db->dbc()->port(),
                                   -user => $core_db->dbc()->username(),
                                   -pass => $core_db->dbc()->password(),
                                   -dbname => $core_db->dbc()->dbname(),
    );

  my $slice_adaptor = $core_db_adaptor->get_SliceAdaptor();
  my @chromosomes   = @{ $slice_adaptor->fetch_all('Chromosome') };

  my $sql = qq(
    SELECT    coord_xref_id, accession,
              txStart, txEnd,
              cdsStart, cdsEnd,
              exonStarts, exonEnds
    FROM      coordinate_xref
    WHERE     species_id = ?
    AND       chromosome = ? AND strand   = ?
    AND       ((txStart BETWEEN ? AND ?)        -- txStart in region
    OR         (txEnd   BETWEEN ? AND ?)        -- txEnd in region
    OR         (txStart <= ? AND txEnd >= ?))   -- region is fully contained
    ORDER BY  accession
  );

  foreach my $chromosome (@chromosomes) {
    my $chr_name = $chromosome->seq_region_name();

    log_progress( "Processing chromsome '%s'\n", $chr_name );

    my @genes = @{ $chromosome->get_all_Genes( undef, undef, 1 ) };

    log_progress( "There are %4d genes on chromosome '%s'\n",
                  scalar(@genes), $chr_name );

    while ( my $gene = shift(@genes) ) {
      my @transcripts = @{ $gene->get_all_Transcripts() };

      my %gene_result;

      foreach my $transcript ( sort { $a->start() <=> $b->start() }
                               @transcripts )
      {

        ################################################################
        # For each Ensembl transcript:                                 #
        #   1. Register all Ensembl exons in a RangeRegistry.          #
        #                                                              #
        #   2. Find all transcripts in the external database that are  #
        #      within the range of this Ensembl transcript.            #
        #                                                              #
        # For each of those external transcripts:                      #
        #   3. Calculate the overlap of the exons of the external      #
        #      transcript with the Ensembl exons using the             #
        #      overlap_size() method in the RangeRegistry.             #
        #                                                              #
        #   4. Register the external exons in their own RangeRegistry. #
        #                                                              #
        #   5. Calculate the overlap of the Ensembl exons with the     #
        #      external exons as in step 3.                            #
        #                                                              #
        #   6. Calculate the match score.                              #
        #                                                              #
        #   7. Decide whether or not to keep the match.                #
        ################################################################

        my @exons = @{ $transcript->get_all_Exons() };

        my %transcript_result;

        # '$rr1' is the RangeRegistry holding Ensembl exons for one
        # transcript at a time.
        my $rr1 = Bio::EnsEMBL::Mapper::RangeRegistry->new();

        my $coding_transcript;
        if ( defined( $transcript->translation() ) ) {
          $coding_transcript = 1;
        } else {
          $coding_transcript = 0;
        }

        foreach my $exon (@exons) {

          #-------------------------------------------------------------
          # Register each exon in the RangeRegistry.  Register both the
          # total length of the exon and the coding range of the exon.
          #-------------------------------------------------------------

          $rr1->check_and_register( 'exon', $exon->start(),
                                    $exon->end() );

          if (    $coding_transcript
               && defined( $exon->coding_region_start($transcript) )
               && defined( $exon->coding_region_end($transcript) ) )
          {
            $rr1->check_and_register(
                                'coding',
                                $exon->coding_region_start($transcript),
                                $exon->coding_region_end($transcript) );
          }
        }

        #---------------------------------------------------------------
        # Get hold of all transcripts from the external database that
        # overlaps with this Ensembl transcript.
        #---------------------------------------------------------------

        my $sth = $xref_dbh->prepare_cached($sql);
        $sth->bind_param( 1, $species_id,          SQL_INTEGER );
        $sth->bind_param( 2, $chr_name,            SQL_VARCHAR );
        $sth->bind_param( 3, $gene->strand(),      SQL_INTEGER );
        $sth->bind_param( 4, $transcript->start(), SQL_INTEGER );
        $sth->bind_param( 5, $transcript->end(),   SQL_INTEGER );
        $sth->bind_param( 6, $transcript->start(), SQL_INTEGER );
        $sth->bind_param( 7, $transcript->end(),   SQL_INTEGER );
        $sth->bind_param( 8, $transcript->start(), SQL_INTEGER );
        $sth->bind_param( 9, $transcript->end(),   SQL_INTEGER );

        $sth->execute();

        my ( $coord_xref_id, $accession, $txStart, $txEnd, $cdsStart,
             $cdsEnd, $exonStarts, $exonEnds );

        $sth->bind_columns(
                        \( $coord_xref_id, $accession, $txStart, $txEnd,
                           $cdsStart, $cdsEnd, $exonStarts, $exonEnds
                        ) );

        while ( $sth->fetch() ) {
          my @exonStarts = split( /,\s*/, $exonStarts );
          my @exonEnds   = split( /,\s*/, $exonEnds );
          my $exonCount = scalar(@exonStarts);

          # '$rr2' is the RangeRegistry holding exons from the external
          # transcript, for one transcript at a time.
          my $rr2 = Bio::EnsEMBL::Mapper::RangeRegistry->new();

          my $exon_match   = 0;
          my $coding_match = 0;

          my $coding_count = 0;

          for ( my $i = 0 ; $i < $exonCount ; ++$i ) {

            #-----------------------------------------------------------
            # Register the exons from the external database in the same
            # was as with the Ensembl exons, and calculate the overlap
            # of the external exons with the previously registered
            # Ensembl exons.
            #-----------------------------------------------------------

            my $overlap =
              $rr1->overlap_size( 'exon', $exonStarts[$i],
                                  $exonEnds[$i] );

            $exon_match +=
              $overlap/( $exonEnds[$i] - $exonStarts[$i] + 1 );

            $rr2->check_and_register( 'exon', $exonStarts[$i],
                                      $exonEnds[$i] );

            if ( !defined($cdsStart) || !defined($cdsEnd) ) {
              # Non-coding transcript.
            } else {
              my $codingStart = (   $exonStarts[$i] > $cdsStart
                                  ? $exonStarts[$i]
                                  : $cdsStart );
              my $codingEnd =
                ( $exonEnds[$i] < $cdsEnd ? $exonEnds[$i] : $cdsEnd );

              if ( $codingStart < $codingEnd ) {
                my $coding_overlap =
                  $rr1->overlap_size( 'coding', $codingStart,
                                      $codingEnd );

                $coding_match +=
                  $coding_overlap/( $codingEnd - $codingStart + 1 );

                $rr2->check_and_register( 'coding', $codingStart,
                                          $codingEnd );

                ++$coding_count;
              }
            }
          } ## end for ( my $i = 0 ; $i < ...

          my $rexon_match   = 0;
          my $rcoding_match = 0;

          my $rcoding_count = 0;

          foreach my $exon (@exons) {

            #-----------------------------------------------------------
            # Calculate the overlap of the Ensembl exons with the
            # external exons.
            #-----------------------------------------------------------

            my $overlap =
              $rr2->overlap_size( 'exon', $exon->start(),
                                  $exon->end() );

            $rexon_match +=
              $overlap/( $exon->end() - $exon->start() + 1 );

            if (    $coding_transcript
                 && defined( $exon->coding_region_start($transcript) )
                 && defined( $exon->coding_region_end($transcript) ) )
            {
              my $coding_overlap =
                $rr2->overlap_size( 'coding',
                                    $exon->coding_region_start(
                                                           $transcript),
                                    $exon->coding_region_end(
                                                            $transcript)
                );

              $rcoding_match +=
                $coding_overlap/
                ( $exon->coding_region_end($transcript) -
                  $exon->coding_region_start($transcript) +
                  1 );

              ++$rcoding_count;
            }
          } ## end foreach my $exon (@exons)

          #-------------------------------------------------------------
          # Calculate the match score.
          #-------------------------------------------------------------

          my $score = ( ( $exon_match + $ens_weight*$rexon_match ) +
                          $coding_weight*(
                              $coding_match + $ens_weight*$rcoding_match
                          )
            )/( ( $exonCount + $ens_weight*scalar(@exons) ) +
                  $coding_weight*(
                              $coding_count + $ens_weight*$rcoding_count
                  ) );

          if ( !defined( $transcript_result{$coord_xref_id} )
               || $transcript_result{$coord_xref_id} < $score )
          {
            $transcript_result{$coord_xref_id} = $score;
          }

        } ## end while ( $sth->fetch() )
        $sth->finish();

        #---------------------------------------------------------------
        # Apply transcript threshold and pick the best match(es) for
        # this transcript.
        #---------------------------------------------------------------

        my $best_score;
        foreach my $coord_xref_id (
             sort( { $transcript_result{$b} <=> $transcript_result{$a} }
                   keys(%transcript_result) ) )
        {
          my $score = $transcript_result{$coord_xref_id};

          if ( $score > $transcript_score_threshold ) {
            $best_score ||= $score;

            if ( sprintf( "%.3f", $score ) eq
                 sprintf( "%.3f", $best_score ) )
            {
              if ( exists( $unmapped{$coord_xref_id} ) ) {
                $mapped{$coord_xref_id} = $unmapped{$coord_xref_id};
                delete( $unmapped{$coord_xref_id} );
                $mapped{$coord_xref_id}{'reason'}      = undef;
                $mapped{$coord_xref_id}{'reason_full'} = undef;
              }

              push( @{ $mapped{$coord_xref_id}{'mapped_to'} }, {
                      'ensembl_id'          => $transcript->dbID(),
                      'ensembl_object_type' => 'Transcript'
                    } );

              # This is now a candidate Xref for the gene.
              if ( !defined( $gene_result{$coord_xref_id} )
                   || $gene_result{$coord_xref_id} < $score )
              {
                $gene_result{$coord_xref_id} = $score;
              }

            } elsif ( exists( $unmapped{$coord_xref_id} ) ) {
              $unmapped{$coord_xref_id}{'reason'} =
                'Was not best match';
              $unmapped{$coord_xref_id}{'reason_full'} =
                sprintf(
                       "Did not top best transcript match score (%.2f)",
                       $best_score );
              if ( !defined( $unmapped{$coord_xref_id}{'score'} )
                   || $score > $unmapped{$coord_xref_id}{'score'} )
              {
                $unmapped{$coord_xref_id}{'score'} = $score;
                $unmapped{$coord_xref_id}{'ensembl_id'} =
                  $transcript->dbID();
              }
            }

          } elsif ( exists( $unmapped{$coord_xref_id} )
                    && $unmapped{$coord_xref_id}{'reason'} ne
                    'Was not best match' )
          {
            $unmapped{$coord_xref_id}{'reason'} =
              'Did not meet threshold';
            $unmapped{$coord_xref_id}{'reason_full'} =
              sprintf( "Match score for transcript "
                         . "lower than threshold (%.2f)",
                       $transcript_score_threshold );
            if ( !defined( $unmapped{$coord_xref_id}{'score'} )
                 || $score > $unmapped{$coord_xref_id}{'score'} )
            {
              $unmapped{$coord_xref_id}{'score'} = $score;
              $unmapped{$coord_xref_id}{'ensembl_id'} =
                $transcript->dbID();
            }
          }
        } ## end foreach my $coord_xref_id (...

      } ## end foreach my $transcript ( sort...

      #-----------------------------------------------------------------
      # Pick the best match(es) for this gene.
      #-----------------------------------------------------------------

      # Do not want them on both Transcripts and Genes. For Biomart purposes. (Ianl)
      # So no need to find the best one for the gene.

#      my $best_score;
#      foreach my $coord_xref_id (
#                         sort( { $gene_result{$b} <=> $gene_result{$a} }
#                               keys(%gene_result) ) )
#      {
#        my $score = $gene_result{$coord_xref_id};
#
#        $best_score ||= $score;
#
#        if (
#           sprintf( "%.3f", $score ) eq sprintf( "%.3f", $best_score ) )
#        {
#          push( @{ $mapped{$coord_xref_id}{'mapped_to'} }, {
#                  'ensembl_id'          => $gene->dbID(),
#                  'ensembl_object_type' => 'Gene'
#                } );
#        }
#      }

    } ## end while ( my $gene = shift(...
  } ## end foreach my $chromosome (@chromosomes)

  # Make all dumps.  Order is important.
  dump_xref( $xref_filename, $xref_id, \%mapped, \%unmapped );
  dump_object_xref( $object_xref_filename, $object_xref_id, \%mapped );
  dump_unmapped_reason( $unmapped_reason_filename, $unmapped_reason_id,
                        \%unmapped, $core_dbh );
  dump_unmapped_object( $unmapped_object_filename, $unmapped_object_id,
                        $analysis_id, \%unmapped );

  if ($do_upload) {
    # Order is important!
    upload_data( 'unmapped_reason', $unmapped_reason_filename,
                 $external_db_id, $core_dbh );
    upload_data( 'unmapped_object', $unmapped_object_filename,
                 $external_db_id, $core_dbh );
    upload_data( 'object_xref', $object_xref_filename, $external_db_id,
                 $core_dbh );
    upload_data( 'xref', $xref_filename, $external_db_id, $core_dbh );
  }

  $sth_stat = $self->xref->dbc->prepare(
                          "INSERT INTO process_status (status, date) "
                        . "VALUES ('coordinate_xref_finished',now())" );
  $sth_stat->execute();
  $sth_stat->finish();

  $self->biomart_fix("UCSC","Translation","Gene");
  $self->biomart_fix("UCSC","Transcript","Gene");
  
} ## end sub run_coordinatemapping

#-----------------------------------------------------------------------

sub dump_xref {
  my ( $filename, $xref_id, $mapped, $unmapped ) = @_;

  ######################################################################
  # Dump for 'xref'.                                                   #
  ######################################################################

  my $fh = IO::File->new( '>' . $filename )
    or croak( sprintf( "Can not open '%s' for writing", $filename ) );

  log_progress( "Dumping for 'xref' to '%s'\n", $filename );

  foreach my $xref ( values( %{$unmapped} ), values( %{$mapped} ) ) {
    # Assign 'xref_id' to this Xref.
    $xref->{'xref_id'} = ++$xref_id;

    my $accession = $xref->{'accession'};

    my ($version) = ( $accession =~ /\.(\d+)$/ );
    $version ||= 0;

    $fh->printf("%d\t%d\t%s\t%s\t%d\t%s\t%s\t%s\n",
                $xref->{'xref_id'},
                $xref->{'external_db_id'},
                $accession,
                $accession,
                $version,
                '\N',
                'COORDINATE_OVERLAP',
                ''                                  # FIXME (possibly)
    );
  }
  $fh->close();

  log_progress("Dumping for 'xref' done\n");

} ## end sub dump_xref

#-----------------------------------------------------------------------

sub dump_object_xref {
  my ( $filename, $object_xref_id, $mapped ) = @_;

  ######################################################################
  # Dump for 'object_xref'.                                            #
  ######################################################################

  my $fh = IO::File->new( '>' . $filename )
    or croak( sprintf( "Can not open '%s' for writing", $filename ) );

  log_progress( "Dumping for 'object_xref' to '%s'\n", $filename );

  foreach my $xref ( values( %{$mapped} ) ) {
    foreach my $object_xref ( @{ $xref->{'mapped_to'} } ) {
      # Assign 'object_xref_id' to this Object Xref.
      $object_xref->{'object_xref_id'} = ++$object_xref_id;

      $fh->printf( "%d\t%d\t%s\t%d\t%s\t%s\n",
                   $object_xref->{'object_xref_id'},
                   $object_xref->{'ensembl_id'},
                   $object_xref->{'ensembl_object_type'},
                   $xref->{'xref_id'},
                   '\N',
                   '0' );
    }
  }
  $fh->close();

  log_progress("Dumping for 'object_xref' done\n");

} ## end sub dump_objexref

#-----------------------------------------------------------------------

sub dump_unmapped_reason {
  my ( $filename, $unmapped_reason_id, $unmapped, $core_dbh ) = @_;

  ######################################################################
  # Dump for 'unmapped_reason'.                                        #
  ######################################################################

  # Create a list of the unique reasons.
  my %reasons;

  foreach my $xref ( values( %{$unmapped} ) ) {
    if ( !exists( $reasons{ $xref->{'reason_full'} } ) ) {
      $reasons{ $xref->{'reason_full'} } = {
                                        'summary' => $xref->{'reason'},
                                        'full' => $xref->{'reason_full'}
      };
    }
  }

  my $fh = IO::File->new( '>' . $filename )
    or croak( sprintf( "Can not open '%s' for writing", $filename ) );

  log_progress( "Dumping for 'unmapped_reason' to '%s'\n", $filename );

  my $sth =
    $core_dbh->prepare(   'SELECT unmapped_reason_id '
                        . 'FROM unmapped_reason '
                        . 'WHERE full_description = ?' );

  foreach my $reason (
            sort( { $a->{'full'} cmp $b->{'full'} } values(%reasons) ) )
  {
    # Figure out 'unmapped_reason_id' from the core database.
    $sth->bind_param( 1, $reason->{'full'}, SQL_VARCHAR );

    $sth->execute();

    my $id;
    $sth->bind_col( 1, \$id );
    $sth->fetch();

    if ( defined($id) ) {
      $reason->{'unmapped_reason_id'} = $id;
    } else {
      $reason->{'unmapped_reason_id'} = ++$unmapped_reason_id;
    }

    $sth->finish();

    $fh->printf( "%d\t%s\t%s\n", $reason->{'unmapped_reason_id'},
                 $reason->{'summary'}, $reason->{'full'} );

  }
  $fh->close();

  log_progress("Dumping for 'unmapped_reason' done\n");

  # Assign reasons to the unmapped Xrefs from %reasons.
  foreach my $xref ( values( %{$unmapped} ) ) {
    $xref->{'reason'}      = $reasons{ $xref->{'reason_full'} };
    $xref->{'reason_full'} = undef;
  }

} ## end sub dump_unmapped_reason

#-----------------------------------------------------------------------

sub dump_unmapped_object {
  my ( $filename, $unmapped_object_id, $analysis_id, $unmapped ) = @_;

  ######################################################################
  # Dump for 'unmapped_object'.                                        #
  ######################################################################

  my $fh = IO::File->new( '>' . $filename )
    or croak( sprintf( "Can not open '%s' for writing", $filename ) );

  log_progress( "Dumping for 'unmapped_object' to '%s'\n", $filename );

  foreach my $xref ( values( %{$unmapped} ) ) {
    # Assign 'unmapped_object_id' to this Xref.
    $xref->{'unmapped_object_id'} = ++$unmapped_object_id;

    $fh->printf(
      "%d\t%s\t%s\t%d\t%s\t%d\t%s\t%s\t%s\t%s\t%s\n",
      $xref->{'unmapped_object_id'},
      'xref',
      $analysis_id || '\N',    # '\N' (NULL) means no analysis exists
                               # and uploading this table will fail.
      $xref->{'external_db_id'},
      $xref->{'accession'},
      $xref->{'reason'}->{'unmapped_reason_id'}, (
        defined( $xref->{'score'} )
        ? sprintf( "%.3f", $xref->{'score'} )
        : '\N'
      ),
      '\N',
      $xref->{'ensembl_id'} || '\N',
      ( defined( $xref->{'ensembl_id'} ) ? 'Transcript' : '\N' ),
      '\N' );
  }
  $fh->close();

  log_progress("Dumping for 'unmapped_object' done\n");

} ## end sub dump_unmapped_object

#-----------------------------------------------------------------------

sub upload_data {
  my ( $table_name, $filename, $external_db_id, $dbh ) = @_;

  ######################################################################
  # Upload data from a file to a table.                                #
  ######################################################################

  if ( !-r $filename ) {
    croak( sprintf( "Can not open '%s' for reading", $filename ) );
  }

  my $cleanup_sql = '';
  if ( $table_name eq 'unmapped_reason' ) {
    $cleanup_sql = qq(
    DELETE  ur
    FROM    unmapped_object uo,
            unmapped_reason ur
    WHERE   uo.external_db_id       = ?
    AND     ur.unmapped_reason_id   = uo.unmapped_reason_id);
  } elsif ( $table_name eq 'unmapped_object' ) {
    $cleanup_sql = qq(
    DELETE  uo
    FROM    unmapped_object uo
    WHERE   uo.external_db_id = ?);
  } elsif ( $table_name eq 'object_xref' ) {
    $cleanup_sql = qq(
    DELETE  ox
    FROM    xref x,
            object_xref ox
    WHERE   x.external_db_id    = ?
    AND     ox.xref_id          = x.xref_id);
  } elsif ( $table_name eq 'xref' ) {
    $cleanup_sql = qq(
    DELETE  x
    FROM    xref x
    WHERE   x.external_db_id    = ?);
  } else {
    croak( sprintf( "Table '%s' is unknown\n", $table_name ) );
  }

  my $load_sql =
    sprintf( "LOAD DATA LOCAL INFILE ? REPLACE INTO TABLE %s",
             $table_name );

  log_progress(
          "Removing old data (external_db_id = '%d') from table '%s'\n",
          $external_db_id, $table_name );

  my $rows = $dbh->do( $cleanup_sql, undef, $external_db_id )
    or croak( $dbh->strerr() );

  log_progress( "Removed %d rows\n", $rows );

  log_progress( "Uploading for '%s' from '%s'\n",
                $table_name, $filename );

  $rows = $dbh->do( $load_sql, undef, $filename )
    or croak( $dbh->errstr() );

  $dbh->do("OPTIMIZE TABLE $table_name") or croak( $dbh->errstr() );

  log_progress( "Uploading for '%s' done (%d rows)\n",
                $table_name, $rows );

} ## end sub upload_data

#-----------------------------------------------------------------------

sub log_progress {
  my ( $fmt, @params ) = @_;
  
  return if (!$verbose);
  printf( STDERR "COORD==> %s", sprintf( $fmt, @params ) );
}

1;
