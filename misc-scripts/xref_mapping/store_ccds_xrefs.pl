#!/usr/bin/env perl

# $Source: /cvsroot/ensembl/ensembl-personal/genebuilders/ccds/scripts/store_ccds_xrefs.pl,v $
# $Revision: 1.13 $

=pod

=head1 NAME

  store_ccds_xrefs.pl

=head1 SYNOPSIS

  Make CCDS Xrefs.

=head1 DESCRIPTION

  Will store the Ensembl transcript stable_id that matches the ccds structure.
  Originally written for Ian Longden. Based on make_enst_to_ccds.pl

=head1 ARGUMENTS

  perl store_ccds_xrefs.pl 
         -ccds_dbname
         -ccds_host
         -ccds_port
         -ccds_user
         -ccds_pass
         -dbname
         -host
         -port
         -user
         -verbose
         -species
         -path
         -write
         -delete_old

=head1 EXAMPLE
 
  perl $ENSEMBL_PERSONAL/genebuilders/ccds/scripts/store_ccds_xrefs.pl -ccds_dbname db8_human_vega_61 \
  -ccds_host genebuild7 -ccds_port 3306 -ccds_user user -ccds_pass password              \
  -dbname homo_sapiens_core_61_37f -host ens-staging1 -port 3306 -user ensro -verbose    \
  -species human -path GRCh37 -write -delete_old

=cut


use warnings;
use strict;

use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(warning throw);

# db of CCDS strcutures
my $ccds_host   = '';
my $ccds_port   = '3306';
my $ccds_user   = 'ensro';
my $ccds_pass   = undef;
my $ccds_dbname = '';

# db of Ensembl (protein_coding) genes
my $host   = 'ens-staging';
my $port   = '';
my $user   = 'ensro';
my $dbname = '';

my $path    = 'GRCh37';
my $species = 'human';
my $verbose;
my $write;
my $delete_old;

&GetOptions( 'ccds_host=s'   => \$ccds_host,
             'ccds_port=s'   => \$ccds_port,
             'ccds_user=s'   => \$ccds_user,
             'ccds_pass=s'   => \$ccds_pass,
             'ccds_dbname=s' => \$ccds_dbname,
             'host=s'        => \$host,
             'port=s'        => \$port,
             'user=s'        => \$user,
             'dbname=s'      => \$dbname,
             'path=s'        => \$path,
             'species=s'     => \$species,
             'verbose'       => \$verbose,
             'delete_old'    => \$delete_old,
             'write'         => \$write, );

if ( !defined $species ) {
  throw("Please define species as human or mouse");
} else {
  $species =~ s/\s//g;
  if ( $species =~ /^human$/i || $species =~ /^mouse$/i ) {
    # we're ok
    print "Species is *$species*\n";
  } else {
    throw("Species must be defined as human or mouse");
  }
}

# we want to keep a record of any polymorphic pseudogenes for havana
# let's not write a file until the end though since they are not
# common
my @polymorphic_pseudogene;


# connect to dbs
my $db =
  new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $host,
                                      -user   => $user,
                                      -port   => $port,
                                      -dbname => $dbname );

my $ccds_db =
  new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $ccds_host,
                                      -user   => $ccds_user,
                                      -pass   => $ccds_pass,
                                      -port   => $ccds_port,
                                      -dbname => $ccds_dbname );
$ccds_db->dnadb($db);
my $ccds_sa = $ccds_db->get_SliceAdaptor;
my $outdea  = $ccds_db->get_DBEntryAdaptor;
my $sa      = $db->get_SliceAdaptor;

###
# delete old ones if delete_old set
###
if($write and $delete_old){
  my $sth = $outdea->prepare('delete ox from xref x, object_xref ox, external_db e where x.xref_id = ox.xref_id and x.external_db_id = e.external_db_id and e.db_name like "Ens_%"');

  $sth->execute || die "Could not delete old object_xrefs";

  $sth = $outdea->prepare('delete x from xref x, external_db e where x.external_db_id = e.external_db_id and e.db_name like "Ens_%"');

  $sth->execute || die "Could not delete ols xrefs";
}

# # #
# Loop thru toplevels
# # #
# maybe should use toplevel instead of chromosome?
foreach my $chr ( @{ $ccds_sa->fetch_all('chromosome') } ) {
  print "Doing chromosome " . $chr->name . "\n" if ($verbose);

  # fetch all CCDS structures on slice
  foreach my $ccds_gene ( @{ $chr->get_all_Genes( undef, undef, 1 ) } ) {

    # make sure genes are al on chr level
    $ccds_gene = $ccds_gene->transform( 'chromosome', $path );

    # loop thru all CCDS transcripts
    foreach my $ccds_trans ( @{ $ccds_gene->get_all_Transcripts() } ) {
      print "=> doing ccds trans "
        . $ccds_trans->dbID
        . ": start "
        . $ccds_trans->start
        . " stop "
        . $ccds_trans->end
        . " strand "
        . $ccds_trans->strand . " \n"
        if ($verbose);

      # find the ccds_id
      my $ccds_id;
      my @db_entries = @{ $ccds_trans->get_all_DBEntries('CCDS') };
      my %xref_hash;

      foreach my $dbe (@db_entries) {
        print "dbe " . $dbe->display_id . " " . $dbe->dbname . "\n";
      }

      # store unique CCDS xrefs for the transcript
      foreach my $entry (@db_entries) {
        $xref_hash{ $entry->display_id() } = 1;
      }

      # we should not have more than one CCDS id
      # associated with a transcript
      if ( scalar keys %xref_hash != 1 ) {
        foreach my $entry ( keys %xref_hash ) {
          print "  Dodgy xref : " . $entry . "\n";
        }
        throw(   "Something odd going on: Transcript dbID "
               . $ccds_trans->dbID . " has "
               . scalar( keys %xref_hash )
               . " xrefs" );
      } else {
        # all is good; CCDS transcript only has 1 CCDS xref
        foreach my $entry ( keys %xref_hash ) {
          $ccds_id = $entry;
          print "=> on ccds $ccds_id\n" if ($verbose);
        }
      }

      # define the genomic location that we're working in
      # ie. where the CCDS transcript is
      my $chr_name = $ccds_trans->slice->seq_region_name;
      my $start    = $ccds_trans->start();
      my $end      = $ccds_trans->end();

      # now fetch the slice out of ensembl db
      my $slice =
        $sa->fetch_by_region( 'chromosome', $chr_name, $start, $end, '1',
                              $path );
      print "  Ensembl slice name " . $slice->name . "\n" if ($verbose);

      # get ccds coding exons
      my @ccds_exons = @{ $ccds_trans->get_all_translateable_Exons() };
      print "  have " . @ccds_exons . " ccds coding exons\n" if ($verbose);

      # get all Ensembl genes overlapping the CCDS regions
      foreach my $gene ( @{ $slice->get_all_Genes( undef, undef, 1 ) } ) {

        # only look at protein_coding genes
        next unless ( $gene->biotype =~ /protein_coding/ || $gene->biotype =~ /polymorphic_pseudogene/);
                                  
        # debug
#        next if  $gene->biotype =~ /protein_coding/ ;
        
        # keep a record if it is a polymorphic pseudogene - these will need to be sent to havana
        if ($gene->biotype =~ /polymorphic_pseudogene/) {
          print STDERR " found a poly pseudo gene\n" if ($verbose);
          push @polymorphic_pseudogene, $ccds_id;
        }

        # make sure ensembl gene also on chr level
        print "  on ensembl gene " . $gene->display_id . "\n" if ($verbose);
        $gene = $gene->transform( 'chromosome', $path );

        # loop thru ensembl transcripts
        foreach my $trans ( @{ $gene->get_all_Transcripts } ) {
          print "  on ensembl trans " . $trans->display_id . "\n"
            if ($verbose);

          # get ensembl coding exons
          my @exons = @{ $trans->get_all_translateable_Exons() };
          print "  have " . @exons . " ensembl coding exons\n" if ($verbose);

        # loop thru ensembl coding exons and make sure they all match the ccds
        # exons exactly
          my $match = 0;

          if ( scalar @exons == scalar @ccds_exons ) {
            for ( my $i = 0 ; $i < scalar(@exons) ; $i++ ) {
              # print "    Ensembl start ".$exons[$i]->start." end ".$exons[$i]->end.
              #       " CCDS start ".$ccds_exons[$i]->start." end ".$ccds_exons[$i]->end."\n";

              if (    $ccds_exons[$i]->start == $exons[$i]->start
                   && $ccds_exons[$i]->end == $exons[$i]->end
                   && $ccds_exons[$i]->strand == $exons[$i]->strand )
              {

                $match++;
              } #else {
                # print "no match ".$ccds_exons[$i]->start." != ".$exons[$i]->start." or ".
                #	$ccds_exons[$i]->end." != ".$exons[$i]->end."\n";
                #}
            }

            if ( $match == scalar @exons ) {
              print "MATCH\t" . $trans->stable_id . "\t" . $ccds_id . "\n"
                if ($verbose);
              store_ensembl_xref( $outdea, $species, $ccds_trans,
                                  $trans->stable_id, $write );
              store_ensembl_xref( $outdea,
                                  $species,
                                  $ccds_trans->translation,
                                  $trans->translation->stable_id,
                                  $write );
            } else {
              print "  no match ($match)\t"
                . $trans->stable_id . "\t"
                . $ccds_id . "\n"
                if ($verbose);
            }
          } ## end if ( scalar @exons == ...)
        } ## end foreach my $trans ( @{ $gene...})
      } ## end foreach my $gene ( @{ $slice...})
    } ## end foreach my $ccds_trans ( @{...})
  } ## end foreach my $ccds_gene ( @{ ...})
} ## end foreach my $chr ( @{ $ccds_sa...})

# report polymorphic pseudogenes
if (@polymorphic_pseudogene) {
  for my $display_id (@polymorphic_pseudogene) {
    print STDERR $display_id." matches a polymorphic pseudogene\n";
  }
} else {
  print STDERR "Found 0 polymorphic pseudogenes\n";
}


sub store_ensembl_xref {
  my ( $dbea, $species, $ccds_trans, $ensembl_trans_stable_id, $write ) = @_;

  if ( ref($ccds_trans) eq "Bio::EnsEMBL::Transcript" ) {

    my $external_db;
    if ( $species =~ /^human$/i ) {
      $external_db = 'Ens_Hs_transcript';
    } elsif ( $species =~ /^mouse$/i ) {
      $external_db = 'Ens_Mm_transcript';
    }
    # make an xref

    my $entry =
      new Bio::EnsEMBL::DBEntry( -adaptor    => $dbea,
                                 -primary_id => $ensembl_trans_stable_id,
                                 -display_id => $ensembl_trans_stable_id,
                                 -version    => 0,
                                 -dbname     => $external_db);
    # store xref
    $dbea->store( $entry, $ccds_trans->dbID, 'Transcript' ) if ($write);

  } elsif ( ref($ccds_trans) eq "Bio::EnsEMBL::Translation" ) {
    my $external_db;
    if ( $species =~ /^human$/i ) {
      $external_db = 'Ens_Hs_translation';
    } elsif ( $species =~ /^mouse$/i ) {
      $external_db = 'Ens_Mm_translation';
    }
    # make an xref
    my $entry =
      new Bio::EnsEMBL::DBEntry( -adaptor    => $dbea,
                                 -primary_id => $ensembl_trans_stable_id,
                                 -display_id => $ensembl_trans_stable_id,
                                 -version    => 0,
                                 -dbname     => $external_db);
    # store xref
    $dbea->store( $entry, $ccds_trans->dbID, 'Translation' ) if ($write);

  } else {
    throw("Not a Transcript or Translation ");
  }

  return;
} ## end sub store_ensembl_xref
