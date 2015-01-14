#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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


=head1 NAME

parse_embl.pl

=head1 SYNOPSIS

 parse_embl -dbname <db name> [-host <db host>] [-port <db port>]
  [-user <db user>] [-pass <db pass>]
  [-genetype <gene type>]
  -codontable <int> -emblfile <file name>

=head1 DESCRIPTION

This is a script originally written by Keith James to load a pathogen
genome from an EMBL flat file into an Otter (modified ensembl) database.
It has been has been stripped down and modified to make it more general
purpose.

This script will peform simple loading of an annotated genome in EMBL
flat file format into an Ensembl database. The script loads a fairly
minimal set of gene information from an EMBL file, but could be
extended to do more.

Certain assumptions are made about the format of the EMBL file such
that it may require pre-processing before loading - these are listed
below. A number of non-standard tags (EMBL qualifiers) are used. I
recommend using BioPerl to write a suitable pre-processor.

=over

=item *

There may be multiple chromosomes/plasmids per EMBL file.

=item *

Each chromosome/plasmid must have a source feature defining the
chromosome/plasmid name and type. They must all have different
names. See EMBL feature table definition for 'chromosome' and
'plasmid' tags.

=item *

Protein coding genes are annotated as CDS features.

=item *

All gene synonyms are defined by multiple 'gene' and/or 'synonym' tags.

=item *

The frame of a CDS feature is defined by a 'codon_start' tag if it is
not in frame 0.

=item *

If a CDS feature has a 'pseudo' tag it is taken to be a pseudogene.

=item *

If a CDS feature has a 'partial' tag but fuzzy ends are not defined in
its location it is assumed the ends have not been found.

=item *

Non-standard translation tables are currently not handled.

=item *

Annotator remarks about genes are taken from 'note' tags of
features. Remarks about transcripts are taken from 'product', 'class'
and 'EC_number' tags of features.

=back

Redirect STDERR to a file to obtain a rough log of the loading
operations.

Some functions in the script are simply stubs which may be implemented
later. Some user assembly may be required.

=head1 AUTHOR - Keith James

Email kdj@sanger.ac.uk

=head1 CONTRIBUTORS

Stephen Searle, sjms@sanger.ac.uk

Graham McVicker

=head1 METHODS

=cut

use strict;
use warnings;
use Getopt::Long;

use Bio::SeqIO;
use Bio::Location::Simple;

use Bio::EnsEMBL::DBSQL::DBAdaptor;

use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Analysis;

use Bio::EnsEMBL::Utils::Exception qw(throw);

###
### Script globals.
###

my $MAX_SEQUENCE_LEN = 2e6;
my $CURATED_LOGIC_NAME = 'Genome Reviews';

###
### Command line options
###
my $emblfile      = undef;
my $host          = '127.0.0.1';
my $user          = undef;
my $pass          = undef;
my $port          = 3306;
my $dbname        = undef;
my $gene_type     = 'known';

GetOptions('emblfile:s'     => \$emblfile,
           'dbname:s'       => \$dbname,
           'user:s'         => \$user,
           'pass:s'         => \$pass,
           'port:n'         => \$port,
           'host:s'         => \$host,
           'genetype:s'     => \$gene_type);

throw("-dbname argument required") if(!defined($dbname));
throw("-emblfile argument required") if(!defined($emblfile));

###
### Get Ensembl DB adaptor
###
my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(-host   => $host,
                                             -user   => $user,
                                             -pass   => $pass,
                                             -port   => $port,
                                             -dbname => $dbname);


###
### Create new analysis object indicating a curated origin
###
my $analysis = Bio::EnsEMBL::Analysis->new(-logic_name => $CURATED_LOGIC_NAME);

###
### Create sequence object and load it
###
my $seqi = Bio::SeqIO->new(-file   => $emblfile,
                           -format => 'embl');
###
### Main loading loop
###
while (my $seq = $seqi->next_seq) {
  # Split sequence and return a slice
  my $slice = load_sequence($dba, $seq);

  foreach my $f ($seq->get_SeqFeatures) {
    if ($f->primary_tag =~ /CDS/) {
      store_gene($dba, $f, $slice);
    }

    if ($f->primary_tag =~ /repeat/) {
      # store_repeat($f, $slice);
    }
  }
}


=head2 store_gene

 Title    : store_gene
 Function : Stores a gene feature (both protein- and RNA-coding)
 Returns  : void
 Argument : Bio::SeqFeatureI, Bio::EnsEMBL::DBSQL::SliceAdaptor

=cut

sub store_gene {
  my ($db, $f, $slice) = @_;

  my $gene_stable_id = get_gene_stable_id($f);

  my $gene = Bio::EnsEMBL::Gene->new();

  $gene->analysis($analysis);
  $gene->type($gene_type);
  $gene->stable_id($gene_stable_id);
  $gene->version(1);
  $gene->slice($slice);

  print STDERR sprintf("Found CDS with ID %s\n", $gene_stable_id);

  my $tcount = 0;
  my $transcript_id = sprintf("%s.%d", $gene_stable_id, $tcount++);

  my $transcript = Bio::EnsEMBL::Transcript->new;
  $transcript->stable_id($transcript_id);
  $transcript->version(1);
  $transcript->slice($slice);

  $gene->add_Transcript($transcript);

  # Add the exons
  my @exons = create_exons($f);
  my $ecount = 0;
  foreach my $exon (@exons) {
    $exon->stable_id(sprintf("%s.%d", $gene_stable_id, $ecount++));
    $exon->version(1);
    $exon->slice($slice);
    $transcript->add_Exon($exon);
  }

  if ($f->primary_tag =~ /RNA/) {
    foreach my $exon (@exons) {
      $exon->phase(-1);
      $exon->end_phase(-1);
    }
  }

  # Handle protein CDS features
  if ($f->primary_tag eq 'CDS') {
    # Based on /codon_start EMBL qualifier
    my $frame = get_initial_frame($f);

    @exons = @{ $transcript->get_all_Exons };

    # This code assumes no UTRs
    my $translation = Bio::EnsEMBL::Translation->new;
    my $rcount = 0;
    $translation->stable_id(sprintf("%s.%d", $gene_stable_id, $rcount++));
    $translation->version(1);
    $translation->start_Exon($exons[0]);
    $translation->start(1 + $frame);
    $translation->end_Exon($exons[$#exons]);
    $translation->end($translation->end_Exon->length);

    set_exon_phases($translation, @exons);

    foreach my $exon (@exons) {
      print STDERR sprintf("Added exon start: %d end: %d strand: %d phase: %d\n",
                           $exon->start, $exon->end, $exon->strand, $exon->phase);
    }

    $transcript->translation($translation);

    my $mrna_seq =
      Bio::Seq->new(-seq      => $transcript->translateable_seq,
                    -moltype  => "dna",
                    -alphabet => 'dna',
                    -id       => $translation->stable_id);

    # Translate args: stop char, unknown aa char, frame,
    # table, full CDS, throw
    my $aa_seq = $transcript->translate()->seq();

    print STDERR sprintf("Translation is: %s\n", $aa_seq);

    if ($aa_seq =~ /\*/) {
      print STDERR sprintf("Failed translating %s after phase setting\n",
                           $translation->stable_id);
    }
  }

  eval {
    $db->get_GeneAdaptor->store($gene);
  };
  throw(sprintf("Failed loading %s\n%s\n", $gene_stable_id, $@)) if($@);
}



=head2 get_gene_stable_id

  Arg [1]    : Bio::SeqFeatureI $f
  Example    : my $stable_id = get_gene_stable_id($f);
  Description: Tries to get a sensible stable identifier from a bioperl feature
               created from an embl flat file.
  Returntype : string
  Exceptions : throw if cannot determine stable identifier
  Caller     : 

=cut

sub get_gene_stable_id {
  my $f = shift;

  my $stable_id;

  if($f->has_tag('locus_tag')) {
    my @vals = $f->get_tag_values('locus_tag');
    ($stable_id) = split(/\s+/,$vals[0]);
    return $stable_id if($stable_id);
  }

  if($f->has_tag('gene')) {
    my @vals = $f->get_tag_values('gene');
    ($stable_id) = split(/\s+/,$vals[0]);
    return $stable_id if($stable_id);
  }

  throw("Could not determine gene identifier\n");
}


=head2 create_exons

 Title    : create_exons
 Function : Returns a list of exons created from the location of
            feature.
 Returns  : List of Bio::EnsEMBL::Exons
 Argument : Bio::SeqFeatureI

=cut

sub create_exons {
  my $f = shift;

  my @exons;

  foreach my $loc ($f->location->each_Location) {
    push(@exons, Bio::EnsEMBL::Exon->new(-start  => $loc->start,
                                         -end    => $loc->end,
                                         -strand => $loc->strand));

    print STDERR sprintf("Creating exon at %d..%d on strand %d\n",
                         $loc->start, $loc->end, $loc->strand);
  }

  if ($f->has_tag('pseudo')) {
    foreach my $exon (@exons) {
      $exon->phase(-1);
      $exon->end_phase(-1);
    }
  }
  else {
    foreach my $exon (@exons) {
      $exon->end_phase(0);
    }
  }

  return @exons;
}

=head2 get_initial_frame

 Title    : get_initial_frame
 Function : Returns the frame specified by the codon_start tag of
            feature.
 Returns  : int
 Argument : Bio::SeqFeatureI

=cut

sub get_initial_frame {
  my $f = shift;

  my $frame = 0;

  if ($f->has_tag('codon_start')) {
    my @vals = $f->get_tag_values('codon_start');
    $frame = $vals[0] - 1;
  }

  return $frame;
}

=head2 set_exon_phases

 Title    : set_exon_phases
 Function : Sets the start and end phases of exons.
 Returns  : void
 Argument : Bio::Otter::AnnotatedTranscript, Bio::Ensembl::Exons

=cut

sub set_exon_phases {
  my $translation = shift;
  my @exons = @_;

  my $found_start = 0;
  my $found_end   = 0;
  my $phase       = 0;

  foreach my $exon (@exons) {
    # Internal and end exons
    if ($found_start && ! $found_end) {
      $exon->phase($phase);
      $exon->end_phase(($exon->length + $exon->phase) % 3);
      $phase = $exon->end_phase;
    }

    if ($translation->start_Exon == $exon) {
      $exon->phase($phase);
      $exon->end_phase((($exon->length - $translation->start + 1)
                        + $exon->phase) % 3);
      $phase = $exon->end_phase;
      $found_start = 1;
    }

    if ($translation->end_Exon == $exon) {
      $found_end = 1;
    }
  }
}



=head2 get_seq_region_data

  Arg [1]    : Bio::SeqI
  Example    : my ($csname, $seq_region_name) = get_seq_region_data($seq);
  Description: Gets the coordinate system name (e.g. 'chromosome', 'plasmid')
               and the name of the seq_region
  Returntype : pair of strings
  Exceptions : none
  Caller     : load_sequence

=cut

sub get_seq_region_data {
  my $seq = shift;
  my $type;
  my $name;

  foreach my $f ($seq->get_SeqFeatures) {
    if ($f->primary_tag eq 'source') {
      if ($f->start == 1 and $f->end == $seq->length) {
        if ($f->has_tag('chromosome')){
          $type = 'chromosome';
          my @vals = $f->get_tag_values('chromosome');
          $name = $vals[0];
          $name =~ s/chr(omosome)\s+//i; #strip off chromosome prefix if present
          last;
        }

        if ($f->has_tag('plasmid')) {
          $type = 'plasmid';
          my @vals = $f->get_tag_values('plasmid');
          $name = $vals[0];
          $name =~ s/plasmid\s+//i; #strip off plasmid prefix if present
          last;
        }
      }
    }
  }

  return ($type, $name);
}


=head2 load_sequence

  Arg [1]    : Bio::EnsEMBL::DBAdaptor $db
  Arg [2]    : Bio::SeqI $seq
  Example    : load_sequence
  Description: Loads the sequence for this organism into the database, creating
               Coordinate system and seq_region as necessary
  Returntype : none
  Exceptions : none
  Caller     : general

=cut

sub load_sequence {
  my $db = shift;
  my $seq = shift;

  my ($csname, $seq_region) = get_seq_region_data($seq);

  my $csa = $db->get_CoordSystemAdaptor();
  my $cs = $csa->fetch_by_name($csname);

  if(!$cs) {
    $cs = Bio::EnsEMBL::CoordSystem->new(-NAME           => $csname,
                                          -SEQUENCE_LEVEL => 1,
                                          -RANK           => 1,
                                          -DEFAULT        => 1);
    $csa->store($cs);
  }

  my $slice = Bio::EnsEMBL::Slice->new(-SEQ_REGION_NAME => $seq_region,
                                       -COORD_SYSTEM    => $cs,
                                       -START           => 1,
                                       -END             => $seq->length(),
                                       -SEQ_REGION_LENGTH => $seq->length());


  my $slice_adaptor = $db->get_SliceAdaptor();

  my $sequence = $seq->seq();
  $slice_adaptor->store($slice, \$sequence);

  return $slice;
}

