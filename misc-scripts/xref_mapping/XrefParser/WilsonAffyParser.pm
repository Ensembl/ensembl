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

package XrefParser::WilsonAffyParser;

use strict;
use warnings;
use Carp;
use base qw( XrefParser::BaseParser );

sub run {
  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $files        = $ref_arg->{files};
  my $verbose      = $ref_arg->{verbose};

  if((!defined $source_id) or (!defined $species_id) or (!defined $files) ){
    croak "Need to pass source_id, species_id and files as pairs";
  }
  $verbose |=0;

  my @xrefs = $self->create_xrefs($source_id, $species_id, @{$files}[0], $verbose);

  if(!@xrefs){
    return 1; #  1 error
  }
  # upload
  if(!defined($self->upload_xref_object_graphs(@xrefs))){
    return 1;
  }
  return 0;

}

sub create_xrefs {

  my ($self, $source_id, $species_id, $file, $verbose) = @_;

  my ($count, $noseq, $direct) = (0,0,0);

  local $| = 1; # don't buffer

  my @xrefs;

  my $file_io = $self->get_filehandle($file);

  if ( !defined $file_io ) {
    print STDERR "ERROR: Could not open $file\n";
    return 1;    # 1 error
  }

  $file_io->getline();    # skip first line

  while ( $_ = $file_io->getline() ) {
    #last if ($count > 200);
    my $xref;

    my @fields = split /\t/;

    # first field (probe_set) is accession
    my $acc = $fields[0];
    chomp($acc);
    $acc =~ s/\"//g;


    # get linked accession (may be RefSeq or EMBL or ensembl)
    my $target = $fields[2];
    chomp($target);
    $target =~ s/\"//g;

    # Create direct xrefs for mappings to Ensembl transcripts
    if ($target =~ /ENSGALT/) {

      # remove version if present
      ($target) = $target =~ /([^.]*)\.([^.]*)/;

      # add xref - not we're assuming it doesn't already exist;
      # may need to check like in CCDS parser
      my $xref_id = $self->add_xref({ acc => $acc,
				      version => 0,
				      label => $acc,
				      desc => "$target direct mapping",
				      source_id => $source_id,
				      species_id => $species_id} );
      $self->add_direct_xref($xref_id, $target, "transcript", "");
      $direct++;

    } else {

      # fetch sequence for others (EMBL ESTs and RefSeqs - pfetch will handle these)
      system ("pfetch -q $target > seq.txt");

      my $seq_io = $self->get_filehandle('seq.txt');

      my $seq = $seq_io->getline();
      $seq_io->close();

      chomp($seq);

      if ($seq && $seq !~ /no match/) {

	$xref->{ACCESSION} = $acc;
	$xref->{SEQUENCE} = $seq;
	$xref->{LABEL} = $acc;
	$xref->{SOURCE_ID} = $source_id;
	$xref->{SPECIES_ID} = $species_id;
	$xref->{SEQUENCE_TYPE} = 'dna';
	$xref->{STATUS} = 'experimental';

	# Add description noting where the mapping came from
	$xref->{DESCRIPTION} = $target . " used as mapping target";

	$count++;

	print "$count " if (($count % 100 == 0) and $verbose);

	push @xrefs, $xref;

      } else {

	print STDERR "Couldn't get sequence for $target\n";
	$noseq++;

      }

    }

  }

  $file_io->close();

  if($verbose){
    print "\n\nParsed $count primary xrefs.\n";
    print "Couldn't get sequence for $noseq primary_xrefs\n" if ($noseq);
    print "Added $direct direct xrefs.\n";
  }

  return \@xrefs;

}

1;
