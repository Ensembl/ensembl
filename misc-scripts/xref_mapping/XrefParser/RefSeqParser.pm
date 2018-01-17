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

# Parse RefSeq files to create xrefs.

package XrefParser::RefSeqParser;

use strict;
use warnings;
use Carp;
use File::Basename;

use base qw( XrefParser::BaseParser );


sub run {

  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $files        = $ref_arg->{files};
  my $release_file = $ref_arg->{rel_file};
  my $verbose      = $ref_arg->{verbose};
  my $dbi          = $ref_arg->{dbi};
  $dbi = $self->dbi unless defined $dbi;

  if((!defined $source_id) or (!defined $species_id) or (!defined $files) ){
    croak "Need to pass source_id, species_id and files as pairs";
  }
  $verbose |=0;

  my @files = @{$files};


  my $peptide_source_id =
    $self->get_source_id_for_source_name('RefSeq_peptide', undef, $dbi);
  my $mrna_source_id =
    $self->get_source_id_for_source_name('RefSeq_mRNA', undef, $dbi);
  my $ncrna_source_id =
    $self->get_source_id_for_source_name('RefSeq_ncRNA', undef, $dbi);

  my $pred_peptide_source_id =
    $self->get_source_id_for_source_name('RefSeq_peptide_predicted', undef, $dbi);
  my $pred_mrna_source_id =
    $self->get_source_id_for_source_name('RefSeq_mRNA_predicted','refseq', $dbi);
  my $pred_ncrna_source_id =
    $self->get_source_id_for_source_name('RefSeq_ncRNA_predicted', undef, $dbi);

  if($verbose){
    print "RefSeq_peptide source ID = $peptide_source_id\n";
    print "RefSeq_mRNA source ID = $mrna_source_id\n";
    print "RefSeq_ncRNA source ID = $ncrna_source_id\n";
    print "RefSeq_peptide_predicted source ID = $pred_peptide_source_id\n";
    print "RefSeq_mRNA_predicted source ID = $pred_mrna_source_id\n" ;
    print "RefSeq_ncRNA_predicted source ID = $pred_ncrna_source_id\n" ;
  }

    my @xrefs;
    foreach my $file (@files) {

        my $xrefs =
          $self->create_xrefs( $peptide_source_id,
                               $pred_peptide_source_id,
			       $mrna_source_id, $ncrna_source_id,
			       $pred_mrna_source_id, $pred_ncrna_source_id,
                               $file,
                               $species_id, $dbi );

        if ( !defined($xrefs) ) {
            return 1;    #error
        }
	print "Read " . scalar(@$xrefs) ." xrefs from $file\n" if($verbose);

        push @xrefs, @{$xrefs};
    }

    if ( !defined( $self->upload_xref_object_graphs( \@xrefs, $dbi ) ) ) {
        return 1;    # error
    }

    if ( defined $release_file ) {
        # Parse and set release info.
        my $release_io = $self->get_filehandle($release_file);
        local $/ = "\n*";
        my $release = $release_io->getline();
        $release_io->close();

        $release =~ s/\s{2,}/ /g;
        $release =~
	  s/.*(NCBI Reference Sequence.*) Distribution Release Notes.*/$1/s;
        # Put a comma after the release number to make it more readable.
        $release =~ s/Release (\d+)/Release $1,/;

        print "RefSeq release: '$release'\n";

        $self->set_release( $source_id,              $release, $dbi );
        $self->set_release( $peptide_source_id,      $release, $dbi );
        $self->set_release( $mrna_source_id,         $release, $dbi );
        $self->set_release( $ncrna_source_id,        $release, $dbi );
        $self->set_release( $pred_peptide_source_id, $release, $dbi );
        $self->set_release( $pred_mrna_source_id,    $release, $dbi );
        $self->set_release( $pred_ncrna_source_id,   $release, $dbi );
    }

  return 0; # successfull

}

# --------------------------------------------------------------------------------
# Parse file into array of xref objects
# There are 2 types of RefSeq files that we are interested in:
# - protein sequence files *.protein.faa
# - mRNA sequence files *.rna.fna
# Slightly different formats

sub create_xrefs {
  my ($self, $peptide_source_id, $pred_peptide_source_id,
      $mrna_source_id, $ncrna_source_id, 
      $pred_mrna_source_id, $pred_ncrna_source_id, $file, $species_id, $dbi ) = @_;

  # Create a hash of all valid names for this species
  my %species2name = $self->species_id2name($dbi);
  my @names   = @{$species2name{$species_id}};
  my %name2species_id     = map{ $_=>$species_id } @names;
  # my %name2species_id = $self->name2species_id();

  my $refseq_io = $self->get_filehandle($file);

  if ( !defined $refseq_io ) {
      print STDERR "ERROR: Can't open RefSeq file $file\n";
      return;
  }

  my @xrefs;

  local $/ = "\n>";

  while ( $_ = $refseq_io->getline() ) {
    my $xref;

    my $entry = $_;
    chomp $entry;
    my ($header, $sequence) = split (/\n/, $entry, 2);
    $sequence =~ s/^>//;
    # remove newlines
    my @seq_lines = split (/\n/, $sequence);
    $sequence = join("", @seq_lines);

    (my $acc, my $description) = split(/\s/, $header, 2);
    $acc =~ s/^>//;
    my ($species, $mrna);
    if ($file =~ /\.faa(\.gz|\.Z)?$/) {

      ($mrna, $description, $species) = $description =~ /(\S*)\s+(.*)\s+\[(.*)\]$/;
      $xref->{SEQUENCE_TYPE} = 'peptide';
      $xref->{STATUS} = 'experimental';
      my $source_id;
      if ($acc =~ /^XP_/) {
          $source_id = $pred_peptide_source_id;
        } else {
          $source_id = $peptide_source_id;
        }
      $xref->{SOURCE_ID} = $source_id;

    } elsif ($file =~ /\.fna(\.gz|\.Z)?$/) {

      ($species, $description) = $description =~ /\s*(\w+\s+\w+)\s+(.*)$/;
      $xref->{SEQUENCE_TYPE} = 'dna';
      $xref->{STATUS} = 'experimental';
      my $source_id;
      if ($acc =~ /^XM_/ ){
	$source_id = $pred_mrna_source_id;
      } elsif( $acc =~ /^XR/) {
	$source_id = $pred_ncrna_source_id;
      } elsif( $acc =~ /^NM/) {
	$source_id = $mrna_source_id;
      } elsif( $acc =~ /^NR/) {
	$source_id = $ncrna_source_id;
      }
      $xref->{SOURCE_ID} = $source_id;

    }

    if (!$species) { next; }

    $species = lc $species;
    $species =~ s/ /_/;

    my $species_id_check = $name2species_id{$species};

    # skip xrefs for species that aren't in the species table
    if (   defined $species_id
        && defined $species_id_check
        && $species_id == $species_id_check )
    {
      my ($acc_no_ver,$ver) = split (/\./,$acc);
      $xref->{ACCESSION} = $acc_no_ver;
      $xref->{VERSION} = $ver;
      $xref->{LABEL} = $acc;
      $xref->{DESCRIPTION} = $description;
      $xref->{SEQUENCE} = $sequence;
      $xref->{SPECIES_ID} = $species_id;
      $xref->{INFO_TYPE} = "SEQUENCE_MATCH";

      # TODO synonyms, dependent xrefs etc

      push @xrefs, $xref;

    }

  }

  $refseq_io->close();

  return \@xrefs;

}

# --------------------------------------------------------------------------------

1;
