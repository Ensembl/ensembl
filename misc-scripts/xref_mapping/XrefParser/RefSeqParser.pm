# Parse RefSeq files to create xrefs.

package XrefParser::RefSeqParser;

use strict;

use File::Basename;

use base qw( XrefParser::BaseParser );

# --------------------------------------------------------------------------------
# Parse command line and run if being run directly

if (!defined(caller())) {

  if (scalar(@ARGV) != 1) {
    print "\nUsage: RefSeqParser.pm file.SPC <source_id> <species_id>\n\n";
    exit(1);
  }

  run($ARGV[0], -1);

}

# --------------------------------------------------------------------------------

sub run {

  my $self = shift if (defined(caller(1)));

  my $source_id = shift;
  my $species_id = shift;
  my @files = @_;

  my $release_file;

    if ( $files[-1] =~ /RefSeq-release/ ) {
        $release_file = pop @files;
    }

    my $peptide_source_id =
      $self->get_source_id_for_source_name('RefSeq_peptide');
    my $dna_source_id =
      $self->get_source_id_for_source_name('RefSeq_dna');

    print "RefSeq_peptide source ID = $peptide_source_id\n";
    print "RefSeq_dna source ID = $dna_source_id\n";

    my $pred_peptide_source_id =
      $self->get_source_id_for_source_name('RefSeq_peptide_predicted');
    my $pred_dna_source_id =
      $self->get_source_id_for_source_name('RefSeq_dna_predicted');

    print "RefSeq_peptide_predicted source ID = "
      . "$pred_peptide_source_id\n";
    print "RefSeq_dna_predicted source ID = $pred_dna_source_id\n";

    my @xrefs;
    foreach my $file (@files) {
        if ( !defined($species_id) ) {
            $species_id = $self->get_species_id_for_filename($file);
        }

        my $xrefs =
          $self->create_xrefs( $peptide_source_id,
                               $dna_source_id,
                               $pred_peptide_source_id,
                               $pred_dna_source_id,
                               $file,
                               $species_id );

        if ( !defined($xrefs) ) {
            return 1;    #error
        }

        push @xrefs, @{$xrefs};
    }

    if ( !defined( $self->upload_xref_object_graphs( \@xrefs ) ) ) {
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

        $self->set_release( $source_id,              $release );
        $self->set_release( $peptide_source_id,      $release );
        $self->set_release( $dna_source_id,          $release );
        $self->set_release( $pred_peptide_source_id, $release );
        $self->set_release( $pred_dna_source_id,     $release );
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
  my $self = shift;

  my ( $peptide_source_id, $dna_source_id, $pred_peptide_source_id,
      $pred_dna_source_id, $file, $species_id ) = @_;

  # Create a hash of all valid names for this species
  my %species2name = $self->species_id2name();
  my @names   = @{$species2name{$species_id}};
  my %name2species_id     = map{ $_=>$species_id } @names;
  # my %name2species_id = $self->name2species_id();

  my $refseq_io = $self->get_filehandle($file);

  if ( !defined $refseq_io ) {
      print "ERROR: Can't open RefSeq file $file\n";
      return undef;
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

    (my $gi, my $n, my $ref, my $acc, my $description) = split(/\|/, $header);
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
      if ($acc =~ /^XM_/ || $acc =~ /^XR/) {
	$source_id = $pred_dna_source_id;
      } else {
	$source_id = $dna_source_id;
      }
      $xref->{SOURCE_ID} = $source_id;

    }

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

      # TODO synonyms, dependent xrefs etc

      push @xrefs, $xref;

    }

  }

  $refseq_io->close();

  print "Read " . scalar(@xrefs) ." xrefs from $file\n";

  return \@xrefs;

}

# --------------------------------------------------------------------------------

1;
