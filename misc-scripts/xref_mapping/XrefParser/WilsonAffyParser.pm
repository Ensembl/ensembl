package XrefParser::WilsonAffyParser;

use strict;

use base qw( XrefParser::BaseParser );

my $xref_sth ;
my $dep_sth;
my $syn_sth;

sub run {

  my $self = shift;
  my $source_id = shift;
  my $species_id = shift;
  my $files       = shift;
  my $release_file   = shift;
  my $verbose       = shift;
#  my ($self, $source_id, $species_id, $file) = @_;

  my @xrefs = $self->create_xrefs($source_id, $species_id, @{$files}[0], $verbose);

  if(!@xrefs){
    return 1; #  1 error
  }
  # upload
  if(!defined(XrefParser::BaseParser->upload_xref_object_graphs(@xrefs))){
    return 1;
  }
  return 0;

}

sub create_xrefs {

  my ($self, $source_id, $species_id, $file, $verbose) = @_;

  my ($count, $noseq, $direct) = (0,0,0);

  $| = 1; # don't buffer

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
      my $xref_id = $self->add_xref($acc, 0, $acc, "$target direct mapping", $source_id, $species_id);
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
