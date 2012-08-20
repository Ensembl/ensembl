package XrefParser::RFAMParser;

use strict;
use warnings;
use Carp;

use base qw( XrefParser::BaseParser );
use Bio::EnsEMBL::Registry;

sub run {
  my ( $self, $ref_arg ) = @_;
  my $source_id  = $ref_arg->{source_id};
  my $species_id = $ref_arg->{species_id};
  my $files      = $ref_arg->{files};
  my $verbose    = $ref_arg->{verbose};
  my $mapper     = $ref_arg->{mapper};

  if ( ( !defined $source_id )
    or ( !defined $species_id )
    or ( !defined $files ) )
  {
    croak "Need to pass source_id, species_id and file as pairs";
  }
  $verbose |= 0;

  my $rfam_transcript_stable_ids =
    $self->get_rfam_transcript_stable_ids($mapper);

  my $file    = @{$files}[0];
  my $file_io = $self->get_filehandle($file);
  if ( !defined $file_io ) {
    print STDERR "ERROR: Could not open $file\n";
    return 1;    # 1 is an error
  }

  my @xrefs;

  local $/ = "//\n";

  my $xref_count;
  my $direct_count;

  while ( $_ = $file_io->getline() ) {

    my $xref;

    my $entry = $_;
    chomp $entry;

    next if ( !$entry );

    my ($accession)   = $entry =~ /\n#=GF\sAC\s+(\w+)/;
    my ($label)       = $entry =~ /\n#=GF\sID\s+([^\n]+)/;
    my ($description) = $entry =~ /\n#=GF\sDE\s+([^\n]+)/;

    if ( exists( $rfam_transcript_stable_ids->{$accession} ) ) {

      #add xref
      my $xref_id = $self->add_xref(
        {
          acc        => $accession,
          version    => 0,
          label      => $label || $accession,
          desc       => $description,
          source_id  => $source_id,
          species_id => $species_id,
          info_type  => "DIRECT"
        }
      );

      my @transcript_stable_ids =
        @{ $rfam_transcript_stable_ids->{$accession} };

      foreach my $stable_id (@transcript_stable_ids) {
        $self->add_direct_xref( $xref_id, $stable_id, "Transcript", "" );
        $direct_count++;
      }

      $xref_count++;

    }

  }

  $file_io->close();

  print "Added $xref_count RFAM xrefs and $direct_count direct xrefs\n"
    if ($verbose);
  if ( !$xref_count ) {
    return 1;    # 1 error
  }

  return 0;      # successfull
}

sub get_rfam_transcript_stable_ids {
  my ( $self, $mapper ) = @_;
  return {};
}

1;
