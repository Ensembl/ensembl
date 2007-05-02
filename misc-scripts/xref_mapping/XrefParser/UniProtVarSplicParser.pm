package XrefParser::UniProtVarSplicParser;

# Parse UniProt alternative splice files

use strict;
use File::Basename;

use base qw( XrefParser::BaseParser );

# UniProtVarSplic file format: fasta, e.g.

#>P48347-2|14310_ARATH Isoform 2 of P48347 - Arabidopsis thaliana (Mouse-ear cress)
#MENEREKQVYLAKLSEQTERYDEMVEAMKKVAQLDVELTVEERNLVSVGYKNVIGARRAS
#WRILSSIEQKEESKGNDENVKRLKNYRKRVEDELAKVCNDILSVIDKHLIPSSNAVESTV
#FFYKMKGDYYRYLAEFSSGAERKEAADQSLEAYKAAVAAAENGLAPTHPVRLGLALNFSV
#FYYEILNSPESACQLAKQAFDDAIAELDSLNEESYKDSTLIMQLLRDNLTLWTSDLNEEG
#DERTKGADEPQDEV

sub run {

  my ( $self, $source_id, $species_id, $file, $release_file ) = @_;

  my @xrefs;

  local $/ = "\n>";

  my $file_io = $self->get_filehandle($file);

  if ( !defined $file_io ) {
    print "ERROR: Could not open $file\n";
    return 1;    # 1 error
  }

  my %swiss = %{ $self->get_valid_codes( "uniprot", $species_id ) };

  my $missed = 0;
  while ( $_ = $file_io->getline() ) {
    my $xref;

    my ($header, $sequence) = $_ =~ /^>?(.+?)\n([^>]*)/s or warn("Can't parse FASTA entry: $_\n");

    # deconstruct header
    my ($accession, @description) = split /\|/, $header;
    my $description = join(" ", @description);
    
    my ($original, $extension) = split/-/, $accession;

    if(defined($swiss{$original})){
      # make sequence into one long string
      $sequence =~ s/\n//g;
      
      # build the xref object and store it
      $xref->{ACCESSION}     = $accession;
      $xref->{LABEL}         = $accession;
      $xref->{DESCRIPTION}   = $description;
      $xref->{SEQUENCE}      = $sequence;
      $xref->{SOURCE_ID}     = $source_id;
      $xref->{SPECIES_ID}    = $species_id;
      $xref->{SEQUENCE_TYPE} = 'peptide';
      $xref->{STATUS}        = 'experimental';
      
      push @xrefs, $xref;
    }
    else{
     $missed++;
    }
  }

  $file_io->close();

  print $missed." ignored as original uniprot not found in database\n";
  print scalar(@xrefs) . " UniProtVarSplic xrefs succesfully parsed\n";

  $self->upload_xref_object_graphs(\@xrefs);

    if ( defined $release_file ) {
        # Parse and apply the Swiss-Prot release info
        # from $release_file.
        my $release_io = $self->get_filehandle($release_file);
        while ( defined( my $line = $release_io->getline() ) ) {
            if ( $line =~ m#(UniProtKB/Swiss-Prot Release .*)# ) {
                print "Swiss-Prot release is '$1'\n";
                $self->set_release( $source_id, $1 );
            }
        }
        $release_io->close();
    }


  print "Done\n";
  return 0;
}

1;
