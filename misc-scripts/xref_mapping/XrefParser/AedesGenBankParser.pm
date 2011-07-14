package XrefParser::AedesGenBankParser;

use strict;
use File::Basename;

use base qw( XrefParser::BaseParser );

#Aedes GenBank protein - because not yet in UniProt
#>EAT48991.1
#MGKSKAHRIKGLTGPKMSLGDQITEGRVSKKPKAPKIRLRAEEEEFVDSRTTKKILQQAR
#KQQAELNLLDDSFGPSLAESAAAASVGKRRHRLGDAASSDESDEEYREEADVDGQDFFDD
#IKINEEDERALEMFQNKDGVKTRTLADLIMDKITEKQTEIQTQFSDTGSLKMEEIDPRVR

sub run {

  my $self = shift;
  my $source_id = shift;
  my $species_id = shift;
  my $files       = shift;
  my $release_file   = shift;
  my $verbose       = shift;

  my $file = @{$files}[0];

  my $cpt = 0 ;

  next if (/^File:/);   # skip header

  my @xrefs;

  local $/ = "\n>";

  my $file_io = $self->get_filehandle($file);
  if ( !defined $file_io ) {
      print "Could not open $file\n";
      return 1;
  }

  while ( $_ = $file_io->getline() ) {

    my $xref;

    my ($header, $sequence) = $_ =~ /^>?(.+?)\n([^>]*)/s or warn("Can't parse FASTA entry: $_\n");
    #print "My header is -$header-\n" ;
    #print "My sequence is -$sequence-\n" ;

    if ($header eq "") {
      $header = "Aedes_GenBank".$cpt ;
      print STDERR "One sequence with a random name ... \n" if($verbose);
      $cpt++ ;
    }

    # deconstruct header - just use first part
    #my ($accession, $description) = split /\|/, $header;  #if description
    my $accession = $header;                               #if no description



    # make sequence into one long string
    $sequence =~ s/\n//g;

    # build the xref object and store it
    #print "ACCESSION & LABEL are $accession\n" ;
    #print "SEQUENCE is $sequence\n" ;
    #print "SOURCE_ID is $source_id\n" ;
    #print "SPECIES_ID is $species_id\n" ;
    #print "SEQUENCE_TYPE is peptide!\n";
    #print "STATUS is experimental!\n" ;

    $xref->{ACCESSION}     = $accession;
    $xref->{LABEL}         = $accession;
    #$xref->{DESCRIPTION}   = $description;
    $xref->{SEQUENCE}      = $sequence;
    $xref->{SOURCE_ID}     = $source_id;
    $xref->{SPECIES_ID}    = $species_id;
    $xref->{SEQUENCE_TYPE} = 'peptide';
    $xref->{STATUS}        = 'experimental';

    push @xrefs, $xref;

  }

  $file_io->close();

  print scalar(@xrefs) . " AedesGenBank xrefs succesfully parsed\n" if($verbose);

  XrefParser::BaseParser->upload_xref_object_graphs(\@xrefs);

  print "Done\n";
  return 0;
}

1;
