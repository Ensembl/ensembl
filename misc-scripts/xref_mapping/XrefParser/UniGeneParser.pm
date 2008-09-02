# Parse UniGene Hs.seq.uniq files to create xrefs.

package XrefParser::UniGeneParser;

use strict;

use File::Basename;

use base qw( XrefParser::BaseParser );

my $verbose;

# --------------------------------------------------------------------------------
# Parse command line and run if being run directly

if (!defined(caller())) {

  if (scalar(@ARGV) != 1) {
    print "\nUsage: UniGeneParser.pm file.SPC <source_id> <species_id>\n\n";
    exit(1);
  }

  run($ARGV[0], -1);

}

# --------------------------------------------------------------------------------

sub run {

  my $self = shift if (defined(caller(1)));

  my $source_id = shift;
  my $species_id = shift;
  my $files       = shift;
  my $release_file   = shift;
  $verbose       = shift;

  my $uniq_file = @{$files}[0];
  my $data_file = @{$files}[1];

  my $unigene_source_id = $self->get_source_id_for_source_name('UniGene');

  print "UniGene source ID = $unigene_source_id.\n" if($verbose);

  if ( !defined($species_id) ) {
    $species_id =
      $self->get_species_id_for_filename($uniq_file);
  }

  my $xrefs =
    $self->create_xrefs( $unigene_source_id, $unigene_source_id,
      $uniq_file, $data_file, $species_id );

  if(!defined($xrefs)){
    return 1; #error
  }
  if(!defined($self->upload_xref_object_graphs($xrefs))){
    return 1; # error
  }

    if ( defined $release_file ) {
        # Get species name from species ID.
        my $species_name;

        my $sth =
          $self->dbi()
          ->prepare("SELECT name FROM species WHERE species_id = ?");

        $sth->execute($species_id);
        $sth->bind_columns( \$species_name );
        $sth->fetchrow_array();

        $species_name =~ tr/_/ /;

        # Parse and set release info.
        my $release;
        my $release_io = $self->get_filehandle($release_file);

        while ( defined( my $line = $release_io->getline() ) ) {
            if ( $line =~ /^(.*$species_name)/i ) {
                $release = $1;
            }
        }
        $release_io->close();

        if ( defined $release ) {
            $release =~ s/\s{2,}/ /g;
            $release =~ s/^(.*) UniGene/$1, UniGene/;

            print "UniGene release: '$release'\n" if($verbose);
            $self->set_release( $unigene_source_id, $release );
        }
    }

  return 0; # successfull

}


my %geneid_2_desc;

sub get_desc{
  my $self = shift;
  my $data_file = shift;

  my $dir = dirname($data_file);

  local $/ = "//";

  my $desc_io = $self->get_filehandle( $data_file );

  if ( !defined $desc_io ) {
    print STDERR "ERROR: Can't open $data_file\n";
    return undef;
  }

  while ( $_ = $desc_io->getline() ) {
    #ID          Hs.159356
    #TITLE       Hypothetical LOC388277
    
    (my $id) = $_ =~ /ID\s+(\S+)/;
    (my $descrip) = $_ =~ /TITLE\s+(.+)\n/;

    if ( defined $id && defined $descrip ) {
        $geneid_2_desc{$id} = $descrip;
    }
    
  }

  $desc_io->close();

  return 1;
}


sub create_xrefs {
  my $self = shift;

  my ( $peptide_source_id, $unigene_source_id, $uniq_file, $data_file,
      $species_id )
    = @_;

  # Create a hash of all valid names for this species. Not used...
  # my %species2name = $self->species_id2name();
  # my @names   = @{$species2name{$species_id}};
  # my %name2species_id     = map{ $_=>$species_id } @names;

  if ( !defined( $self->get_desc($data_file) ) ) {
    return undef;
  }

  my $unigene_io = $self->get_filehandle($uniq_file);

  if ( !defined $unigene_io ) {
    print STDERR "Can't open RefSeq file $uniq_file\n";
    return undef;
  }

#>gnl|UG|Hs#S19185843 Homo sapiens N-acetyltransferase 2 (arylamine N-acetyltransferase)
  # , mRNA (cDNA clone MGC:71963 IMAGE:4722596), complete cds /cds=(105,977) /gb=BC067218 /gi=45501306 /ug=Hs.2 /len=1344
#GGGGACTTCCCTTGCAGACTTTGGAAGGGAGAGCACTTTATTACAGACCTTGGAAGCAAG


  my @xrefs;

  local $/ = "\n>";

  while ( $_ = $unigene_io->getline() ) {

    my $xref;

    my $entry = $_;
    chomp $entry;
    my ($header, $sequence) = split (/\n/, $entry, 2);
    $sequence =~ s/^>//;
    # remove newlines
    my @seq_lines = split (/\n/, $sequence);
    $sequence = join("", @seq_lines);

#    (my $gnl, my $n, my $rest) = split(/\|/, $header,3);

    (my $acc_no_ver) = $header =~ /\/ug=(\S*)/;

    if(!defined($geneid_2_desc{$acc_no_ver})){
      print "****$_\n";
      $geneid_2_desc{$acc_no_ver} = "";
      warn "No desc for $acc_no_ver\n";
    }


    $xref->{SEQUENCE_TYPE} = 'dna';
    $xref->{STATUS} = 'experimental';
    $xref->{SOURCE_ID} = $unigene_source_id;
   

    ##No species check as files contain data  fro only one species.
    
    $xref->{ACCESSION} = $acc_no_ver;
    $xref->{LABEL} = $acc_no_ver;
    $xref->{DESCRIPTION} = $geneid_2_desc{$acc_no_ver};
    $xref->{SEQUENCE} = $sequence;
    $xref->{SPECIES_ID} = $species_id;
    
    push @xrefs, $xref;
    
  }

  $unigene_io->close();

  %geneid_2_desc=();
  print "Read " . scalar(@xrefs) ." xrefs from $uniq_file\n" if($verbose);

  return \@xrefs;

}

1;
