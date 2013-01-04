package XrefParser::VBCommunitySynonymParser;
 
use strict;
use warnings;
use Carp;
use POSIX qw(strftime);
use File::Basename;
 
use base qw( XrefParser::BaseParser );


sub run {
  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $files        = $ref_arg->{files};
  my $release_file = $ref_arg->{rel_file};
  my $verbose      = $ref_arg->{verbose};

  if ((!defined $source_id) or (!defined $species_id) or (!defined $files)) {
    croak "Need to pass source_id, species_id, files and rel_file as pairs";
  }
  $verbose |=0;

  my $file = @{$files}[0];
  #The Synonyms are for gene symboles. Therefore using the source_id of symboels
  my $symboel_source_id;
  open IN,"<","source_id_file" or die "could not open source_id_file";
  while(my $line = <IN>){
  	  chomp $line;
  	  my($organism_id,$source_id) = split /\t/ , $line;
  	  if($organism_id eq $species_id){
  	  	$symboel_source_id = $source_id;
  	  }
  }
  close IN;
  #my $symboel_source_id = $self->get_source_id_for_source_name("VB_Community_Annotation");
  print "source_id = $symboel_source_id, species= $species_id, file = $file\n" if $verbose;

  my $added = 0;
  my $count = 0;

  my $file_io = $self->get_filehandle($file);

  if ( !defined $file_io ) {
    print STDERR "ERROR: Could not open file $file\n";

    return 1;
  }

  while ( my $line = $file_io->getline() ) {
    next unless $line =~ /^\w+/;

    chomp $line;
    my ($stable_id,$synonym) = split "\t", $line;
    $count++;

    
    #Add synonym for an xref given by accession and source_id
    $self->add_to_syn($stable_id, $symboel_source_id,$synonym,$species_id);

  }

  $file_io->close();

  print "Added $count synonyms to synonym  for VBCommunitySynonyms\n" if $verbose;

  return 0;
}

1;
