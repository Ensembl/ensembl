package XrefParser::CGNCParser;

use strict;
use warnings;
use Carp;
use DBI;

use base qw(XrefParser::BaseParser);

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

  my $file = @{$files}[0];

  my $file_io = $self->get_filehandle($file);
  if ( !defined $file_io ) {
    print STDERR "ERROR: Could not open $file\n";
    return 1;    # 1 is an error
  }

  $source_id = $self->get_source_id_for_source_name("CGNC");

  
  my $count = 0;
  while ( $_ = $file_io->getline() ) {
#48941  ENSGALG00000002652      FZD10

    chomp;
    my @array = split /\t/x, $_;

    my $ensid = $array[1];
    my $acc = $array[2];
    my $desc = $array[3];

    if($ensid =~ /ENSGAL/){
      my $xref_id = $self->add_xref({ acc        => $acc,
				      version    => 0,
				      label      => $acc,
				      desc       => $desc,
				      source_id  => $source_id,
				      species_id => $species_id,
				      info_type  => "DIRECT"} );

      $self->add_direct_xref( $xref_id, $ensid, "Gene", '');
      $count++;
    } else{
      print STDERR "No match for $acc\n";
    }
  }
  print "$count direct CGNC xrefs added\n";
  return 0;

}

1;
