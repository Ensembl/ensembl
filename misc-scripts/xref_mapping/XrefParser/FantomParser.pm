package XrefParser::FantomParser;

use strict;
use warnings;
use Carp;
use POSIX qw(strftime);
use File::Basename;
use File::Spec::Functions;

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

  my $file = @{$files}[0];

  my $dir = dirname($file);


  my (%embl) = %{$self->get_valid_codes("embl",$species_id)};

  my $fantom_io =
    $self->get_filehandle( $file  );

  if ( !defined $fantom_io ) {
    print STDERR "ERROR: Could not open " . $file . "\n" ;
    return 1;    # 1 error
  }

  my $ecount =0;

  my $mismatch=0;

  $fantom_io->getline(); # remove header

  while ( $_ = $fantom_io->getline() ) {
    chomp;
    my ($master, $label, $acc) = split (/\s+/,$_);
    if(defined($embl{$master})){
      foreach my $xref_id (@{$embl{$master}}){
	$self->add_dependent_xref({ master_xref_id => $xref_id,
				    acc            => $label,
				    label          => $label,
				    source_id      => $source_id,
				    species_id     => $species_id} );
	$ecount++;
      }
    }
    else{
      if($mismatch < 10){
	print "Could not find master $master\n";
      }
      $mismatch++;
    }
  }

  $fantom_io->close();

  if($verbose){
    print "\t$ecount xrefs from EMBL and\n";
    print "\t$mismatch xrefs ignored as no master found\n";
  }
  return 0;
}

1;
