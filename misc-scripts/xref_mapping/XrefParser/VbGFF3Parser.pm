package XrefParser::VbGFF3Parser;

use strict;
use warnings;
use Carp;
use File::Basename;
use Bio::SeqIO;

# perl xref_parser.pl -user ensadmin -pass PASS -host genebuild2 -dbname snr_anopheles_gambiae_48_xref_test -species anopheles_gambiae -source ARRAY_JHSPH_AGGAMBER_15k_v1 -download_path $mywork/VB_xref/downloads/ -checkdownload


use base qw( XrefParser::CoordinateParser );

# Parser for GFF3-format probe mappings from Vectorbase
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

  open (my $FH, "<", $file) || return 1;
  my $i=0; my $type = "transcript";
  while (my $ln = <$FH>) {
#	# parse GFF line:
    chomp($ln); $i++;
    my ($seqid, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split("\t",$ln);
    if ($strand eq "."){
      $strand = 0;
    }
    elsif ($strand eq "?"){
      $strand = 0;
    }# nb should more properly be NULL;
    my %attributes;
    unless ($attributes eq ".") {
      my @attributes = split(";",$attributes);
      foreach my $pair ( @attributes ){
	my ($key, $value) = split("=",$pair);
	
	if ($value =~ m/,/i ) {
	  my @values = split(",",$value);
	  foreach my $splitval (@values) {
	    $splitval =~ s/%59/;/;
	    $splitval =~ s/%44/,/;
	    $splitval =~ s/%61/=/;
	    push (@{$attributes{uc($key)}},$splitval);
	  }
	}
	else {
	  $value =~ s/%59/;/;
	  $value =~ s/%44/,/;
	  $value =~ s/%61/=/;
	  push(@{$attributes{uc($key)}},$value);
	}
      }
    }


    #		# assess and build oligo_feature
    if (($seqid eq ".") || ($start eq ".") || ($end eq ".")) {
      #				_parse_error(\%attributes,"incomplete location (".$seqid.":".$start."-".$end.")");
    }
    else {
      foreach my $name (@{$attributes{"NAME"}}) {
	my %xref = ( 'accession'  => $name,
		     'chromosome' => $seqid,
		     'strand'     => $strand,
		     'txStart'    => $start,
		     'txEnd'      => $end,
		     'cdsStart'   => $start,
		     'cdsEnd'     => $end,
		     'exonStarts' => $start,
		     'exonEnds'   => $end );
	print STDERR "$name, $start, $end\n" if($verbose);
	$self->add_xref( $source_id, $species_id, \%xref );
      }
    }
  }
  print STDERR $i." VB GFF3 xrefs succesfully parsed\n" if($verbose);
  close($FH);

  return 0;
} ## end sub run

1;
