
=pod

=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 NAME

GFFParser - simple gff3 parser.


=head1 AUTHOR

Monika Komorowska, 2012 - monika@ebi.ac.uk

=head1 SYNOPSIS

use strict;
use Bio::EnsEMBL::Utils::IO::GFFParser;
use Bio::EnsEMBL::Utils::Scalar qw/wrap_array/;
use FileHandle;

my $file_name = "features.gff";
my $fh = FileHandle->new;
$fh->open("< $file_name");
my $parser = Bio::EnsEMBL::Utils::IO::GFFParser->new($fh);

my @header_lines = @{$parser->parse_header()};
#do something with the header lines array, e.g. print array elements

foreach my $header_line (@header_lines) {
    print $header_line . "\n";
}
print "\n\n";
my $feature = $parser->parse_next_feature();

while (defined($feature) ) {

    my %feature = %{$feature};

    #do something with the feature, e.g. print hash keys and values 
    foreach my $key (keys %feature) {
	if ($key ne 'attribute') {
	    print $key . " " . $feature{$key} ."\n";
	} else {
	    print $key . "\n";
	    my %attribs =  %{$feature{$key}};
	    foreach my $attrib_key (keys %attribs) {
	      my $values = $attribs{$attrib_key};
		    printf("\t%s %s\n", $attrib_key, join(q{, }, wrap_array($values)));
	    }
	}
    }
    print "\n\n";
    $feature = $parser->parse_next_feature();
}

$parser->close();

$fh->close();



=head1 DESCRIPTION

GFF3 format as defined in http://www.sequenceontology.org/gff3.shtml.

Use parse_header method to parse a GFF3 file header, and parse_next_feature to parse the next feature line in the file.

This class can be extended to convert a feature hash into a feature object reversing
the processing done by GFFSerializer.

=cut

package Bio::EnsEMBL::Utils::IO::GFFParser;
use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception;
use IO::File;
use URI::Escape;

my %strand_conversion = ( '+' => '1', '?' => '0', '-' => '-1' );

=head2 new

    Constructor
    Arg [1]    : File handle
    
    Returntype : Bio::EnsEMBL::Utils::IO::GFFParser

=cut

sub new {
  my $class = shift;
  my $self = { filehandle => shift, };
  bless $self, $class;
  if ( !defined( $self->{'filehandle'} ) ) {
    throw("GFFParser requires a valid filehandle to a GFF3 formatted file");
  }
  return $self;

}

=head2 parse_header

    Arg [1]    : File handle 
    Description: Returns a arrayref with each header line stored in array element
    Returntype : Arrayref of GFF3 file header lines

=cut

sub parse_header {

  my $self = shift;

  my $next_line;
  my @header_lines;

  while ( ( $next_line = $self->_read_line() ) && ( $next_line =~ /^[\#|\s]/ ) )
  {

    #header lines start with ##
    if ( $next_line =~ /^[\#]{2}/ ) {
      push @header_lines, $next_line;
      if ( $next_line =~ /gff-version\s+(\d+)/ ) {
        if ( $1 != 3 ) {
          warning(
"File has been formatted in GFF version $1. GFFParser may return unexpected results as it is designed to parse GFF3 formatted files."
          );
        }
      }
    }
  }

  if ( defined($next_line) && ( $next_line !~ /^[\#|\s]/ ) ) {
    $self->{'first_feature_line'} = $next_line;
  }
  return \@header_lines;

}

=head2 parse_next_feature

    Arg [1]    : File handle
    Description: Returns a hashref in the format -
                 {
                   seqid => scalar,
                   source => scalar,
                   type => scalar,
                   start => scalar,
                   end => scalar,
                   score => scalar,
                   strand => scalar,
                   phase => scalar,
                   attribute => hashref, 
                   
		 }
		If the attribute value held more than one value then we hold an arrayref
		not a scalar
    Returntype : Hashref of a GFF3 feature line

=cut

sub parse_next_feature {

  my $self = shift;

  #    my $next_line;
  my $feature_line;

  while ( my ($next_line) = $self->_read_line() ) {
    next
      if ( $next_line =~ /^\#/
      || $next_line =~ /^\s*$/
      || $next_line =~ /^\/\// );
    $feature_line = $next_line;
    last;
  }

  return undef unless $feature_line;

  my %feature;
  my %attribute;

  #strip off trailing comments
  $feature_line =~ s/\#.*//;

  my @chunks = split( /\t/, $feature_line );

  %feature = (
    'seqid'  => uri_unescape( $chunks[0] ),
    'source' => uri_unescape( $chunks[1] ),
    'type'   => uri_unescape( $chunks[2] ),
    'start'  => $chunks[3],
    'end'    => $chunks[4],
    'score'  => $chunks[5],
    'strand' => $strand_conversion{ $chunks[6] },
    'phase'  => $chunks[7]
  );

  if ( $chunks[8] ) {
    my @attributes = split( /;/, $chunks[8] );
    my %attributes;
    foreach my $attribute (@attributes) {
      my ( $name, $value ) = split( /=/, $attribute );
      $name = uri_unescape($name);
      my @split_values = map { uri_unescape($_) } split(/,/, $value);
      if(scalar(@split_values) > 1) {
        $attributes{$name} = \@split_values;
      }
      else {
        $attributes{$name} = $split_values[0];
      }
    }
    $feature{'attribute'} = \%attributes;
  }

  return \%feature;
}

sub _read_line {

  my $self = shift;
  my $fh   = $self->{'filehandle'};

  my $line;

  if ( defined( $self->{'first_feature_line'} ) ) {
    $line = $self->{'first_feature_line'};
    $self->{'first_feature_line'} = undef;
  }
  else {
    $line = <$fh>;
    if ( defined($line) ) {
      chomp $line;
    }
  }

  return $line;
}

sub close {
  my $self = shift;
  $self->{"filehandle"} = undef;
}

1;
