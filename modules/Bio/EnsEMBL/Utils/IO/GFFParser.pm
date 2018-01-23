=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

=pod


=head1 NAME

GFFParser - simple gff3 parser.


=head1 SYNOPSIS

use strict;
use Bio::EnsEMBL::Utils::IO::GFFParser;
use IO::File;

my $file_name = "features.gff";
my $fh = IO::File->new($file_name, 'r');
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
        printf("\t%s %s\n", $attrib_key, join(q{, }, @{wrap_array($values)}));

        }
    }
    }
    print "\n\n";
    $feature = $parser->parse_next_feature();
}

my $sequence = $parser->parse_next_sequence();

while (defined($sequence)) {
    my %sequence = %{$sequence};

    foreach my $key (keys %sequence) {      
        print $key . " " . $sequence{$key} ."\n";
    }
    print "\n\n";   

    $sequence = $parser->parse_next_sequence();
}

$parser->close();

$fh->close();



=head1 DESCRIPTION

GFF3 format as defined in http://www.sequenceontology.org/gff3.shtml

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
use Bio::EnsEMBL::Utils::Scalar qw/wrap_array/;


my %strand_conversion = ( '+' => '1', '?' => '0', '-' => '-1');

=head2 new

    Constructor
    Arg [1]    : File handle
    
    Returntype : Bio::EnsEMBL::Utils::IO::GFFParser

=cut

sub new {
    my $class = shift;
    my $self = {
        filehandle => shift,
    };
    bless $self, $class;
    if (!defined($self->{'filehandle'})) {
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
    
    while (($next_line = $self->_read_line()) && ($next_line =~ /^[\#|\s]/) )  {
    
        #stop parsing features if ##FASTA directive encountered
        last if ($next_line =~ /\#\#FASTA/ );
    
        #header lines start with ## (except for the ##FASTA directive indicating sequence section)
        if ($next_line =~ /^[\#]{2}/ ) {
            push @header_lines, $next_line;
            if ($next_line =~ /gff-version\s+(\d+)/) {
                if ($1 != 3) {
                    warning("File has been formatted in GFF version $1. GFFParser may return unexpected results as it is designed to parse GFF3 formatted files.");  
                }
            }
        }
    }

    if (defined($next_line)) {
        $self->{'first_non_header_line'} = $next_line;
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
    Returntype : Hashref of a GFF3 feature line

=cut

sub parse_next_feature {

    my $self = shift;

    my $next_line;
    my $feature_line;
    
    while (($next_line = $self->_read_line() ) && defined($next_line) ) {

        #stop parsing features if ##FASTA directive
        last if ($next_line =~ /\#\#FASTA/);
    
    
        next if ($next_line =~ /^\#/ || $next_line =~ /^\s*$/ ||
            $next_line =~ /^\/\//);
    
        $feature_line = $next_line;
        last;
    }

    return undef unless $feature_line;

    my %feature;
    my %attribute;


    #strip off trailing comments
    $feature_line =~ s/\#.*//;
    
    my @chunks = split(/\t/, $feature_line);

    %feature = (
        'seqid' => uri_unescape($chunks[0]),
        'source' => uri_unescape($chunks[1]),
        'type' => uri_unescape($chunks[2]),
        'start' => $chunks[3],
        'end' => $chunks[4],
        'score' => $chunks[5],
        'strand' => $strand_conversion{$chunks[6]},
        'phase' => $chunks[7] 
    );
    
    if ($chunks[8]) {
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

=head2 parse_next_sequence

    Arg [1]    : File handle
    Description: Returns a hashref in the format -
                 {
                   header => scalar,
                   sequence => scalar,
                   
         }
    Returntype : Hashref of a GFF3 sequence line

=cut

sub parse_next_sequence {

    my $self = shift;

    my $next_line;
    my $sequence;
    my $header;
    
    while (($next_line = $self->_read_line() ) && defined($next_line) ) {

        next if ($next_line =~ /^\#/ || $next_line =~ /^\s*$/ ||
            $next_line =~ /^\/\//);
    
        if ($next_line =~ /^>/) {
            if ($header) {
                #next fasta header encountered
                $self->{'next_fasta_header'} = $next_line; 
                last;
            
            } else {
                $header = $next_line;
            }
        } else {
            $sequence .= $next_line;
        }
    }

    return undef unless ($sequence || $header);

    my %sequence = (header => $header , sequence => $sequence );

    return \%sequence;    
}


sub _read_line {

    my $self = shift;
    my $fh = $self->{'filehandle'};

    my $line;
    
    if (defined($self->{'first_non_header_line'})) {
        $line = $self->{'first_non_header_line'};
        $self->{'first_non_header_line'} = undef;
    } elsif ( defined($self->{'next_fasta_header'} )) {
        $line = $self->{'next_fasta_header'};
        $self->{'next_fasta_header'} = undef;
    }
    else {
        $line = <$fh>;
        if (defined($line)) {
            chomp $line;
            if (!$line) {
            #parse next line if current line is empty
            $line = $self->_read_line();
            }
        }
    }

    return $line;
}

sub close {

    my $self = shift;
    $self->{"filehandle"} = undef;

}

1;
