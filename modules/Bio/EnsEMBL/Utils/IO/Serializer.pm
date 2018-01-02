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

Serializer - An abstract serializer for turning EnsEMBL data into other formats

=head1 SYNOPSIS

my $serializer = new Serializer( $filehandle );
$serializer->print_feature_list( \@list_of_features );

=head1 DESCRIPTION

Subclass this class to create a format-specific serializer.
Be sure to implement print_feature at the bare minimum

=cut

package Bio::EnsEMBL::Utils::IO::Serializer;
use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception;
use Bio::EnsEMBL::Utils::SeqDumper;


=head2 new

	Constructor
	Arg [1]    : Optional File handle
	Returntype : Bio::EnsEMBL::Utils::IO::Serializer

=cut

sub new {
	my $class = shift;
	my $self = {
		'filehandle' => shift,
		'achieved_something' => 0,
	};
	bless $self, $class;
	if (!defined ($self->{'filehandle'})) {
		# no file handle, let the handle point to a copy of STDOUT instead
		open $self->{'filehandle'}, ">&STDOUT";
		$self->{'stdout'} = 1;
	}
	return $self;
}

=head2 DESTROY

	Destructor
	Description: Restores default state of the STDOUT filehandle as it is a copy
	             and may not flush correctly.
=cut

sub DESTROY {
	my $self = shift;
	if ($self->{'stdout'}) {
		close $self->{'filehandle'};
	}
}

=head2 print_metadata 
	
	Arg [1]    : String
	Description: Pipes a custom string into the filehandle that the serializer is using
        Returntype:  None
=cut

sub print_metadata {
	my $self = shift;
	my $text = shift;
	my $fh = $self->{'filehandle'};
	print $fh "\n".$text."\n";
}

=head2 print_main_header

	Arg [1]    : Data for header, depends on serializer
	Description: Abstract method for printing the header text or metadata required for this file format,
	             Re-implement in the serializer.
	Returntype : None
=cut

sub print_main_header {
	my $self = shift;
	warning("No writer for headers in this format. Nothing done" );
}

=head2 printed_something
    Description: Check if serializer has printed any useful data. Not accurate with FASTA
                 due to non-reporting dumper.
    Returntype : Boolean
=cut

sub printed_something {
	my $self = shift;
	if ($self->{'achieved_something'}) { return 1;}
	else {return 0;}
}

=head2 formatted_write

   Arg [1]    : Line format, see Perldoc of formline()
   Arg [2]    : Array of arguments to suit the line format in Arg [1]
   Description: Writes data to the filehandle and rigidly formats it.
                Refer to Perldoc on formline() to specify valid formats.
                Useful for fixed-width file formats.
                Suicides in the event of file system issues.
   Example    : my $FORMAT = '^<<<<<<<<<<<<<<<<<<<|<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n';
                $serializer->formatted_write($FORMAT,@text_fields);
   Returntype : None
=cut

sub formatted_write {
  my ($self, $FORMAT, @values) = @_;
  my $fh = $self->{'filehandle'};
  
  #while the last value still contains something
  while(defined($values[-1]) and $values[-1] ne '') {
    formline($FORMAT, @values);
    print $fh $^A or die "Failed write to filehandle";
    $^A = '';
  }
}

1;
