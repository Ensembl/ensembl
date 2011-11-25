=pod

=head1 LICENSE

  Copyright (c) 1999-2011 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

  http://www.ensembl.org/info/about/code_licence.html

=head1 NAME

Serializer - An abstract serializer for turning EnsEMBL data into other formats

=head1 AUTHOR

Kieron Taylor, 2011 - ktaylor@ebi.ac.uk

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

sub print_feature {
	throw( "print_feature method not implemented.");
}

=head2 print_feature_list

	Arg [1]    : Listref of features
	Description: Run print_feature on every feature in the list

=cut

sub print_feature_list {
	my $self = shift;
	my $feature_list = shift;
	if (ref($feature_list) eq 'ARRAY') {
		if (scalar(@$feature_list) > 0) {$self->{'achieved_something'} = 1;} 
		foreach my $feature (@{$feature_list}) {
			$self->print_feature($feature);
		}
	}
	else {
		throw( "print_feature_list requires a listref as argument" );
	}
}

=head2 print_feature_Iterator

	Arg [1]    : Bio::EnsEMBL::Utils::Iterator
	Description: Automatically spools through an iterator for convenience
	Returntype : None
=cut

sub print_feature_Iterator {
	my $self = shift;
	my $iterator = shift;
	if ($iterator->can('has_next')) {
		$iterator->each(sub {$self->print_feature($_); $self->{'achieved_something'} = 1;});
	}
	else {
		throw("Supplied iterator does not look like Bio::EnsEMBL::Utils::Iterator");
	}
}

=head2 print_metadata 
	
	Arg [1]    : String
	Description: Pipes a custom string into the filehandle that the serializer is using

=cut

sub print_metadata {
	my $self = shift;
	my $text = shift;
	my $fh = $self->{'filehandle'};
	print $fh "\n".$text."\n";
}

=head2 print_main_header

	Arg [1]    : Arrayref of slices going into the file.
	Description: Printing the header text or metadata required for this file format,
	             Re-implement in the serializer.
	Returntype : None
=cut

sub print_main_header {
	my $self = shift;
#	my $arrayref_of_slices = shift;
#	my $fh = $self->{'filehandle'};
	warning("No writer for headers in this format. Nothing done" );
}

=head2 print_sequence 
	Arg [1]    : Bio::EnsEMBL::Slice
	Description: By default, prints a block of FASTA format sequence from the given slice
=cut

sub print_sequence {
	my $self = shift;
	my $slice = shift;
	my $fh = $self->{'filehandle'};
	Bio::EnsEMBL::Utils::SeqDumper->dump_fasta( $slice, $fh);
	$self->{'achieved_something'} = 1;
}

=head2 printed_something
    Description: Check if serializer has printed any useful data. Not accurate with FASTA
                 due to non-reporting dumper.
    Returntype : Boolean
=cut
#TODO: Find a way for SeqDumper to indicate whether it printed anything or just the headers.
sub printed_something {
	my $self = shift;
	if ($self->{'achieved_something'}) { return 1;}
	else {return 0;}
}

1;
