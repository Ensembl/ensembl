=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::IdMapping::MappingList - object holding a list of Entries

=head1 SYNOPSIS

  # create a new MappingList
  my $mappings = Bio::EnsEMBL::IdMapping::MappingList->new(
    -DUMP_PATH  => $dump_path,
    -CACHE_FILE => 'gene_mappings.ser',
  );

  # add entries
  my $mappings->add_Entry($entry1);
  my $mappings->add_all( $entry2, $entry3 );

  # serialise to file
  $mappings->write_to_file;

  # later, read these mappings from file
  my $mappings1 = Bio::EnsEMBL::IdMapping::MappingList->new(
    -DUMP_PATH  => $dump_path,
    -CACHE_FILE => 'gene_mappings.ser',
  );
  $mappings1->read_from_file;

=head1 DESCRIPTION

This object represents a list of Bio::EnsEMBL::IdMapping::Entry
objects. It's essentially an OO wrapper for an array with some type
checking and convenience methods.

=head1 METHODS

  new
  add_Entry
  get_all_Entries
  add_all
  get_entry_count
  log
  to_string

=cut

package Bio::EnsEMBL::IdMapping::MappingList;

use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::IdMapping::Serialisable;
our @ISA = qw(Bio::EnsEMBL::IdMapping::Serialisable);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::ScriptUtils qw(path_append);


=head2 new

  Arg[1-N]    : see superclass
  Example     : my $gene_mappings = Bio::EnsEMBL::IdMapping::MappingList->new(
                  -DUMP_PATH   => $dump_path,
                  -CACHE_FILE  => 'gene_mappings.ser',
                );
  Description : Constructor.
  Return type : Bio::EnsEMBL::IdMapping::MappingList
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);

  # initialise internal datastructure unless we loaded a serialised object
  unless ($self->loaded) {
    $self->{'cache'}->{'entries'} = [];
  }

  return $self;
}


=head2 add_Entry

  Arg[1]      : Bio::EnsEMBL::IdMapping::Entry - Entry to add
  Example     : $mappings->add_Entry($entry);
  Description : Adds an Entry to the MappingList.
  Return type : none
  Exceptions  : thrown on wrong or missing argument
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub add_Entry {
  my $self = shift;
  my $entry = shift;

  unless ($entry and $entry->isa('Bio::EnsEMBL::IdMapping::Entry')) {
    throw("Need a Bio::EnsEMBL::IdMapping::Entry");
  }

  push @{ $self->{'cache'}->{'entries'} }, $entry;
}


=head2 get_all_Entries

  Example     : foreach my $entry (@{ $mappings->get_all_Entries }) {
                  # do something with the entry
                }
  Description : Gets all Entries in the MappingList.
  Return type : Arrayref of Bio::EnsEMBL::IdMapping::Entry
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub get_all_Entries {
  my $self = shift;
  return $self->{'cache'}->{'entries'};
}


=head2 add_all

  Arg[1]      : List of Bio::EnsEMBL::IdMapping::Entry objects
  Example     : my @entries = ($entry1, $entry2);
                $mappings->add_all(@entries);
  Description : Adds a list of Entries to the MappingList.
  Return type : none
  Exceptions  : thrown on wrong argument
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub add_all {
  my $self = shift;
  my @mappings = @_;

  foreach my $mapping (@mappings) {
    
    unless ($mapping->isa('Bio::EnsEMBL::IdMapping::MappingList')) {
      throw("Need a Bio::EnsEMBL::IdMapping::MappingList");
    }

    push @{ $self->{'cache'}->{'entries'} }, @{ $mapping->get_all_Entries };
  }
}


=head2 get_entry_count

  Example     : my $num_entries = $mappings->get_entry_count;
  Description : Returns the number of Entries in the MappingList.
  Return type : Int
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub get_entry_count {
  my $self = shift;
  return scalar(@{ $self->{'cache'}->{'entries'} });
}


=head2 log

  Arg[1]      : String $type - object type (e.g. 'gene')
  Arg[2]      : String $dump_path - path for writing output
  Example     : $mappings->log('gene', $conf->param('basedir'));
  Description : Logs all Entries in the MappingList to a file. Used for
                debugging.
  Return type : none
  Exceptions  : thrown on I/0 error
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub log {
  my $self = shift;
  my $type = shift;
  my $dump_path = shift;
  
  my $debug_path = path_append($dump_path, 'debug');
  my $logfile = "$debug_path/${type}_final_scores.txt";
  
  open(my $fh, '>', $logfile) or
    throw("Unable to open $logfile for writing: $!");

  foreach my $entry (@{ $self->get_all_Entries }) {
    print $fh ($entry->to_string."\n");
  }

  close($fh);
}


=head2 to_string

  Example     : print LOG $mappings->to_string, "\n";
  Description : Returns a string representation of the MappingList. This is
                simply a multi-line string, where each line is a stringified
                Entry.
                Useful for debugging and logging.
  Return type : String
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub to_string {
  my $self = shift;
  
  my $string = '';
  
  foreach my $entry (@{ $self->get_all_Entries }) {
    $string .= $entry->to_string."\n";
  }

  return $string;
}


1;

