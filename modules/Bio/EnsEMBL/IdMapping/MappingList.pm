package Bio::EnsEMBL::IdMapping::MappingList;

=head1 NAME


=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS


=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http://www.ensembl.org/info/about/code_licence.html for details.

=head1 AUTHOR

Patrick Meidl <meidl@ebi.ac.uk>, Ensembl core API team

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<ensembl-dev@ebi.ac.uk>

=cut


use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::IdMapping::Serialisable;
our @ISA = qw(Bio::EnsEMBL::IdMapping::Serialisable);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::ScriptUtils qw(path_append);


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


sub add_Entry {
  my $self = shift;
  my $entry = shift;

  unless ($entry and $entry->isa('Bio::EnsEMBL::IdMapping::Entry')) {
    throw("Need a Bio::EnsEMBL::IdMapping::Entry");
  }

  push @{ $self->{'cache'}->{'entries'} }, $entry;
}


sub get_all_Entries {
  my $self = shift;
  return $self->{'cache'}->{'entries'};
}


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


sub get_entry_count {
  my $self = shift;
  return scalar(@{ $self->{'cache'}->{'entries'} });
}


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


sub to_string {
  my $self = shift;
  
  my $string = '';
  
  foreach my $entry (@{ $self->get_all_Entries }) {
    $string .= $entry->to_string."\n";
  }

  return $string;
}


1;

