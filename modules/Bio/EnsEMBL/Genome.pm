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


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Genome - A generic Genome class.

=head1 SYNOPSIS

  use Bio::EnsEMBL::Genome;

  my $genome = Bio::EnsEMBL::Genome->new
       (-NAME => 'myName',
        -CODE => 'MyCode',
        -DESCRIPTION => 'This is my statistics description.',
        -VALUE => '10023');

  print $statistic->name(), "\n";
  print $statistic->code(), "\n";
  print $statistic->description(), "\n";
  print $statistic->value(), "\n";

=head1 DESCRIPTION

This is a generic genome class used to represent statistics
associated with a genome and species

=head1 SEE ALSO

Bio::EnsEMBL::DBSQL::GenomeContainer

=cut

package Bio::EnsEMBL::Genome;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Scalar::Util qw(weaken isweak);

=head2 new

  Arg [-NAME]        : string - the name for this statistic
  Arg [-CODE]        : string - the code for this statistic
  Arg [-DESCRIPTION] : string - a description for this statistic
  Arg [-VALUE]       : value  - the value of this statistic
  Example            :   my $genome = Bio::EnsEMBL::Genome->new
                         (-name => 'myName',
                          -CODE => 'MyCode',
                          -DESCRIPTION => 'This is my statistics description.',
                          -VALUE => '10023');  
  Description        : Constructor.  Instantiates a Bio::EnsEMBL::Genome object.
  Returntype         : Bio::EnsEMBL::Genome
  Exceptions         : none
  Caller             : general
  Status             : Stable

=cut


sub new {
  my $caller = shift;

  # allow to be called as class or object method
  my $class = ref($caller) || $caller;

  my ($dbID, $statistic, $value, $species_id, $code, $name, $desc) =
    rearrange([qw(DBID STATISTIC VALUE SPECIES_ID CODE NAME DESCRIPTION)], @_);

  return bless {'dbID'        => $dbID,
                'statistic'   => $statistic,
                'value'       => $value,
                'species_id'  => $species_id,
                'code'        => $code,
                'name'        => $name,
                'description' => $desc}, $class;
}

=head2 new_fast

  Arg [1]    : hashref to be blessed
  Description: Construct a new Bio::EnsEMBL::Genome using the hashref.
  Exceptions : none
  Returntype : Bio::EnsEMBL::Genome
  Caller     : general, subclass constructors
  Status     : Stable

=cut


sub new_fast {
  my $class = shift;
  my $hashref = shift;
  my $self = bless $hashref, $class;
  weaken($self->{adaptor})  if ( ! isweak($self->{adaptor}) );
  return $self;
}

=head2 statistic

  Arg [1]    : string $statistic (optional)
  Example    : $statistic = $statistic->name();
  Description: Getter/Setter for the statistic name
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub statistic {
  my $self = shift;
  $self->{'statistic'} = shift if(@_);
  return $self->{'statistic'};
}

=head2 code

  Arg [1]    : string $code (optional)
  Example    : $code = $genome->code();
  Description: Getter/Setter for code statistic
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub code {
  my $self = shift;
  $self->{'code'} = shift if(@_);
  return $self->{'code'};
}


=head2 name

  Arg [1]    : string $name (optional)
  Example    : $name = $statistic->name();
  Description: Getter/Setter for the name associated to the statistic attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub name {
  my $self = shift;
  $self->{'name'} = shift if(@_);
  return $self->{'name'};
}

=head2 description

  Arg [1]    : string $description (optional)
  Example    : $description = $statistic->description();
  Description: Getter/Setter for description statistic
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub description {
  my $self = shift;
  $self->{'description'} = shift if(@_);
  return $self->{'description'};
}


=head2 value

  Arg [1]    : string $value (optional)
  Example    : $value = $statistic->value();
  Description: Getter/Setter for value statistic
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub value {
  my $self = shift;
  $self->{'value'} = shift if(@_);
  return $self->{'value'};
}

=head2 dbID

  Arg [1]    : string $value (optional)
  Example    : $dbID = $statistic->dbID();
  Description: Getter/Setter for statistic internal id
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub dbID {
  my $self = shift;
  $self->{'dbID'} = shift if(@_);
  return $self->{'dbID'};
}

=head2 species_id

  Arg [1]    : string $value (optional)
  Example    : $species_id = $statistic->dbID();
  Description: Getter/Setter for statistic species id
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub species_id {
  my $self = shift;
  $self->{'species_id'} = shift if(@_);
  return $self->{'species_id'};
}


1;
