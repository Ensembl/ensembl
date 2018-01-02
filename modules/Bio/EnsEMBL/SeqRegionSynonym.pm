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

Bio::EnsEMBL::SeqRegionSynonym -
Object representing an alternatice name.

=head1 SYNOPSIS

=head1 DESCRIPTION

This object holds information about alternative name to
Ensembl seq regions.

=head1 METHODS

=cut

package Bio::EnsEMBL::SeqRegionSynonym;

use strict;
use warnings;
no warnings qw(uninitialized);

use Bio::EnsEMBL::Storable;
use Bio::Annotation::DBLink;

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(deprecate);

our @ISA = qw(Bio::EnsEMBL::Storable);

=head2 new

  Args [...] : list of named parameters 
  Example    : my $srs = new Bio::EnsEMBL::SeqRegionSynonym(
                    -adaptor => $adaptor,
                    -synonym => $alt_name,
                    -external_db_id => 1234
                    -seq_region_id  => 12);
  Description: Creates a new SeqRegionSynonym object
  Returntype : Bio::EnsEMBL::SeqRegionSynonym
  Exceptions : none
  Caller     : Bio::EnsEMBL::SeqRegionSynonymAdaptor
  Status     : At Risk
=cut

sub new {
  my ($class, @args) = @_;
  
  my $self = bless {},$class;

  my ( $adaptor, $synonym, $ex_db, $seq_region_id, $dbid, $dbname, $db_display_name) =
    rearrange ( ['ADAPTOR','SYNONYM','EXTERNAL_DB_ID','SEQ_REGION_ID','DBID', 'DBNAME', 'DB_DISPLAY_NAME'], @args );

  $self->adaptor($adaptor);

  if( defined $ex_db ) { $self->external_db_id( $ex_db ) } ;
  if( defined $seq_region_id ) { $self->seq_region_id( $seq_region_id ) } ;
  if (defined $dbid) { $self->{'dbID'} = $dbid} ;
  if (defined $dbname) { $self->{'dbname'} = $dbname };
  if (defined $db_display_name) { $self->{'db_display_name'} = $db_display_name };

  if( defined $synonym ) { 
    $self->name( $synonym ) ;
  } else {
    warn "No alternative name given\n";
    return undef;
  }

  return $self;
}

sub name{
 my $self = shift;
  $self->{'name'} = shift if(@_);
  return $self->{'name'};
}

sub external_db_id{
 my $self = shift;
  $self->{'ex_db'} = shift if(@_);
  return $self->{'ex_db'};
}

sub seq_region_id{
  my $self = shift;
  $self->{'seq_region_id'} = shift if(@_);
  return $self->{'seq_region_id'};
}

sub dbname {
  my $self = shift;
  $self->{'dbname'} = shift if(@_);
  return $self->{'dbname'};
}

sub db_display_name {
  my $self = shift;
  $self->{'db_display_name'} = shift if(@_);
  return $self->{'db_display_name'};
}

sub summary_as_hash {
  my $self = shift;
  my %summary;
  $summary{name} = $self->name;
  $summary{dbname} = $self->dbname;

  return \%summary;
}

1;
