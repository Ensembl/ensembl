#
# Ensembl module for Bio::EnsEMBL::DBSQL::MiscSetAdaptor
#
# Copyright (c) 2003 EnsEMBL
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::MiscSetAdaptor - Provides database interaction for
Bio::EnsEMBL::MiscSet objects.


=head1 SYNOPSIS

  #$db is a Bio::EnsEMBL::DBSQL::DBAdaptor object:
  $msa = $db->get_MiscSetAdaptor();

  $misc_set = $msa->fetch_by_dbID(1234);

  $misc_set = $msa->fetch_by_code('clone');


=head1 DESCRIPTION

This class provides database interactivity for MiscSet objects.  MiscSets
are used to classify MiscFeatures into groups.

=head1 CONTACT

This modules is part of the Ensembl project http://www.ensembl.org

Questions can be posted to the ensembl-dev mailing list:
ensembl-dev@ebi.ac.uk

=cut

use strict;
use warnings;

package Bio::EnsEMBL::DBSQL::MiscSetAdaptor;

use Bio::EnsEMBL::MiscSet;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);


=head2 new

  Arg [...]  : Superclass args.  See Bio::EnsEMBL::DBSQL::BaseAdaptor
  Description: Instantiates a Bio::EnsEMBL::DBSQL::MiscSetAdaptor and
               caches the contents of the MiscSet table.
  Returntype : Bio::EnsEMBL::MiscSet
  Exceptions : none
  Caller     : MiscFeatureAdaptor

=cut


sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  #cache the entire contents of the misc set table
  #the table is small and it removes the need to repeatedly query the
  #table or join to the table

  $self->{'_id_cache'} = {};
  $self->{'_code_cache'} = {};

  my $sth = $self->prepare
    ('SELECT misc_set_id, code, name, description, max_length FROM misc_set');

  $sth->execute();

  my ($dbID, $code, $name, $desc, $max_len);
  $sth->bind_columns(\$dbID, \$code, \$name, \$desc, \$max_len);

  while($sth->fetch()) {
    my $ms =
      Bio::EnsEMBL::MiscSet->new($dbID, $self, $code, $name, $desc, $max_len);

    $self->{'_id_cache'}->{$dbID} = $ms;
    $self->{'_code_cache'}->{lc($code)} = $ms;
  }

  $sth->finish();

  return $self;
}




=head2 fetch_all

  Arg [1]    : none
  Example    : foreach my $ms (@{$msa->fetch_all()}) {
                 print $ms->code(), ' ', $ms->name(), "\n";
               }
  Description: Retrieves every MapSet defined in the DB
  Returntype : listref of Bio::EnsEMBL::MapSets
  Exceptions : none
  Caller     : general

=cut

sub fetch_all {
  my $self = shift;

  my @all = values(%{$self->{'_id_cache'}});

  return \@all;
}



=head2 fetch_by_dbID

  Arg [1]    : int $dbID
               The internal identifier of the misc set to retrieve
  Example    : my $ms = $msa->fetch_by_dbID($dbID);
  Description: Retrieves a misc set via its internal identifier
  Returntype : Bio::EnsEMBL::MiscSet
  Exceptions : none
  Caller     : general

=cut

sub fetch_by_dbID {
  my $self = shift;
  my $dbID = shift;

  return $self->{'_id_cache'}->{$dbID};
}



=head2 fetch_by_code

  Arg [1]    : string $code
               The unique code of the MiscSet to retrieve
  Example    : my $ms = $msa->fetch_by_code('clone');
  Description: Retrieves a MiscSet via its code
  Returntype : Bio::EnsEMBL::MiscSet
  Exceptions : none
  Caller     : general

=cut

sub fetch_by_code {
  my $self = shift;
  my $code = shift;

  return $self->{'_code_cache'}->{lc($code)};
}



#
# Called during db destruction to clean up internal cache structures
# that result in circular references
#
sub deleteObj {
  my $self = shift;

  #break circular db <-> adaptor references
  $self->SUPER::deleteObj();

  #break circular object <-> adaptor references
  delete $self->{'_id_cache'};
  delete $self->{'_code_cache'};
}


1;
