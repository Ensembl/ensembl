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

Bio::EnsEMBL::Storable

=head1 SYNOPSIS

  my $dbID    = $storable_object->dbID();
  my $adaptor = $storable_object->adaptor();
  if ( $storable_object->is_stored($db_adaptor) ) { ... }

=head1 DESCRIPTION

This is a storable base class.  All objects which are storable
in the database should inherit from this class.  It provides two
getter/setters: dbID() adaptor().  And a is_stored() method that can be
used to determine if an object is already stored in a database.

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Storable;


use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Scalar::Util qw(weaken isweak);

=head2 new

  Arg [-ADAPTOR] : Bio::EnsEMBL::DBSQL::BaseAdaptor
  Arg [-dbID]    : database internal id
  Caller         : internal calls
  Description    : create a new Storable object 
  Returntype     : Bio::EnsEMBL::Storable
  Exceptions     : Adaptor not a Bio::EnsEMBL::DBSQL::BaseAdaptor
  Status         : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($adaptor, $dbID) = rearrange(['ADAPTOR', 'dbID'],@_);

  if($adaptor) {
    if(!ref($adaptor) || !$adaptor->isa('Bio::EnsEMBL::DBSQL::BaseAdaptor')) {
      throw('-ADAPTOR argument must be a Bio::EnsEMBL::DBSQL::BaseAdaptor');
    }
  }

  my $self = bless({'dbID' => $dbID}, $class);
  $self->adaptor($adaptor);
  return $self;
}

=head2 new_fast

  Arg [1]    : hashref to be blessed
  Description: Construct a new Bio::EnsEMBL::Storable object using the hashref.
               This is a very quick constructor that requires internal knowledge 
               of the class. This is used in speed critical sections of the code
               where many objects need to be created quickly.
  Exceptions : none
  Returntype : Instance of Bio::EnsEMBL::Storable subclass
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

=head2 dbID

  Arg [1]    : int $dbID
  Description: getter/setter for the database internal id
  Returntype : int
  Exceptions : none
  Caller     : general, set from adaptor on store
  Status     : Stable

=cut

sub dbID {
  my $self = shift;
  $self->{'dbID'} = shift if(@_);
  return $self->{'dbID'};
}



=head2 adaptor

  Arg [1]    : Bio::EnsEMBL::DBSQL::BaseAdaptor $adaptor
  Description: get/set for this objects Adaptor
  Returntype : Bio::EnsEMBL::DBSQL::BaseAdaptor
  Exceptions : none
  Caller     : general, set from adaptor on store
  Status     : Stable

=cut

sub adaptor {
  my ($self, $adaptor) = @_;
  if(scalar(@_) > 1) {
    if(defined $adaptor) {
      assert_ref($adaptor, 'Bio::EnsEMBL::DBSQL::BaseAdaptor', 'adaptor');
      $self->{adaptor} = $adaptor;
      weaken($self->{adaptor});
    }
    else {
      $self->{adaptor} = undef;
    }
  }
  return $self->{adaptor}
}



=head2 is_stored

  Arg [1]    : Bio::EnsEMBL::DBSQL::DBConnection 
             : or Bio::EnsEMBL::DBSQL::DBAdaptor
  Example    : do_something if($object->is_stored($db));
  Description: Returns true if this object is stored in the provided database.
               This works under the assumption that if the adaptor and dbID are
               set and the database of the adaptor shares the port, dbname and
               hostname with the provided database, this object is stored in
               that database.
  Returntype : 1 or 0
  Exceptions : throw if dbID is set but adaptor is not
               throw if adaptor is set but dbID is not
               throw if incorrect argument is passed
  Caller     : store methods
  Status     : Stable

=cut

my $message_only_once =1;

sub is_stored {
  my $self = shift;
  my $db = shift;

  if($db and $db->isa('Bio::EnsEMBL::DBSQL::DBAdaptor')) {
    $db = $db->dbc();
  }
  if(!$db || !ref($db) || !$db->isa('Bio::EnsEMBL::DBSQL::DBConnection')) {
    throw('db argument must be a Bio::EnsEMBL::DBSQL::DBConnection');
  }

  my $adaptor = $self->{'adaptor'};
  my $dbID = $self->{'dbID'};

  if($dbID && !$adaptor) {
    if($message_only_once){
      warning("Storable object has a dbID but not an adaptor.\n" .
	      'Storable objects must have neither OR both.');
      $message_only_once = 0;
    }
    return 0;
  }

  if($adaptor && !$dbID) {
  	if($message_only_once){
      warning("Storable object has an adaptor but not a dbID.\n".
            "Storable objects must have neither OR both.");
      $message_only_once = 0;
  	}
    return 0;
  }

  return 0 if (!$adaptor && !$dbID);

  my $cur_db = $adaptor->dbc();

  #
  # Databases are the same if they share the same port, host and username
  #
  if ( ($db->port()||'') eq ($cur_db->port()||'')
    && ($db->host()||'') eq ($cur_db->host()||'')
    &&  $db->dbname()    eq  $cur_db->dbname() )
  {
    return 1;
  }

  return 0;
}

sub get_all_DAS_Features{
  my ($self, $slice) = @_;

  $self->{_das_features} ||= {}; # Cache
  $self->{_das_styles} ||= {}; # Cache
  $self->{_das_segments} ||= {}; # Cache
  my %das_features;
  my %das_styles;
  my %das_segments;

  foreach my $dasfact( @{$self->get_all_DASFactories} ){
    my $dsn = $dasfact->adaptor->dsn;
    my $name = $dasfact->adaptor->name;
    my $url = $dasfact->adaptor->url;

    # Construct a cache key : SOURCE_URL/TYPE
    # Need the type to handle sources that serve multiple types of features

    my ($type) = ref($dasfact->adaptor->mapping) eq 'ARRAY' ? @{$dasfact->adaptor->mapping} : $dasfact->adaptor->mapping;
    $type ||=$dasfact->adaptor->type;
    my $key = join('/', $name, $type);

    if( $self->{_das_features}->{$key} ){ # Use cached
        $das_features{$name} = $self->{_das_features}->{$key};
        $das_styles{$name} = $self->{_das_styles}->{$key};
        $das_segments{$name} = $self->{_das_segments}->{$key};
    } else { # Get fresh data
  
        my ($featref, $styleref, $segref) = ($type =~ /^ensembl_location/) ?  ($dasfact->fetch_all_Features( $slice, $type )) : $dasfact->fetch_all_by_ID( $self );
   
        $self->{_das_features}->{$key} = $featref;
        $self->{_das_styles}->{$key} = $styleref;
        $self->{_das_segments}->{$key} = $segref;
        $das_features{$name} = $featref;
        $das_styles{$name} = $styleref;
        $das_segments{$name} = $segref;
    }
  }

  return (\%das_features, \%das_styles, \%das_segments);
}

1;
