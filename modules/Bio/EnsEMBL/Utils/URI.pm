=head1 LICENSE

  Copyright (c) 1999-2011 The European Bioinformatics Institute and
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

Bio::EnsEMBL::Utils::URI

=head1 SYNOPSIS

  use Bio::EnsEMBL::Utils::URI qw/parse_uri/;
  my $up = Bio::EnsEMBL::Utils::URIParser->new();

  my $db_uri = parse_uri('mysql://user@host:3157/db');
  my $http_uri = parse_uri('http://www.google.co.uk:80/search?q=t');

=head1 DESCRIPTION

This object is a generic URI parser which is primarily used in the
parsing of database URIs into a more managable data structure. We also provide
the resulting URI object

=head1 METHODS

=cut

package Bio::EnsEMBL::Utils::URI;

use strict;
use warnings;

use Scalar::Util qw/looks_like_number/;
use Bio::EnsEMBL::Utils::Exception qw(throw);
use File::Spec;

use base qw/Exporter/;
our @EXPORT_OK;
@EXPORT_OK = qw/parse_uri/;

####URI Parsing

=head2 parse_uri

  Arg[1]      : Scalar; URI to parse
  Example     : my $uri = parse_uri('mysql://user:pass@host:415/db');
  Description : A URL parser which attempts to convert many different types
                of URL into a common data structure.
  Returntype  : Bio::EnsEMBL::Utils::URI
  Caller      : General
  Status      : Beta

=cut

sub parse_uri {
  my ($url) = @_;

  my $SCHEME = qr{ ([^:]*) :// }xms;
  my $USER = qr{ ([^/:\@]+)? :? ([^/\@]+)? \@ }xms;
  my $HOST = qr{ ([^/:]+)? :? ([^/]+)? }xms;
  my $DB = qr{ / ([^/?]+)? /? ([^/?]+)? }xms;
  my $PARAMS = qr{ \? (.+)}xms;

  my $p;

  if($url =~ qr{ $SCHEME ([^?]+) (?:$PARAMS)? }xms) {
    my $scheme = $1;
    $p = Bio::EnsEMBL::Utils::URI->new($scheme);
    my ($locator, $params) = ($2, $3);

    if($scheme eq 'file') {
      $p->path($locator);
    }
    elsif($scheme eq 'sqlite') {
      $p->path($locator);
    }
    else {
      if($locator =~ s/^$USER//) {
        $p->user($1);
        $p->pass($2);
      }
      if($locator =~ s/^$HOST//) {
        $p->host($1);
        $p->port($2);
      }

      if($p->is_db_scheme() || $scheme eq q{}) {
        if($locator =~ $DB) {
          $p->db_params()->{dbname} = $1;
          $p->db_params()->{table} = $2;
        }
      }
      else {
        $p->path($locator);
      }
    }

    if(defined $params) {
      my @kv_pairs = split(/;|&/, $params);
      foreach my $kv_string (@kv_pairs) {
        my ($key, $value) = split(/=/, $kv_string);
        $p->add_param($key, $value);
      }
    }
  }

  return $p;
}

####URI Object

=pod

=head2 new()

  Arg[1]      : String; scheme the URI will confrom to
  Description : New object call
  Returntype  : Bio::EnsEMBL::Utils::URIParser::URI
  Exceptions  : Thrown if scheme is undefined.
  Status      : Stable

=cut

sub new {
  my ($class, $scheme) = @_;
  $class = ref($class) || $class;
  throw "Scheme cannot be undefined. Empty string is allowed" if ! defined $scheme;

  my $self = bless ({
    params => {},
    param_keys => [],
    db_params => {},
    scheme => $scheme,
  }, $class);

  return $self;
}

=head2 db_schemes()

  Description: Returns a hash of scheme names known to be databases
  Returntype : HashRef
  Exceptions : None
  Status     : Stable

=cut

sub db_schemes {
  my ($self) = @_;
  return {map { $_ => 1 } qw/mysql ODBC sqlite Oracle Sybase/};
}


=head2 is_db_scheme()

  Description: Returns true if the code believes the scheme to be a Database
  Returntype : Boolean
  Exceptions : None
  Status     : Stable

=cut

sub is_db_scheme {
  my ($self) = @_;
  return ( exists $self->db_schemes()->{$self->scheme()} ) ? 1 : 0;
}

=head2 scheme()

  Description : Getter for the scheme attribute
  Returntype  : String
  Exceptions  : None
  Status      : Stable

=cut

sub scheme {
  my ($self) = @_;
  return $self->{scheme};
}

=head2 path()

  Arg[1]      : Setter argument
  Description : Getter/setter for the path attribute
  Returntype  : String
  Exceptions  : None
  Status      : Stable

=cut

sub path {
  my ($self, $path) = @_;
  $self->{path} = $path if defined $path;
  return $self->{path};
}

=head2 user()

  Arg[1]      : Setter argument
  Description : Getter/setter for the user attribute
  Returntype  : String
  Exceptions  : None
  Status      : Stable

=cut

sub user {
  my ($self, $user) = @_;
  $self->{user} = $user if defined $user;
  return $self->{user};
}

=head2 pass()

  Arg[1]      : Setter argument
  Description : Getter/setter for the password attribute
  Returntype  : String
  Exceptions  : None
  Status      : Stable

=cut

sub pass {
  my ($self, $pass) = @_;
  $self->{pass} = $pass if defined $pass;
  return $self->{pass};
}

=head2 host()

  Arg[1]      : Setter argument
  Description : Getter/setter for the host attribute
  Returntype  : String
  Exceptions  : None
  Status      : Stable

=cut

sub host {
  my ($self, $host) = @_;
  $self->{host} = $host if defined $host;
  return $self->{host};
}

=head2 port()

  Arg[1]      : Setter argument
  Description : Getter/setter for the port attribute
  Returntype  : Integer
  Exceptions  : If port is not a number, less than 1 or not a whole integer
  Status      : Stable

=cut

sub port {
  my ($self, $port) = @_;
  if(defined $port) {
    if(! looks_like_number($port) || $port < 1 || int($port) != $port) {
      throw "Port $port is not a number, less than 1 or not a whole integer";
    }
    $self->{port} = $port if defined $port;
  }
  return $self->{port};
}

=head2 param_keys()

  Description : Getter for the paramater map keys in the order they were first
                seen. Keys should only appear once in this array
  Returntype  : ArrayRef
  Exceptions  : None
  Status      : Stable

=cut

sub param_keys {
  my ($self) = @_;
  return [@{$self->{param_keys}}];
}

=head2 param_exists_ci()

  Arg[1]      : String; Key
  Description : Performs a case-insensitive search for the given key
  Returntype  : Boolean; returns true if your given key was seen
  Exceptions  : None
  Status      : Stable

=cut

sub param_exists_ci {
  my ($self, $key) = @_;
  my %keys = map { uc($_) => 1 } @{$self->param_keys()};
  return ($keys{uc($key)}) ? 1 : 0;
}

=head2 add_param()

  Arg[1]      : String; key
  Arg[1]      : Scalar; value
  Description : Add a key/value to the params map. Multiple inserts of the same
                key is allowed
  Returntype  : None
  Exceptions  : None
  Status      : Stable

=cut

sub add_param {
  my ($self, $key, $value) = @_;
  if(!exists $self->{params}->{$key}) {
    $self->{params}->{$key} = [];
    push(@{$self->{param_keys}}, $key);
  }
  push(@{$self->{params}->{$key}}, $value);
  return;
}

=head2 get_params()

  Arg[1]      : String; key
  Description : Returns the values which were found to be linked to the given
                key. Arrays are returned because one key can have many
                values in a URI
  Returntype  : ArrayRef[Scalar]
  Exceptions  : None
  Status      : Stable

=cut

sub get_params {
  my ($self, $key) = @_;
  return [] if ! exists $self->{params}->{$key};
  return [@{$self->{params}->{$key}}];
}

=head2 db_params()

  Description : Storage of parameters used only for database URIs since
                they require
  Returntype  : HashRef; Database name is keyed under C<dbname> and the
                table is keyed under C<table>
  Exceptions  : None
  Status      : Stable

=cut

sub db_params {
  my ($self) = @_;
  return $self->{db_params};
}

=head2 generate_dbsql_params()

  Arg[1]      : boolean $no_table alows you to avoid pushing -TABLE as an 
                output value
  Description : Generates a Hash of Ensembl compatible parameters to be used
                to construct a DB object. We combine those parameters
                which are deemed to be part of the C<db_params()> method
                under C<-DBNAME> and C<-TABLE>. We also search for a number
                of optional parameters which are lowercased equivalents
                of the construction parameters available from a
                L<Bio::EnsEMBL::DBSQL::DBAdaptor>,
                L<Bio::EnsEMBL::DBSQL::DBConnection> as well as C<verbose>
                being supported.

                We also convert the scheme type into the driver attribute

  Returntype  : Hash (not a reference). Output can be put into a C<DBConnection>
                constructor.
  Exceptions  : None
  Status      : Stable

=cut

sub generate_dbsql_params {
  my ($self, $no_table) = @_;
  my %db_params;

  $db_params{-DRIVER} = $self->scheme();
  $db_params{-HOST}   = $self->host() if $self->host();
  $db_params{-PORT}   = $self->port() if $self->port();
  $db_params{-USER}   = $self->user() if $self->user();
  $db_params{-PASS}   = $self->pass() if $self->pass();
  
  my $dbname;
  my $table;
  if($self->scheme() eq 'sqlite') {
    ($dbname, $table) = $self->_decode_sqlite();
  }
  else {
    $dbname = $self->db_params()->{dbname};
    $table = $self->db_params()->{table};
  }
  
  $db_params{-DBNAME} = $dbname if $dbname;
  $db_params{-TABLE}  = $table if ! $no_table && $table;

  foreach my $boolean_param (qw/disconnect_when_inactive reconnect_when_connection_lost is_multispecies no_cache verbose/) {
    if($self->param_exists_ci($boolean_param)) {
      $db_params{q{-}.uc($boolean_param)} = 1;
    }
  }
  foreach my $value_param (qw/species group species_id wait_timeout/) {
    if($self->param_exists_ci($value_param)) {
      $db_params{q{-}.uc($value_param)} = $self->get_params($value_param)->[0];
    }
  }

  return %db_params;
}

=head2 _decode_sqlite

  Description : Performs path gymnastics to decode into a number of possible
                options. The issue with SQLite is that the normal URI scheme
                looks like sqlite:///my/path.sqlite/table but how do we know
                that the DB name is C</my/path.sqlite> and the table is 
                C<table>?
                
                The code takes a path, looks for the full path & if it cannot
                be found looks for the file a directory back. In the above
                example it would have looked for C</my/path.sqlite/table>,
                found it to be non-existant, looked for C</my/path.sqlite>
                and found it. 
                
                If the path splitting procdure resulted in just 1 file after
                the first existence check e.g. C<sqlite://db.sqlite> it assumes
                that should be the name. If no file can be found we default to
                the full length path.
  Caller      : internal

=cut

sub _decode_sqlite {
  my ($self) = @_;
  my $dbname;
  my $table;
  my $path = $self->path();
  if(-f $path) {
    $dbname = $path;
  }
  else {
    my ($volume, $directories, $file) = File::Spec->splitpath($path);
    my @splitdirs = File::Spec->splitdir($directories);
    if(@splitdirs == 1) {
      $dbname = $path;
    }
    else {
      my $new_file = pop(@splitdirs);
      $new_file ||= q{};
      my $new_path = File::Spec->catpath($volume, File::Spec->catdir(@splitdirs), $new_file);
      if($new_path ne File::Spec->rootdir() && -f $new_path) {
        $dbname = $new_path;
        $table = $file;
      }
      else {
        $dbname = $path;
      }
    }
  }
  
  $self->db_params()->{dbname} = $dbname if $dbname;
  $self->db_params()->{table} = $table if $table;
  
  return ($dbname, $table);
}

=head2 generate_uri()

  Description : Generates a URI string from the paramaters in this object
  Returntype  : String
  Exceptions  : None
  Status      : Stable

=cut

sub generate_uri {
  my ($self) = @_;
  my $scheme = sprintf('%s://', $self->scheme());
  my $user_credentials = q{};
  my $host_credentials = q{};
  my $location = q{};

  if($self->user() || $self->pass()) {
    $user_credentials = sprintf('%s%s@',
      ( $self->user() ? $self->user() : q{} ),
      ( $self->pass() ? q{:}.$self->pass() : q{} )
    );
  }

  if($self->host() || $self->port()){
    $host_credentials = sprintf('%s%s',
      ($self->host() ? $self->host() : q{}),
      ($self->port() ? q{:}.$self->port() : q{})
    );
  }

  if($self->is_db_scheme() || $self->scheme() eq '') {
    if($self->scheme() eq 'sqlite') {
      if(! $self->path()) {
        my $tmp_loc = $self->db_params()->{dbname};
        throw "There is no dbname available" unless $tmp_loc;
        $tmp_loc .= q{/}.$self->db_params()->{table} if $self->db_params()->{table};
        $self->path($tmp_loc);
      }
      $location = $self->path();
    }
    else {
      if($self->db_params()->{dbname} || $self->db_params()->{table}) {
        $location = sprintf('/%s%s',
          ($self->db_params()->{dbname} ? $self->db_params()->{dbname} : q{}),
          ($self->db_params()->{table} ? q{/}.$self->db_params()->{table} : q{})
        );
      }
    }
  }
  else {
    $location = $self->path() if $self->path();
  }

  my $param_string = q{};
  if(@{$self->param_keys()}) {
    $param_string = q{?};
    my @params;
    foreach my $key (@{$self->param_keys}) {
      my $values_array = $self->get_params($key);
      foreach my $value (@{$values_array}) {
        push(@params, (defined $value) ? "$key=$value" : $key);
      }
    }
    $param_string .= join(q{;}, @params);
  }

  return join(q{}, $scheme, $user_credentials, $host_credentials, $location, $param_string);
}

1;