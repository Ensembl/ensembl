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

Bio::EnsEMBL::Utils::URIParser

=head1 SYNOPSIS

  #Using OO
  use Bio::EnsEMBL::Utils::URIParser;
  my $up = Bio::EnsEMBL::Utils::URIParser->new();
  
  my $db_uri = $up->parse('mysql://user@host:3157/db');
  my $http_uri = $up->parse('http://www.google.co.uk:80/search?q=t');
  
  #Going functional
  use Bio::EnsEMBL::Utils::URIParser qw/parse_uri/;
  my $db_uri = parse_uri('mysql://user@host:3157/db');
  my $http_uri = parse_uri('http://www.google.co.uk:80/search?q=t');

=head1 DESCRIPTION

This object is a generic URI parser which is primarily used in the
parsing of database URIs into a more managable data structure. It's aim
is this conversion and nothing else. The code produces URI objects
which contain more methods for introspection of the data points.

=head1 METHODS

=cut

package Bio::EnsEMBL::Utils::URIParser;

use strict;
use warnings;

use base qw/Exporter/;
our @EXPORT_OK;
@EXPORT_OK = qw/parse_uri/;

sub new {
  my ($class) = @_;
  return bless({}, ref($class) || $class);
}

=head2 parse_uri

  Arg[1]      : Scalar; URI to parse
  Example     : my $p = parse_uri('mysql://user:pass@host:415/db');
  Description : Functional form of the C<parse()> call.
  Returntype  : Bio::EnsEMBL::Utils::URIParser::URI
  Caller      : General
  Status      : Stable
  
=cut

sub parse_uri {
  my ($uri) = @_;
  return __PACKAGE__->new()->parse($uri);
}

=head2 parse

  Arg[1]      : Scalar; URI to parse
  Example     : my $uri = $p->parse('mysql://user:pass@host:415/db');
  Description : A URL parser which attempts to convert many different types
                of URL into a common data structure.
  Returntype  : Bio::EnsEMBL::Utils::URIParser::URI
  Caller      : General
  Status      : Beta
  
=cut

sub parse {
  my ($self, $url) = @_;
  
  my $SCHEME = qr{ ([^:]+) :// }xms;
  my $USER = qr{ ([^/:\@]+)? :? ([^/\@]+)? \@ }xms;
  my $HOST = qr{ ([^/:]+)? :? ([^/]+)? }xms;
  my $DB = qr{ / ([^/?]+)? /? ([^/?]+)? }xms;
  my $PARAMS = qr{ \? (.+)}xms;
  
  my $p;

  if($url =~ qr{ $SCHEME ([^?]+) (?:$PARAMS)? }xms) {
    my $scheme = $1;
    $p = Bio::EnsEMBL::Utils::URIParser::URI->new($scheme);
    my ($locator, $params) = ($2, $3);
    
    if($scheme eq 'file') {
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
      
      if($p->is_db_scheme()) {
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

1;