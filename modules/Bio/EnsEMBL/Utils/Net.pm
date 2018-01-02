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

package Bio::EnsEMBL::Utils::Net;

=pod


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=pod

=head1 NAME

Bio::EnsEMBL::Utils::Net

=head1 SYNOPSIS

  use Bio::EnsEMBL::Utils::Net qw/do_GET do_FTP/;

  #Doing a HTTP get
  my $google_contents = do_GET('http://www.google.co.uk/');
  
  #Doing a FTP request; delegates onto LWP
  my $ftp_contents = do_GET('ftp://ftp.ensembl.org/pub/current_README');

=head1 DESCRIPTION

A collection of subroutines aimed to helping network based operations. The code
will use HTTP::Tiny for all HTTP operations if available. Otherwise it
will delegate to LWP. LWP is currently the only supported target for FTP.

=head1 METHODS

See subroutines.

=head1 VERSION

$Revision$

=cut

use strict;
use warnings;

use base qw/Exporter/;
use Time::HiRes;

our @EXPORT_OK;

@EXPORT_OK = qw(
  do_GET
  do_FTP
  do_FTP_to_file
);

use Bio::EnsEMBL::Utils::Exception qw(throw);

our $HTTP_TINY = 0;
our $LWP = 0;
eval {
  require HTTP::Tiny;
  $HTTP_TINY = 1;
};
eval {
  require LWP::UserAgent;
  $LWP = 1;
};

=head2 do_GET

  Arg [1]     : string $url The URL to fetch including all parameters
  Arg [2]     : int; $total_attempts The number of times to try the URL 
                before throwing an exception
  Arg [3]     : number; $sleep Amount of time to sleep between attempts. 
                Delegates onto Time::HiRes so floating point numbers are 
                supported
  Description : Performs a HTTP GET method call to return the specified remote
                resource.
  Returntype  : Scalar of the contents of the remote URL. Do not use to
                retrieve very large amounts of data.
  Example     : my $contents = do_GET('http://www.google.co.uk/');
  Exceptions  : If we could not retrieve the resource after the specified 
                number of attempts.
  Status      : Stable

=cut

sub do_GET {
  my ($url, $total_attempts, $sleep) = @_;
  return _retry_sleep(sub {
    if($HTTP_TINY) {
      return _get_http_tiny($url);
    }
    elsif($LWP) {
      return _get_lwp($url);
    }
    else {
      throw "Cannot continue. You do not have HTTP::Tiny or LWP available."
    }
  }, $total_attempts, $sleep);
}

=head2 do_FTP

  Arg [1]     : string $uri
  Arg [2]     : int; $total_attempts The number of times to try the URI 
                before throwing an exception
  Arg [3]     : number; $sleep Amount of time to sleep between attempts. 
                Delegates onto Time::HiRes so floating point numbers are 
                supported
  Description : Performs a FTP fetch using a non-authenticated connection (
                however some servers will allow you to encode this in the URI).
  Returntype  : Scalar of the contents of the remote URL. Do not use to
                retrieve very large amounts of data.
  Example     : my $contents = do_GET('http://www.google.co.uk/');
  Exceptions  : If we could not retrieve the resource after the specified 
                number of attempts.
  Status      : Stable

=cut

sub do_FTP {
  my ($url, $total_attempts, $sleep) = @_;
  return _retry_sleep(sub {
    return _get_lwp($url);
  }, $total_attempts, $sleep);
}

=head2 do_FTP_to_file

  Arg [1]     : string $uri
  Arg [2]     : int; $total_attempts The number of times to try the URI 
                before throwing an exception
  Arg [3]     : number; $sleep Amount of time to sleep between attempts. 
                Delegates onto Time::HiRes so floating point numbers are 
                supported
  Description : Performs a FTP fetch using a non-authenticated connection (
                however some servers will allow you to encode this in the URI).
  Returntype  : Boolean true if download was successful.
  Example     : my $contents = do_GET('http://www.google.co.uk/');
  Exceptions  : If we could not retrieve the resource after the specified 
                number of attempts.
  Status      : Stable

=cut

sub do_FTP_to_file {
  my ($url, $total_attempts, $sleep, $filename) = @_;
  return _retry_sleep(sub {
    return _get_lwp_to_file($url, $filename);
  }, $total_attempts, $sleep);
}

sub _retry_sleep {
  my ($callback, $total_attempts, $sleep) = @_;
  $total_attempts ||= 1;
  $sleep ||= 0;
  my $response;
  my $retries = 0;
  my $fail = 1;
  while($retries <= $total_attempts) {
    $response = $callback->();
    if(defined $response) {
      $fail = 0;
      last;
    }
    $retries++;
    Time::HiRes::sleep($sleep);
  }
  if($fail) {
    throw "Could not request remote resource after $total_attempts attempts";
  }
  return $response;
}

sub _get_http_tiny {
  my ($url) = @_;
  my $response = HTTP::Tiny->new->get($url);
  return unless $response->{success};
  return $response->{content} if length $response->{content};
  return;
}

sub _get_lwp {
  my ($url) = @_;
  throw "Cannot perform action as LWP::UserAgent is not available" unless $LWP;
  my $ua = LWP::UserAgent->new();
  $ua->env_proxy;
  my $response = $ua->get($url);
  return $response->decoded_content if $response->is_success;
  return;
}

sub _get_lwp_to_file {
  my ($url, $filename) = @_;
  throw "Cannot perform action as LWP::UserAgent is not available" unless $LWP;
  throw "Filename required for download to proceed." unless $filename;
  my $ua = LWP::UserAgent->new();
  $ua->env_proxy;
  my $response = $ua->get($url, ":content_file" => $filename);
  print 'UA Response: '.$response->is_success."\n";
  if ($response->is_success) {return 1}
  throw $response->status_line;
}

1;
