use strict;
use warnings;

use Test::More;
use Test::Exception;

require_ok('Bio::EnsEMBL::Utils::Net');
Bio::EnsEMBL::Utils::Net->import('do_GET');

my $fake_httpd = 0;
eval {
  require Test::Fake::HTTPD;
  $fake_httpd = 1;
};

my $http_test_count = 5;
my $retry_count = 0;

note 'Performing HTTP::Tiny tests';
if($Bio::EnsEMBL::Utils::Net::HTTP_TINY) {
  local $Bio::EnsEMBL::Utils::Net::LWP = 0;
  run_http_test();
}

note 'Performing LWP::UserAgent tests';
if($Bio::EnsEMBL::Utils::Net::LWP) {
  local $Bio::EnsEMBL::Utils::Net::HTTP_TINY = 0;
  run_http_test();
}

done_testing();

sub run_http_test {
  SKIP: {
    skip 'Test::Fake::HTTPD is not available. Cannot perform tests', $http_test_count unless $fake_httpd;
    
    my $httpd = get_server();
    my $endpoint = $httpd->endpoint;
    
    is(do_GET($endpoint), '/', 'Basic successful fetch');
    
    $retry_count = 0;
    dies_ok { do_GET($endpoint.'/retry', 1, 0.1) }  '/retry needs more than 1 attempt';
    $retry_count = 0;
    is(do_GET($endpoint.'/retry', 2, 0.1), '/retry', '/retry needs more than 1 fetch');
    
    dies_ok { do_GET($endpoint.'/404') } 'Sending a bad URL makes it all die very quickly';
  };
}

sub get_server {
  my ($self) = @_;
  my $httpd = Test::Fake::HTTPD->new(
    timeout => 5,
  );
  
  $httpd->run(sub {
    my ($req) = @_;
    my $uri = $req->uri;
    return do {
      if( $uri->path eq '/' ) {
        [ 200, [ 'Content-Type', 'text/plain' ], [ '/' ] ];
      }
      elsif( $uri->path eq '/retry' ) {
        if($retry_count <= 1) {
          $retry_count++;
          [400, [ 'Content-Type', 'text/plain' ], ['']];
        }
        else {
          [ 200, [ 'Content-Type', 'text/plain' ], [ '/retry' ] ];
        }
      }
      else {
        [404, [ 'Content-Type', 'text/plain' ], ['Unsupported URL']];
      }
    }
  });
  ok( defined $httpd, 'Got a web server' );
  diag( sprintf "You can connect to your server at %s.\n", $httpd->host_port );
  return $httpd;
}

