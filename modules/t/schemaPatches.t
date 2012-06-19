use strict;
use warnings;

use Test::More;

use Bio::EnsEMBL::ApiVersion qw/software_version/;
use Bio::EnsEMBL::Test::MultiTestDB;
use File::Find;
use File::Spec::Functions qw/updir catfile catdir/;
use File::Temp qw/tempfile/;
use FindBin qw/$Bin/;

#Using the best available HTTP client
my $tiny = 0;
my $lwp = 0;
eval {
  require HTTP::Tiny;
  note 'Using HTTP::Tiny';
  $tiny = 1;
};
eval {
  if(!$tiny) {
    require LWP::UserAgent;
    note 'Using LWP';
    $lwp = 1;
  }
};

sub get_url {
  my ($url) = @_;
  my $contents;
  if($tiny) {
    my $response = HTTP::Tiny->new->get($url);
    return unless $response->{success};
    return $response->{content} if length $response->{content};
    return;
  }
  elsif($lwp) {
    my $ua = LWP::UserAgent->new;
    $ua->agent("EnsemblTest/1.0.0 ");
    $ua->env_proxy;
    my $response = $ua->get($url);
    return $response->decoded_content if $response->is_success;
    return;
  }
  fail("We do not have access to HTTP::Tiny or LWP. Cannot continue");
  return;
}

my $db = Bio::EnsEMBL::Test::MultiTestDB->new();
my $dba = $db->get_DBAdaptor('core');
my $dbc = $dba->dbc();

#Get last DB version and download the last SQL schema
my $current_release = software_version();
my $last_release = $current_release - 1;
my $cvs_url = "http://cvs.sanger.ac.uk/cgi-bin/viewvc.cgi/ensembl/sql/table.sql?root=ensembl&view=co&pathrev=branch-ensembl-${last_release}";
my $table_sql = get_url($cvs_url);

# Get patch location
my $sql_dir = catdir($Bin, updir(), updir(), 'sql');
my @patches;
find(sub {
  if($_ =~ /^patch_${last_release}_${current_release}_\w+\.sql$/) {
    push(@patches, $File::Find::name);
  }
}, $sql_dir);

my $new_db_name = $db->create_db_name('schemapatchestemp');

SKIP: {
  
  skip 'Skipping DB patch tests as we cannot find the SQL at URL '.$cvs_url, (scalar(@patches)+1) unless defined $table_sql;
  
  my ($fh, $sql_schema_file) = tempfile();
  print $fh $table_sql;
  close $fh;
  
  #Create DB
  
  note 'Creating database '.$new_db_name;
  $dba->dbc()->do("create database $new_db_name");
  
  # Load SQL subroutine
  my $load_sql = sub {
    my ($path) = @_;
    my %args = ( host => $dbc->host(), port => $dbc->port(), user => $dbc->username(), password => $dbc->password());
    my $cmd_args = join(q{ }, map { "--${_}=$args{$_}" } keys %args);
    my $cmd = "mysql $cmd_args $new_db_name < $path 2>&1";
    my $output = `$cmd`;
    my $ec = ($? >> 8);
    if($ec != 0) {
      note($output);
      return fail("MySQL command failed with error code '$ec'");
    }
    return pass("MySQL was able to load the file $path core schema");
  };
  
  #Load last release's schema
  my $loaded_schema = $load_sql->($sql_schema_file);
  
  skip 'Skipping DB patch tests as we cannot find the SQL at URL '.$cvs_url, scalar(@patches) unless $loaded_schema;
  
  #Now apply all current patches
  foreach my $patch (@patches) {
    note "Applying patch $patch";
    $load_sql->($patch);
  }  
}

note 'Dropping database '.$new_db_name;
$dba->dbc()->do("drop database if exists $new_db_name");

done_testing();