=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

package XrefParser::ArrayExpressParser;

## Parsing format looks like (so we extract the species name):
# anopheles_gambiae.A-AFFY-102.tsv
#

use strict;
use warnings;
use Carp;
use base qw( XrefParser::BaseParser );
use Bio::EnsEMBL::Registry;
use Net::FTP;

my $default_ftp_server = 'ftp.ebi.ac.uk';
my $default_ftp_dir = 'pub/databases/microarray/data/atlas/bioentity_properties/ensembl';

sub run_script {

  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $file         = $ref_arg->{file};
  my $verbose      = $ref_arg->{verbose};

  if((!defined $source_id) or (!defined $species_id) or (!defined $file) ){
    croak "Need to pass source_id, species_id and file as pairs";
  }
  $verbose |=0;

  my $project;
  my $user  ="ensro";
  my $host;
  my $port = 3306;
  my $dbname;
  my $pass;

  if($file =~ /project[=][>](\S+?)[,]/){
    $project = $1;
  }
  if($file =~ /host[=][>](\S+?)[,]/){
    $host = $1;
  }
  if($file =~ /port[=][>](\S+?)[,]/){
    $port =  $1;
  }
  if($file =~ /dbname[=][>](\S+?)[,]/){
    $dbname = $1;
  }
  if($file =~ /pass[=][>](\S+?)[,]/){
    $pass = $1;
  }
  if($file =~ /user[=][>](\S+?)[,]/){
    $user = $1;
  }

  my %species_id_to_names = $self->species_id2name();
  my $species_id_to_names = \%species_id_to_names;
  my $names = $species_id_to_names->{$species_id};
  my $species_lookup = $self->_get_species($verbose);
  my $active = $self->_is_active($species_lookup, $names, $verbose);
 
  if (!$active) {
      return;
  }

  my $species_name = $species_id_to_names{$species_id}[0];

  #get stable_ids from core and create xrefs 

  my $registry = "Bio::EnsEMBL::Registry";
  my ($gene_adaptor);
  if ($host) {
    my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
      '-host'     => $host,
      '-port'     => $port,
      '-user'     => $user,
      '-pass'     => $pass,
      '-dbname'   => $dbname,
      '-species'  => $species_name,
      '-group'    => 'core',
        );
    $gene_adaptor = $db->get_GeneAdaptor();
  } elsif ($project eq 'ensembl') {
    print "Loading the Registry\n" if $verbose;
    $registry->load_registry_from_multiple_dbs( 
      {
        '-host'    => 'ens-staging1',
        '-user'    => 'ensro',
      },
      {
        '-host'     => 'ens-staging2',
        '-user'     => 'ensro',
      }
        );
    $gene_adaptor = $registry->get_adaptor($species_name, 'core', 'Gene');
  } elsif ($project eq 'ensemblgenomes') {
    $registry->load_registry_from_multiple_dbs( 
      {
        '-host'     => 'mysql-eg-staging-1.ebi.ac.uk',
        '-port'     => 4160,
        '-user'     => 'ensro',
      },
      {
        '-host'     => 'mysql-eg-staging-2.ebi.ac.uk',
        '-port'     => 4275,
        '-user'     => 'ensro',
      },
        );
    $gene_adaptor = $registry->get_adaptor($species_name, 'core', 'Gene');
  } else {
      die("Missing or unsupported project value. Supported values: ensembl, ensemblgenomes");
  }
  print "Finished loading the registry\n" if $verbose;

  my @stable_ids = map { $_->stable_id } @{$gene_adaptor->fetch_all()};

  my $xref_count = 0;
  foreach my $gene_stable_id (@stable_ids) {
      
      my $xref_id = $self->add_xref({ acc        => $gene_stable_id,
					 label      => $gene_stable_id,
					 source_id  => $source_id,
					 species_id => $species_id,
					 info_type => "DIRECT"} );
	
      $self->add_direct_xref( $xref_id, $gene_stable_id, 'gene', '');
      if ($xref_id) {
	 $xref_count++;
      }
  }
	   
  print "Added $xref_count DIRECT xrefs\n" if($verbose);
  if ( !$xref_count ) {
      return 1;    # 1 error
  }

  return 0; # successfull	   

}

sub _get_species {
  my ($self, $verbose) = @_;
  $verbose = (defined $verbose) ? $verbose : 0;
  
  my $ftp = Net::FTP->new($default_ftp_server, Debug => $verbose) or confess "Cannot connect to $default_ftp_server: $@";
  $ftp->login("anonymous",'-anonymous@') or confess "Cannot login ", $ftp->message;
  $ftp->cwd($default_ftp_dir);
  my @files = $ftp->ls() or confess "Cannot change to $default_ftp_dir: $@";
  $ftp->quit;

  my %species_lookup;
  foreach my $file (@files) {
    my ($species) = split(/\./, $file);
    $species_lookup{$species} = 1;
  }
  return \%species_lookup;
}

sub _is_active {
  my ($self, $species_lookup, $names, $verbose) = @_;
  #Loop through the names and aliases first. If we get a hit then great
  my $active = 0;
  foreach my $name (@{$names}) {
    if($species_lookup->{$name}) {
      printf('Found ArrayExpress has declared the name "%s". This was an alias'."\n", $name) if $verbose;
      $active = 1;
      last;
    }
  }
  return $active;
}


1;
