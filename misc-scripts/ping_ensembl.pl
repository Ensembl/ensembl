#!/usr/bin/env perl

#######################################################
# This script is used to detect if you can contact the Ensembl database 
# server with your current setup. The program will attempt to print out
# helpful hints about how to resolve your problems. If they still persist
# then please contact helpdesk@ensembl.org.
#######################################################

use strict;
use warnings;

use File::Temp qw/tempfile/;
use Net::FTP;

my $host = 'ensembldb.ensembl.org';
my $user = 'anonymous';
my $port = 5306;
my $db_version = '-';

eval {
  require Bio::EnsEMBL::Registry;
  require Bio::EnsEMBL::ApiVersion;
  require Bio::Perl;
  Bio::EnsEMBL::Registry->load_registry_from_db(
    -host    => $host,
    -user    => $user
  );
  $db_version = Bio::EnsEMBL::ApiVersion::software_version();
  my $human = Bio::EnsEMBL::Registry->get_DBAdaptor('homo_sapiens', 'core');
  my $name = $human->get_MetaContainer()->get_scientific_name();
  if('Homo sapiens' eq $name) {
    print "Installation is good. Connection to Ensembl works and you can query the human core database\n";
  }
  else {
    print "Installation is good. Connection to Ensembl works but the species name '$name' was not as expected. Please contact the developers at dev\@ensembl.org to update this assertion\n";
  }
};
my $error = $@;

#If no error found then see if we've got all of our external modules available
if(!$error) {
  $error = '';
  eval {
    foreach my $module (qw/Compara Variation Funcgen/) {
      my $full_module = "Bio::EnsEMBL::${module}::DBSQL::DBAdaptor";
      eval "require $full_module;";
      if($@) {
        $error .= "\tMissing the checkout $module\n"; 
      }
    }
  };
}

#Check the current release of datafiles from the FTP site
my $ftp_version = -1;
eval {
  my $ftp = Net::FTP->new('ftp.ensembl.org', Debug => 0);
  $ftp->login("anonymous",'-anonymous@');
  $ftp->cwd('/pub');
  my ($fh, $filename) = tempfile();
  close($fh);
  $ftp->get('current_README', $filename);
  $ftp->quit();
  open($fh, '<', $filename);
  local $/ = undef;
  my $ftp_readme = <$fh>;
  close($fh);
  ($ftp_version) = $ftp_readme =~ /Ensembl Release (\d+) Databases/;
};

#Print all the errors which could have occured 
if($error) {
  print "ERROR: Error detcted when connecting to Ensembl!\n";
  if($error =~ /Can't locate Bio\/E/) {
    print "\tLooks like you need to setup your PERL5LIB with the Ensembl API. Please consult http://www.ensembl.org/info/docs/api/api_installation.html\n";
  }
  if($error =~ /Can't locate Bio\/Perl/) {
    print "\tLooks like you need to setup your PERL5LIB with BioPerl. Please consult http://www.ensembl.org/info/docs/api/api_installation.html\n";
  }
  if($error =~ /Cannot connect to/) {
    print "\tCannot seem to contact EnsemblDB at '$host' with the username '$user'. Try running 'ping $host' or asking your systems about firewalls against port $port\n";
  }
  if($error =~ /internal name/) {
    print "\tSpecies was not found. You may have accidentally download the HEAD API version (found API release $db_version & public FTP release is $ftp_version). Please consult http://www.ensembl.org/info/docs/api/api_installation.html\n";
  }
  if($error =~ /Species not defined/) {
    print "\tSpecies was not found. You may have accidentally download the HEAD API version (found API release $db_version & public FTP release is $ftp_version). Please consult http://www.ensembl.org/info/docs/api/api_installation.html\n";
  }
  if($error =~ /Missing the checkout/) {
    print "\tYour core installation was good but supplementary modules cannot be found. If you wish to access these other Ensembl resources add the libraries to your PERL5LIB:\n";
    print $error;
    #bail early
    exit 0;
  }
  print '='x80, "\n";
  print "If the problem persists please send the following error message to helpdesk\@ensembl.org\n";
  print $error;
  print '='x80, "\n";
  exit 1;
}