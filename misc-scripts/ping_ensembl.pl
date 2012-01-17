#!/usr/bin/env perl

#######################################################
# This script is used to detect if you can contact the Ensembl database 
# server with your current setup. The program will attempt to print out
# helpful hints about how to resolve your problems. If they still persist
# then please contact helpdesk@ensembl.org.
#######################################################

use strict;
use warnings;

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

if($@) {
  print "ERROR: Cannot connect to Ensembl!\n";
  if($@ =~ /Can't locate Bio\/E/) {
    print "\tLooks like you need to setup your PERL5LIB with the Ensembl API. Please consult http://www.ensembl.org/info/docs/api/api_installation.html\n";
  }
  if($@ =~ /Can't locate Bio\/Perl/) {
    print "\tLooks like you need to setup your PERL5LIB with BioPerl. Please consult http://www.ensembl.org/info/docs/api/api_installation.html\n";
  }
  if($@ =~ /Cannot connect to/) {
    print "\tCannot seem to contact EnsemblDB at '$host' with the username '$user'. Try running 'ping $host' or asking your systems about firewalls against port $port\n";
  }
  if($@ =~ /internal name/) {
    print "\tSpecies was not found. You may have accidentally download the HEAD API version (found API release $db_version). Please consult http://www.ensembl.org/info/docs/api/api_installation.html\n";
  }
  if($@ =~ /Species not defined/) {
    print "\tSpecies was not found. You may have accidentally download the HEAD API version (found API release $db_version). Please consult http://www.ensembl.org/info/docs/api/api_installation.html\n";
  }
  print '='x80, "\n";
  print "If the problem persists please send the following error message to helpdesk\@ensembl.org\n";
  print $@;
  print '='x80, "\n";
  exit 1;
}