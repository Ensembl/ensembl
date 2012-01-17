#!/usr/bin/env perl

use Bio::EnsEMBL::Registry;
eval {
  Bio::EnsEMBL::Registry->load_registry_from_db(
    -host    => 'ensembldb.ensembl.org',
    -user    => 'anonymous'
  );
  my $human = Bio::EnsEMBL::Registry->get_DBAdaptor();
  my $name = $human->get_MetaContainer()->get_scientific_name();
  if('Homo sapiens' == $name) {
    print "Connection to Ensembl works and you can query the human core database\n";
  }
};

if($@) {
  print "ERROR: Cannot connect to Ensembl!\n";
  print "Please send the following error message to helpdesk\@ensembl.org\n";
  print $@;
}