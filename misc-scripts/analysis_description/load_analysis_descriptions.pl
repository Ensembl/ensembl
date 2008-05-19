#!/usr/local/ensembl/bin/perl -w

# POD documentation - main docs before the code

=pod

=head1 NAME

  load_analysis_descriptions.pl

=head1 SYNOPSIS

 script loads the analysis description table

=head1 DESCRIPTION

 The script reads the analysis description file also found in this directory
 analysis.descriptions and loads the descriptions which match the logic names in
 the analysis table

 It will warn about analyses present in the database which don't have descriptions
 in the file

=head1 OPTIONS

     Database options  

    -dbhost      host name for database (gets put as host= in locator)
    -dbport      For RDBs, what port to connect to (port= in locator)
    -dbname      For RDBs, what name to connect to (dbname= in locator)
    -dbuser      For RDBs, what username to connect as (dbuser= in locator)
    -dbpass      For RDBs, what password to use (dbpass= in locator)
    -description_file path to file containing descriptions. The file 
                      analysis.descriptions in this directory can be used and is an 
                      example of the format
    -help print out documentation

=head1 EXAMPLES

 perl load_analysis_descriptions.pl -dbhost my_host -dbuser user -dbpass ***** 
 -dbname my_db -description_file analysis.descriptions

=cut

use strict;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(warning throw);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Gene;

$! = 1;

my $dbhost = '';
my $dbuser;
my $dbpass;
my $dbport = 3306;
my $dbname = '';
my $file;
my $help = 0;

&GetOptions( 
            'host|dbhost=s'      => \$dbhost,
            'dbname=s'      => \$dbname,
            'user|dbuser=s'      => \$dbuser,
            'pass|dbpass=s'      => \$dbpass,
            'port|dbport=s'      => \$dbport,
            'file|descriptions=s' => \$file,
            'h|help!' => \$help,
           );

if(!$dbhost || !$dbname){
  print ("Need to pass in -dbhost $dbhost and -dbname $dbname\n");
  $help = 1;
}

if(!$file){
  $file = shift;
  if(!$file){
    $help = 1;
    print "Need to specify a description file on the commandline using -file\n";
  }
}

if($help){
  useage();
}



my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor
  (
   -host   => $dbhost,
   -user   => $dbuser,
   -dbname => $dbname,
   -pass   => $dbpass,
   -port   => $dbport,
  );


open(FH, $file) or throw("Failed to open $file $@");
my $aa = $db->get_AnalysisAdaptor;
my $analyses = $aa->fetch_all;
my %hash;
foreach my $analysis(@$analyses){
  $hash{lc($analysis->logic_name)} = $analysis;
}
LINE:while(<FH>){
  /^(\S+)\s+(\S+)\s+(.+)/ and do {
    my ($displayable, $logic_name, $description) = ($1, $2, $3);
    #print "Parsed ".$displayable." ".$logic_name." ".$description."\n";
    $description =~ s/^\s+//;
    $description =~ s/\s+$//;

    next if not $description;

    if (exists $hash{lc($logic_name)}) {
      my $analysis = $hash{lc($logic_name)};
      
      $analysis->description($description);
      $analysis->displayable($displayable);
      
      my $display_label = $logic_name;
      $display_label =~ s/_//g;
      $analysis->display_label($display_label);

      $aa->update($analysis);

      delete $hash{lc($logic_name)};
    }
  }
}

close(FH) or throw("Failed to close $file $@");

if ( scalar(keys %hash)==0) { 
  print "\nAll analysis descriptions have been updated, every analysis has a description now\n" ; 
}else{
  foreach my $k (keys %hash) {
    warning "No description was found for logic name $k ( ".$hash{$k}->dbID." ) \n";
  }
}


sub useage{
  exec('perldoc', $0);
  exit;
}
