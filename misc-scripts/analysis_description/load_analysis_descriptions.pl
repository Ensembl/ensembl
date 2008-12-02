#!/opt/local/bin/perl -w
###!/usr/local/ensembl/bin/perl -w

# POD documentation - main docs before the code

=pod

=head1 NAME

  load_analysis_descriptions.pl

=head1 SYNOPSIS

 script loads the analysis description table

=head1 DESCRIPTION

 The script reads the analysis description file also found in this directory
 analysis.descriptions and loads the descriptions which match the logic names in
 the analysis table. Display labels are also set from this file.

 It will warn about analyses present in the database which don't have descriptions
 in the file.

 To actually update analyses in the database you need to pass the update option.

=head1 OPTIONS

     Database options

    -dbhost      host name for database (gets put as host= in locator)
    -dbport      For RDBs, what port to connect to (port= in locator)
    -dbname      For RDBs, what name to connect to (dbname= in locator)
    -dbuser      For RDBs, what username to connect as (dbuser= in locator)
    -dbpass      For RDBs, what password to use (dbpass= in locator)
    -file        Path to file containing descriptions. The file 
                      analysis.descriptions in this directory can be used and is an 
                      example of the format
	-update      Perform actual updates of analyses
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
my ($file, $update);
my $help = 0;

&GetOptions (
	'host|dbhost=s'       => \$dbhost,
	'dbname=s'            => \$dbname,
	'user|dbuser=s'       => \$dbuser,
	'pass|dbpass=s'       => \$dbpass,
	'port|dbport=s'       => \$dbport,
	'file|descriptions=s' => \$file,
	'update'              => \$update,
	'h|help!'             => \$help
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

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
                                            -host   => $dbhost,
                                            -user   => $dbuser,
                                            -dbname => $dbname,
                                            -pass   => $dbpass,
                                            -port   => $dbport
                                            );

open(FH, $file) or throw("Failed to open $file $@");
my $aa = $db->get_AnalysisAdaptor();
my $analyses = $aa->fetch_all();
my %hash;
foreach my $analysis(@$analyses){
    $hash{lc($analysis->logic_name())} = $analysis;
}

LINE: while(my $row = <FH>){
    
    chomp($row);
    
    next if ($row =~ /^\#/);   # skip comments
    next if ($row =~ /^$/);    # and blank lines
    next if ($row =~ /^\s+$/); # and whitespace-only lines
    
    my ($nr, $logic_name, $description, $display_label, $displayable, $web_data) = split(/\t/, $row);
    #print join("\t", $logic_name, $description, $display_label, $displayable, $web_data), "\n";
    
    $description =~ s/^\s+//;
    $description =~ s/\s+$//;
    
    next if not $description;
    
    if (exists $hash{lc($logic_name)}) {

        my $analysis = $hash{lc($logic_name)};
        
        $analysis->description($description);
        $analysis->displayable($displayable);
        $analysis->display_label($display_label);
        $analysis->web_data($web_data) if $web_data;
        
        $aa->update($analysis) if $update;
        
        delete $hash{lc($logic_name)};
    }
}

close(FH) or throw("Failed to close $file $@");

if ( scalar(keys %hash)==0) {
	if ($update) {
		print STDERR "\nAll analysis descriptions have been updated, every analysis has a description now\n";
	} else {
		print STDERR "\nAll analysis descriptions can been updated, every analysis has a description in the file\n";
	}
}else{
    foreach my $k (keys %hash) {
        
        warning "[$dbname] No description was found for logic name $k ( dbID=".$hash{$k}->dbID."; displayable=".$hash{$k}->displayable." ) \n";

    }
}


sub useage{
  exec('perldoc', $0);
  exit;
}
