#!/usr/local/ensembl/bin/perl -w

# POD documentation - main docs before the code

=pod

=head1 NAME

 update_misc_set_descripitons.pl 

=head1 SYNOPSIS

 script loads the misc_set descriptions in the misc_set table

=head1 DESCRIPTION

 The script reads the misc_set description file also found in this directory
 misc_set.descriptions and loads the descriptions which match the names in
 the misc_set table.

 It will warn about names present in the database which don't have descriptions
 in the file.

 To not update descriptions in the database you need to pass the -noupdate option.

=head1 OPTIONS

     Database options

    -dbhost      host name for database (gets put as host= in locator)
    -dbport      For RDBs, what port to connect to (port= in locator)
    -dbname      For RDBs, what name to connect to (dbname= in locator)
    -dbuser      For RDBs, what username to connect as (dbuser= in locator)
    -dbpass      For RDBs, what password to use (dbpass= in locator)
    -file        Path to file containing descriptions.
    -noupdate    Do not perform actual updates of analyses
    -help print out documentation

=head1 EXAMPLES

 perl update_misc_set_descripitons.pl -dbhost my_host -dbuser user -dbpass ***** 
 -dbname my_db -file misc_set.descriptions

=cut

use strict;
use Data::Dumper;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(warning throw);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Gene;

$| = 1;

my $dbhost = '';
my $dbuser;
my $dbpass;
my $dbport = 3306;
my $dbname = '';
my @files = ();
my $noupdate;
my $help = 0;

&GetOptions (
	'host|dbhost=s'       => \$dbhost,
	'dbname=s'            => \$dbname,
	'user|dbuser=s'       => \$dbuser,
	'pass|dbpass=s'       => \$dbpass,
	'port|dbport=s'       => \$dbport,
	'file|descriptions=s' => \@files,
	'noupdate'              => \$noupdate,
	'h|help!'             => \$help
             );

if(!$dbhost || !$dbname){
  print ("Need to pass in -dbhost $dbhost and -dbname $dbname\n");
  $help = 1;
}

unless(@files){
  @files = @_;
  unless(@files){
    $help = 1;
    print "Need to specify a description file on the command line using -file\n";
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

# Pre-fetch all misc_set names in the database
my $msa = $db->get_MiscSetAdaptor;;
my $misc_sets = $msa->fetch_all();

my (%hash,%reference);
foreach my $misc_set(@$misc_sets){
  $hash{lc($misc_set->name())} = $misc_set;
}

# Parse the description files
foreach my $file( @files ){
  open(FH, $file) or throw("Failed to open $file $@");

 LINE: while(my $row = <FH>){
    
    chomp($row);
    
    next if ($row =~ /^\#/);   # skip comments
    next if ($row =~ /^$/);    # and blank lines
    next if ($row =~ /^\s+$/); # and whitespace-only lines
    
    my ($name, $description) = split(/\t/, $row);
    #print join("\t", $logic_name, $description, $display_label, $displayable, $web_data), "\n";
    
    $reference{lc($name)} = {
      description   => $description,
    };
    
    $description =~ s/^\s+//;
    $description =~ s/\s+$//;
    
    next if not $description;
    
    if (exists $hash{lc($name)}) {
	
        my $misc_set = $hash{lc($name)};
        
        $misc_set->description($description);

        unless ( $noupdate ) { 
          $msa->update($misc_set) ; 
        }
        
        delete $hash{lc($name)};
      }
  }
  close(FH) or throw("Failed to close $file $@");
}

if ( scalar(keys %hash)==0) {
	unless  ($noupdate) {
		print STDERR "\nAll misc_set descriptions have been updated, every misc_set has a description now\n";
	} else {
		print STDERR "\nEvery misc_set has a description in the file, all misc_set descriptions can be updated.\n".
                                "To write misc_set descriptions to the misc_set table in your DB,\n".
                                "please run this script including the -update option on the commandline\n";
	}
} else {
    foreach my $ln (keys %hash) {
        warning ("Misc_set '$ln' doesn't exist in reference file(s) '"
               . join( "','", @files )
               . "'! It needs to be added first")
			unless (exists $reference{$ln});
        warning "[$dbname] No description was found for name '$ln':\n";

    }
}


sub useage{
  exec('perldoc', $0);
  exit;
}
