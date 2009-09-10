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
 the analysis table. Display labels are also set from this file.

 It will warn about analyses present in the database which don't have descriptions
 in the file.

 To not update analyses in the database you need to pass the -noupdate option.

=head1 OPTIONS

     Database options

    -dbhost      host name for database (gets put as host= in locator)
    -dbport      For RDBs, what port to connect to (port= in locator)
    -dbname      For RDBs, what name to connect to (dbname= in locator)
    -dbuser      For RDBs, what username to connect as (dbuser= in locator)
    -dbpass      For RDBs, what password to use (dbpass= in locator)
    -file        Path to file containing descriptions. The file 
                 analysis.descriptions in this directory can be used and is 
                 an example of the format. Multiple -file args can be specified
    -noupdate    Do not perform actual updates of analyses
    -pattern     check databases matching this PATTERN
                 Note that this is a database pattern of the form %core_53_%
    -help print out documentation

=head1 EXAMPLES

 perl load_analysis_descriptions.pl -dbhost my_host -dbuser user -dbpass ***** 
 -dbname my_db -file analysis.descriptions -file myanalysis.descriptions

if you want to update all databases for a type

perl load_analysis_descriptions.pl -dbhost my_host -dbuser user -dbpass ***** 
 -pattern '%_55_%' -file analysis.descriptions > & load_analysis_descriptions.log

syntax errors found as the definition file is parsed (usually the use of spaces rather
than tabs), cause the script to exit - fix each one and rerun.

=cut

use strict;
use Data::Dumper;
use Getopt::Long;
use DBI;

use Bio::EnsEMBL::Utils::Exception qw(warning throw);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Gene;

$| = 1;

my ($dsn,$dbh);

my $dbhost = '';
my $dbuser;
my $dbpass;
my $dbport = 3306;
my $dbname = '';
my @files = ();
my $noupdate;
my $help = 0;
my $pattern;


&GetOptions (
  'host|dbhost=s'       => \$dbhost,
  'dbname=s'            => \$dbname,
  'user|dbuser=s'       => \$dbuser,
  'pass|dbpass=s'       => \$dbpass,
  'port|dbport=s'       => \$dbport,
  'file|descriptions=s' => \@files,
  'noupdate'            => \$noupdate,
  'pattern=s'           => \$pattern,
  'h|help!'             => \$help
);

if (!$dbhost){
  print ("Need to pass a dbhost\n");
  $help =1;
}
if (!$dbname and !$pattern){
  print("Need to enter either a database name in -dbname or a pattern in -pattern\n");
  $help = 1;
}

unless(@files){
  @files = @_;
  unless(@files){
    $help = 1;
    print "Need to specify a description file on the commandline using -file\n";
  }
}

if($help){
  usage();
}

#connect to database
$dsn = "DBI:mysql:host=" . $dbhost . ";port=" . $dbport;

eval{
  $dbh = DBI->connect($dsn, $dbuser, $dbpass, 
		      {'RaiseError' => 1,
		       'PrintError' => 0});
};

# get all database names that match pattern
my ($sth, $sql);
my $sql_pattern = $pattern || $dbname;
$sql = "SHOW DATABASES LIKE '". $sql_pattern ."'";
$sth = $dbh->prepare($sql);
$sth->execute;

while (my ($dbname) = $sth->fetchrow_array){
  next unless $dbname =~ /core|cdna|otherfeatures/;
  next if $dbname =~ /coreexpression/;
  print "\n\nLooking at ... $dbname\n";
  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host   => $dbhost,
    -user   => $dbuser,
    -dbname => $dbname,
    -pass   => $dbpass,
    -port   => $dbport
  );

  # Pre-fetch all analyses in the database
  my $aa = $db->get_AnalysisAdaptor();
  my $analyses = $aa->fetch_all();
  my (%hash,%reference);
  foreach my $analysis(@$analyses){
    $hash{lc($analysis->logic_name())} = $analysis;
  }

# Parse the description files and check formatting
  foreach my $file( @files ){
    open(FH, $file) or throw("Failed to open $file $@");
	
  LINE: while(my $row = <FH>){
      chomp($row);
      next if ($row =~ /^\#/);   # skip comments
      next if ($row =~ /^$/);    # and blank lines
      next if ($row =~ /^\s+$/); # and whitespace-only lines

      my ($nr, $logic_name, $description, $display_label, $displayable, $web_data) = split(/\t/, $row);
#	  print join("\t", $logic_name, $description, $display_label, $displayable, $web_data), "\n";

      unless ($logic_name && defined($displayable)) {
	throw("Please check description file entry for logic_name $logic_name (" . join("\n", $logic_name, $description, $display_label, $displayable, $web_data) . ")");
	exit;
      }
      unless (defined $displayable){
	throw("In the analysis_description file, logic name '$logic_name' should contain, at least, 5 columns: Number, logic_name, description, display_label and displayable. Fix it !!");
      }
      unless ($displayable =~ m/^[01]$/) {
	throw("Displayable flag for analysis '$logic_name' has to be either 0 or 1, but not '$displayable'!");
      }

      $reference{lc($logic_name)} = {
	nr            => $nr,
	description   => $description   || '',
	display_label => $display_label || '',
	displayable   => $displayable   || '', 
	web_data      => $web_data      || '',
      };
	  
      $description =~ s/^\s+//;
      $description =~ s/\s+$//;
	  
      next if not $description;
	  
      if (exists $hash{lc($logic_name)}) {
	      
	my $analysis = $hash{lc($logic_name)};
	      
	$analysis->description($description);
	$analysis->displayable($displayable);
	$analysis->display_label($display_label);
	$web_data ? $analysis->web_data($aa->get_dumped_data($web_data)) : $analysis->{_web_data} = undef;
#	print Dumper  $analysis->web_data();
	      
	unless ( $noupdate ) {
	  $aa->update($analysis) ; 
	}
	      
	delete $hash{lc($logic_name)};
      }
    }
    close(FH) or throw("Failed to close $file $@");
  }
    
  if ( scalar(keys %hash)==0) {
    unless  ($noupdate) {
      print STDERR "\nAll analysis descriptions have been updated, every analysis has a description now\n";
    } else {
      print STDERR "\nEvery analysis has a description in the file, all analysis descriptions can be updated.\n".
	"To write analysis descriptions to the analysis_description table in your DB,\n".
	  "please run this script excluding the -noupdate option on the commandline\n";
    }
  } 
  else {
    foreach my $ln (keys %hash) {
      unless (exists $reference{$ln}) {
	warning ("[$dbname]: Analysis '$ln' doesn't exist in reference file(s) '"
		   . join( "','", @files )
		     . "'! It needs to be added first");
      }
    }
  }
}

sub usage{
  exec('perldoc', $0);
  exit;
}
