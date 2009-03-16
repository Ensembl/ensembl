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

 To actually update analyses in the database you need to pass the update option.

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
    -update      Perform actual updates of analyses
    -help print out documentation

=head1 EXAMPLES

 perl load_analysis_descriptions.pl -dbhost my_host -dbuser user -dbpass ***** 
 -dbname my_db -file analysis.descriptions -file myanalysis.descriptions -update

=cut

use strict;
use Data::Dumper;
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
my @files = ();
my $update;
my $help = 0;

&GetOptions (
	'host|dbhost=s'       => \$dbhost,
	'dbname=s'            => \$dbname,
	'user|dbuser=s'       => \$dbuser,
	'pass|dbpass=s'       => \$dbpass,
	'port|dbport=s'       => \$dbport,
	'file|descriptions=s' => \@files,
	'update'              => \$update,
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

# Pre-fetch all analyses in the database
my $aa = $db->get_AnalysisAdaptor();
my $analyses = $aa->fetch_all();
my (%hash,%reference);
foreach my $analysis(@$analyses){
  $hash{lc($analysis->logic_name())} = $analysis;
}

# Parse the description files
foreach my $file( @files ){
  open(FH, $file) or throw("Failed to open $file $@");

 LINE: while(my $row = <FH>){
    
    chomp($row);
    
    next if ($row =~ /^\#/);   # skip comments
    next if ($row =~ /^$/);    # and blank lines
    next if ($row =~ /^\s+$/); # and whitespace-only lines
    
    my ($nr, $logic_name, $description, $display_label, $displayable, $web_data) = split(/\t/, $row);
    #print join("\t", $logic_name, $description, $display_label, $displayable, $web_data), "\n";
    
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

        $aa->update($analysis) if $update;
        
        delete $hash{lc($logic_name)};
      }
  }
  close(FH) or throw("Failed to close $file $@");
}

if ( scalar(keys %hash)==0) {
	if ($update) {
		print STDERR "\nAll analysis descriptions have been updated, every analysis has a description now\n";
	} else {
		print STDERR "\nEvery analysis has a description in the file, all analysis descriptions can be updated.\n".
                                "To write analysis descriptions to the analysis_description table in your DB,\n".
                                "please run this script including the -update option on the commandline\n";
	}
} else {
    foreach my $ln (keys %hash) {
        warning ("Analysis '$ln' doesn't exist in reference file(s) '"
               . join( "','", @files )
               . "'! It needs to be added first")
			unless (exists $reference{$ln});
        warning "[$dbname] No description was found for logic_name '$ln':\n".
			"\tref:\t display_label='".$reference{$ln}{display_label}."'; displayable=".$reference{$ln}{displayable}."; nr=".$reference{$ln}{nr}."\n".
			"\tdb: \t display_label='".$hash{$ln}->display_label."'; displayable=".$hash{$ln}->displayable."; dbID=".$hash{$ln}->dbID."\n";

    }
}


sub useage{
  exec('perldoc', $0);
  exit;
}
