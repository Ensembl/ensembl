#!/usr/local/bin/perl -w

use strict;
use warnings;

use lib 't';

use Getopt::Std;
use Test::Harness;
use MultiTestDB;

use vars qw($opt_l $opt_h);

#read command line options
&usage unless getopts('lh');


#print usage on '-h' command line option
if($opt_h) {
  &usage;
  exit;
}

#list test files on '-l' command line option
if($opt_l) {
  foreach my $file (@{&get_all_tests('.', \@ARGV )}) {
    print "$file\n";
  }
  exit;
}


#set environment var
$ENV{'RUNTESTS_HARNESS'} = 1;

#make sure proper cleanup is done if the user interrupts the tests
$SIG{HUP} = $SIG{KILL} = $SIG{INT} = 
  sub {warn "\n\nINTERRUPT SIGNAL RECEIEVED\n\n"; &clean;};

#create a multitest db, its destruction will clean up after scripts
my $clean_up = new MultiTestDB;

#run all specified tests
eval {
  runtests(@{&get_all_tests('.', \@ARGV)});
};

&clean;

sub clean {
  #unset env var indicating final cleanup should be performed
  delete $ENV{"RUNTESTS_HARNESS"};
  exit;
}

=head2 get_all_tests

  Arg [1]    : string $dir
               the name of the directory retrieve a list of tests from
  Arg [2]    : (optional) listref $input_files 
               testfiles or directories to retrieve. If not specified all 
               ".t" files in $dir are taken.
  Example    : @test_files = read_test_dir('t');
  Description: Returns a list of testfiles in the directories specified by
               the @tests argument.  The relative path is given as well as
               with the testnames returned.  Only files ending with .t are
               returned.  Subdirectories are recursively entered and the test
               files returned within them are returned as well.
  Returntype : listref of strings.
  Exceptions : none
  Caller     : general

=cut

sub get_all_tests {
  my ($dir, $input_files) = @_;
  
  my @files;
  my @out = ();
  local *DIR;
  
  unless(opendir(DIR, $dir)) {
    warn("WARNING: cannot open directory $dir\n");
    return [];
  }

  if($input_files && @$input_files) {
    #input files were specified so use them
    @files = @$input_files; 
  } else {
    #otherwise use every file in the directory
    @files = readdir DIR;
  }     
     
  #filter out CVS files, files beginning with '.' and files ending in ~
  @files = grep !/(^\.)|(^CVS$)|(~$)/, @files;

  foreach my $file (@files) {
    $file = "$dir/$file";

    if(-d $file) {
      #do a recursive call on directories
      push @out, @{get_all_tests("$file")};
    } elsif ($file =~ /\.t$/) {
      #files ending with a '.t' are considered test files
      unless(-r $file && -f $file) {
	warn("WARNING: cannot read test file $file\n");
      }
      push @out, $file;
    } 
  }
  
  closedir DIR;

  return \@out;
}
  


sub usage {
  print "usage:\n";
  print "\tlist tests:        run_tests.pl -l [<testfiles or dirs> ...]\n";
  print "\trun tests:         run_tests.pl [<testfiles or dirs> ...]\n";
} 


