#!/usr/bin/env perl

# Command line front end to ExportProgramWriter.pm.
#
# Use --help for help.

use ExportProgramWriter;
use Carp;
use Getopt::Long;

# Command line param variables. e.g. --newtable value stored in $newtable
my $verbose;
my $help;
my $outdir;
my $newtable;
my $oldtable;
my $print;

# Created during compare_two_schema()
my @unchanged = [];
my @changed = [];
my @novel = [];


main();

sub main {
  GetOptions(
	     "verbose" => \$verbose,
	     "help" => \$help,
	     "outdir:s" => \$outdir,
	     "newtable:s" => \$newtable,
	     "oldtable:s" => \$oldtable,
	     "print" => \$print
	    );
  
  if (scalar($help || @ARGV)==0) {
    print "usage: $0 <--print> <--verbose> <--newtable NEW_TABLE_NAME <--oldtable OLD_TABLE_NAME>> <--outdir OUTPUT_DIR> new_schema_filename old_schema_filename \n";
    print "usage: $0 --help\n";
    exit 1;
    
  }
  
  elsif ($newtable) {
    if ( !$oldtable ) { $oldtable = $newtable; }
    create_specific_export_program($newtable, @ARGV[0], $oldtable, @ARGV[1]);
  }

  elsif (scalar(@ARGV)==2) {
    compare_two_schemas(@ARGV[1], @ARGV[0]);
    if ( !$print ) {
      create_many_export_programs("Creating export program for Changed table : ", @changed);
      create_many_export_programs("Creating export program for Novel table : ", @novel);
    }
  }
  

}


sub file_to_string {
  (my $filename) = @_;
  my $str = "";
  open(INFILE, $filename)|| confess "cannot open \"$filename\": $!";
  $str = join("", <INFILE>);
  close(INFILE);
  return $str;
}


sub save {
  (my $package_name, 
   my $new_create_table, 
   my $old_create_table) =@_;

  my $filename = $package_name;
  $filename =~ s/::/\//g;
  $filename .= ".pm";
  $filename = $outdir . "/" . $filename;

  open(OUTFILE, ">" . $filename) || die "failed to open file $filename : $!";
  print OUTFILE ExportProgramWriter::create_export_program($package_name, 
							   $new_create_table, 
							   $old_create_table);
  close(OUTFILE);
}


sub create_many_export_programs {
  (my $description, my @new_old_pairs) = @_;

  foreach $pair (@new_old_pairs) {
    my $new_create_table = @{$pair}[0];
    my $old_create_table = @{$pair}[1];
    my $package_name = ucfirst($new_create_table->get_table_name());
    if ($verbose) {print $description . $new_create_table->get_table_name() . "\n"};
    save($package_name, $new_create_table, $old_create_table);
  }
}


sub create_specific_export_program {
  ( my $new_table_name,
    my $new_schema_file,
    my $old_table_name,
    my $old_schema_file, ) = @_;

  my $package_name = ucfirst($new_table_name);
  
  my $new_schema = file_to_string($new_schema_file);
  my $old_schema =  file_to_string($old_schema_file);
  
  my $new_create_table = CreateTable::create_from_string($new_table_name, $new_schema);
  my $old_create_table = CreateTable::create_from_string($old_table_name, $old_schema);
  if ($old_create_table==0) {
    $old_create_table = CreateTable::new();
  }
  save($package_name, $new_create_table, $old_create_table);
}



sub compare_two_schemas {
  (my $new_schema_filename, my $old_schema_filename) = @_;
  
  # Read in new schema file and create CreateTable objects for each..
  
  my $new_schema = file_to_string($new_schema_filename);
  my $old_schema = file_to_string($old_schema_filename);
  my @create_tables = CreateTable::create_all_from_string($new_schema);
  
  # Identify different tables.
  my $len = scalar(@create_tables);
  
  my $unchanged_index = 0;    
  my $changed_index = 0;
  my $novel_index = 0;

  my $empty_create_table = CreateTable::new();
  $empty_create_table->set_table_name("****UNKNOWN*****");
  my @empty = [];
  $empty_create_table->set_fields(\@empty);
  
  print "EMPTY TABLE : " . $empty_create_table->to_string() . "\n";

  for(my $i=0; $i<$len; ++$i) {
    # look for a correspondingly named table in the old file
    my $new_create_table = @create_tables[$i];
    my $new_table_name = $new_create_table->get_table_name();
    my $equivalent_old_create_table = CreateTable::create_from_string( $new_table_name, $old_schema);
    
    if ( $verbose ) {print "$i/$len " . $new_create_table->get_table_name() . " ";}
    
    if ( !$equivalent_old_create_table ) {
      #@novel[$novel_index++] = $new_create_table;
      @novel[$novel_index++] = ([$new_create_table, $empty_create_table]);
      if ( $verbose ) {print "UNKOWN\n";}
    }
    elsif ( $new_create_table->equals($equivalent_old_create_table) ) {
      #@unchanged[$unchanged_index++] = $new_create_table;
      @unchanged[$unchanged_index++] = ([$new_create_table, $equivalent_old_create_table]);
      if ( $verbose ) {print "SAME\n" . join(",", @{$new_create_table->get_fields()}) . "\n" . join(",", @{$equivalent_old_create_table->get_fields()}) . "\n";}
    } 
    else {
      #@changed[$changed_index++] = $new_create_table;
      @changed[$changed_index++] = ([$new_create_table, $equivalent_old_create_table]);
      if ( $verbose ) {print "CHANGED\n" . join(",", @{$new_create_table->get_fields()}) . "\n" . join(",", @{$equivalent_old_create_table->get_fields()}) . "\n";}
    }
  }
  
  print "Unchanged : " . scalar(@unchanged) ."\n";
  print "Changed : " . scalar(@changed) ."\n";
  print "Novel : " . scalar(@novel) ."\n";

  
}
