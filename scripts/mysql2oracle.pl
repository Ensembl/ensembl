#!/usr/bin/perl -w

use strict;
use vars qw/$DEBUG %output %types $filename $outfile $loading/;
use Getopt::Long;

$DEBUG = 0;
%output = ();
%types = ();

GetOptions('infile=s' =>\$filename,
	   'outfile=s'=>\$outfile,
	   'load'=>\$loading);

defined $filename or die "MySQL Ensemble table file not given...\n";
defined $outfile or $outfile = "oratable.sql";

open(SQL, "<$filename");
while (<SQL>) {
    chomp;
    /^ *\#/ and next;
    # bugger, synonym is not a valid oracle column name
    /create table/i and parseCreateTable(\*SQL ,$_);
}

buildSchema();


sub buildSchema {
    open (SCHEMA, ">$outfile");
    open (CONSTRAINTS, ">$outfile.con") if defined $loading;
    
    while (my ($table, $rh_definition) = each %output) {
	print SCHEMA "\nCREATE TABLE ", $table, " (";
	my $columnDefinition = " ";
	foreach my $column (@{$rh_definition->{'columns'}}) {
	    $columnDefinition .= "\n\t" . $column . ",";
	}
	# last , needs removing
	$columnDefinition =~ s|,$|\n\);\n|;
	print SCHEMA $columnDefinition;
	
	#output the primary key data
	defined $output{$table}->{'primarykey'} and do {
	    my $out = "\nalter table " . $table . " add constraint " . $table . "_pk ".  $output{$table}->{'primarykey'} . ";\n";
	    defined $loading ? print CONSTRAINTS $out : print SCHEMA $out;
	};
	
	
	# output the normal key data
	defined $output{$table}->{'key'} and do {
	    foreach my $rl_key (@{$output{$table}->{'key'}}) {
		my $out =  "\ncreate index  " . $rl_key->[0] . " on " . $table . " " . $rl_key->[1] . ";\n";
		defined $loading ? print CONSTRAINTS $out : print SCHEMA $out;
	    }
	};
	defined $output{$table}->{'unique'} and do {
	    foreach my $rl_unique (@{$output{$table}->{'unique'}}) {
		my $out =  "\ncreate unique index " . $rl_unique->[0] . " on " . $table . " " . $rl_unique->[1] . ";\n";
		defined $loading ? print CONSTRAINTS $out : print SCHEMA $out;
	    }
	};
    }
}
	    


sub parseCreateTable {
    my ($fh, $firstLine) = @_;
    my $tableName;
    
    if ($firstLine =~ /(create table) +(\w+) *\(/i) {
	$tableName = $2;
    } else {
	die "Line does not match create table: $firstLine";
    }
    while (<$fh>) {
	chomp;
	$_ =~ s/\s+//i;
	$_ =~ s/synonym/synonymname /gi;
	# return if it matches close bracket
	/\);/ and return;	
	
	# there is a problem with the MAX_ROWS stuff. Look at it later
	/\) *MAX_ROWS.+;/ and return;
	# skip empty lines
	next if /^ *$/;
	next if /^ *\#/;

	# process table keys and constraints
	/primary key/i and do {parsePrimaryKey($tableName, $_); next};
	/^ *key/i and do {parseKey($tableName, $_);next};
	/unique/i and do {parseUniqueConstraint($tableName, $_); next};
	parseColumnDefinition($tableName, $_);

    }
}

sub parsePrimaryKey{
    my ($table, $def) = @_;
    $def =~ s/^ *//;
    $def =~ s/, *?$//;
    $DEBUG and print "\n--->DEBUG: $def\n";
    $output{$table}->{'primarykey'} = $def;
}

sub parseKey {
    my ($table, $def) = @_;
    $def =~ s/^ *//;
    $def =~ s/, *$//;
    if ($def =~ /(KEY) +(\w+) *(\(.+\))/i) {
	push @{$output{$table}->{'key'}}, [$table . $2, $3];
    } elsif ($def =~ /(KEY) *(\(.+\))/i) {
	my $index = 1;
	$index = scalar @{$output{$table}->{'key'}} + 1 if 
	    defined $output{$table}->{'key'};
	push @{$output{$table}->{'key'}}, [$table . "_key_" . $index, $2];
    } else {
	die "Unknown key definition: ", $def;
    }
}

sub parseUniqueConstraint{
    my ($table, $def) = @_;

    $def =~ s/^ *//;
    $def =~ s/, *$//;

    my $index = 1;
    if ($def =~ /UNIQUE KEY (.*?) *(\(.+?\))/i) {
      push  @{$output{$table}->{'unique'}}, [$table.$1, $2];
    } elsif ($def =~ /UNIQUE +(.*?) *(\(.+?\))/i) {
      push @{$output{$table}->{'unique'}}, [$table.$1, $2];
    } elsif ($def =~ /UNIQUE *(\(.+?\))/i) {
      push  @{$output{$table}->{'unique'}}, [$table . "_" . $index, $1];
    } else {
      die "unique constraint doesn\'t parse: ", $def;
    }
    
  }

sub parseColumnDefinition {
    my ($table, $def) = @_;

    # throw out , if it is in there
    $def =~ s/,\s*$//;
    $def =~ s/^ +//;
    
    my $columnName = "";
    # get column name
    $def =~ /^([^ ]+)/ and $columnName = $1;

    #convert typedefinitions
    $def =~ s/mediumtext/long/gi;
    $def =~ s/int\([0-9]*\)/number/gi;
    $def =~ s/integer\(.*\)/number/gi;
    $def =~ s/INT /number /gi;
    $def =~ s/double\(.*\)/number/gi;
    $def =~ s/tinynumber/number/gi;
    $def =~ s/bignumber/number/gi;
    $def =~ s/varchar/varchar2/gi;
    $def =~ s/datetime/date/gi;
    $def =~ s/ time /varchar2(40)/gg;
    $def =~ s/auto_increment//gi;          # this can lead to bugs
    $def =~ s/unsigned//gi;
    $def =~ s/^start /start_point /gi;
    $def =~ s/^end /end_point /gi;
    #this is an interesting one; replace enum with check (value in
    $def =~ s/enum *(\(.*\))(.*)/varchar2(20) $2 CHECK \($columnName IN $1\)/i;
    # there is a nasty problem with the time things; this is a way to solve it
    # but is probably not correct....
    $def =~ s/\'0000-00-00 00\:00\:00\'/to_date\( \'01-01-0001\', 'dd-MM-yyyy' \)/;

    $DEBUG and print "--->DEBUG: $def\n";
    push @{$output{$table}->{'columns'}}, $def;
}
    
