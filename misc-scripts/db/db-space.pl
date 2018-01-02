#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


use strict;
use warnings;


use Getopt::Long;
my @options;

# Get output immediately. It won't hurt performance.
use FileHandle;
autoflush STDERR;
autoflush STDOUT;

my $debug = 0;

my $pw;
push(@options, "password=s", \$pw);

my $host = "localhost";
push(@options, "host=s", \$host);

my $user = "localuser";
push(@options, "user=s", \$user);

my $port = "3306";
push(@options, "port=s", \$port);

die "Couldn't parse options" if !GetOptions(@options);

my $cmd = mysql_cmd("show databases");
open(CMD, $cmd) or die "Couldn't $cmd: $!\n";
my @databases;
my $header = <CMD>;
while ( <CMD> ) {
  s/[\r\n]$//g;
  #print "$_\n";
  push (@databases, $_);
}
close(CMD);

#print "@databases";

my %colmap = ( 'Data_length' => 6,
               'Index_length' => 8,
               'Engine' => 1,
               'Comment' => 17 );

my %size;
my %total_size;
my %engine_map;
my %top_tables;

my $inno_db_free;
foreach my $db (@databases) {
  print STDERR ".";
  $cmd = mysql_cmd("use $db; show table status");

  open(CMD, $cmd) or die "Couldn't $cmd: $!\n";
  my $header = <CMD>;
  my $total_size = 0;
  if (defined($header)) {
    $header =~ s/[\r\n]$//g;
    my @head = split("\t", $header);

    foreach my $col (keys %colmap) {
      die "$db: Expected '$col', found '" . $head[$colmap{$col}] . "'"
        if $head[$colmap{$col}] ne $col;
    }

    while (<CMD>) {
      my @data = split("\t");
      print STDERR "== $db\n" if ($debug);
      my ($data_length, $index_length) = @data[6,8];
      my ($engine, $comment) = @data[1,17];
      my $tbl_name = $data[0];
      $engine_map{$engine}++;
      $index_length = 0 if (!defined($index_length) || $index_length eq 'NULL');
      $data_length = 0 if (!defined($data_length) || $data_length eq 'NULL');
      $engine = 'Unknown' if (!defined($engine) || $engine eq 'NULL');
      $db = 'Unknown' if (!defined($db) || $db eq 'NULL');
      $size{$db}{$engine} += $data_length + $index_length;
      $total_size{$db} += $data_length + $index_length;
      my $length = $data_length + $index_length;
      $top_tables{$length}{$tbl_name}{$db} = 1;

      if ( $comment =~ /InnoDB free: (\d+) kB/ ) {
        warn "Found two different inno DB free sized - $db - $tbl_name.\n"
          if defined($inno_db_free) && $inno_db_free != $1;
        $inno_db_free = $1;
      }
    }
    close(CMD);
  }
}
print STDERR "\n";

print "NOTE: All numbers are in megabytes (M).\n";
printf("Inno DB free: %.1f\n", $inno_db_free / 1024)
  if defined($inno_db_free);

printf("%-40s ", "database");
foreach my $engine (sort keys(%engine_map)) {
  printf "%7s ", $engine;
}
printf "%8s", "total";
print "\n";

foreach my $db (sort {$total_size{$b} <=> $total_size{$a}} keys %total_size) {
  printf("%-40s ", $db);
  foreach my $engine (sort keys(%engine_map)) {
    my $size= $size{$db}{$engine};
    $size = 0 if !defined($size);
    printf("%7.1f ", $size / 1024 / 1024);
  }
  printf("%8.1f\n", $total_size{$db} / 1024 / 1024);
}

my $total_bytes;
map {$total_bytes+=$_} values %total_size;
print "================================\n";
printf("TOTAL SPACE USED %7.1f\n", $total_bytes / 1024 / 1024);
print "================================\n";

my $count++;
print "==================\n";
print "Top tables by size\n";
print "==================\n";
foreach my $size (sort {$b<=>$a} keys %top_tables) {
  last if ($count > 5);
  my @tbl_names = keys %{$top_tables{$size}};
  my @dbs = keys %{$top_tables{$size}{$tbl_names[0]}};
  printf("%-40s ", $dbs[0]);
  printf("%-40s ", $tbl_names[0]);
  printf("%7.1f ", $size / 1024 / 1024);
  print("\n");
  $count++;
}

sub mysql_cmd {
  my $mysql_cmd = shift;

  my $pw_args = $pw ? "-p$pw" : '';
  return "mysql -uroot -h$host -u$user $pw_args -P$port -e '$mysql_cmd'|";
}

