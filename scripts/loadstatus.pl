#!/usr/local/bin/perl

=head1 NAME - loadstatus.pl

 compares loading status of RDB and TimDB

=head1 SYNOPSIS - 

    loadstatus.pl

=head1 OPTIONS

=cut

use strict;

# comment out ace as most people don't have aceperl.
# need a run-time loader really.
#use Bio::EnsEMBL::AceDB::Obj;
use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::TimDB::Obj;
use Bio::EnsEMBL::EMBL_Dump;
use Bio::AnnSeqIO;
use Bio::SeqIO;

use Getopt::Long;

# defaults for msql (rdb) access
# msql was 'croc'
my $host      = 'obi-wan';
# msql was 'ensdev'
my $dbname    = 'ens2';
my $dbuser    = 'humpub';
my $dbpass    = 'ens2pass';

my $help;
my $verbose;

&GetOptions( 'help|h'    => \$help,
	     'dbuser:s'  => \$dbuser,
	     'dbpass:s'  => \$dbpass,
	     'dbname:s'  => \$dbname,
	     'host:s'    => \$host,
	     'verbose'   => \$verbose,
	     );

if($help){
    exec('perldoc', $0);
}

my $db1 = Bio::EnsEMBL::DBSQL::Obj->new( -user => $dbuser, -db => $dbname , 
				      -host => $host, -password => $dbpass );
# just need a dummy clone here (but must exist!)
my $raclones=['AP000228'];
my $db2 = Bio::EnsEMBL::TimDB::Obj->new($raclones);

print STDERR "Making clone list from DBSQL\n";
my @clones1 = $db1->get_all_Clone_id();
print STDERR scalar(@clones1)." clones found\n\n";
my %clones1=map{ $_, 1 } @clones1;

print STDERR "Making clone list from TimDB\n";
my @clones = $db2->get_all_Clone_id();
print STDERR scalar(@clones)." clones found\n\n";
# for entries from TimDB, need to convert to accession numbers
my %clones2=map{ ($db2->get_id_acc($_))[0], $_ } @clones;
my @clones2 = (keys %clones2);

# check what clones are only in TimDB
my @clone2only;
foreach my $clone (@clones2){
    if(!$clones1{$clone}){
	push(@clone2only,$clone);
    }
}
my $nclone2only=scalar(@clone2only);
if($nclone2only){
    print STDERR "$nclone2only clones not loaded into DBSQL\n";
    foreach my $clone (@clone2only){
	print STDERR "$clone [$clones2{$clone}]\n" if $verbose;
    }
}
print STDERR "\n";

# check what clones are only in DBSQL
my @clone1only;
foreach my $clone (@clones1){
    if(!$clones2{$clone}){
	push(@clone1only,$clone);
    }
}
my $nclone1only=scalar(@clone1only);
if($nclone1only){
    print STDERR "$nclone1only clones only in DBSQL\n";
    foreach my $clone (@clone1only){
	print STDERR "$clone\n" if $verbose;
    }
}
print "\n";

