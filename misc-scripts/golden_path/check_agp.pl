#!/usr/local/bin/perl
use strict;
use Bio::EnsEMBL::DBSQL::Obj;
if ((scalar (@ARGV)) ==0) {
    print STDERR "No files given!\n";
    print STDERR "Usage: check_agp.pl coords all.agp nt.length\n";
    exit;
}
my $coords = shift (@ARGV);
my $agp = shift (@ARGV);
my $ntfile = shift (@ARGV);
my $cfile = shift (@ARGV);

my $db=Bio::EnsEMBL::DBSQL::Obj->new(-dbname=>'arne_oct07_tim',-host=>'ecs1b',-user=>'lijnzaad');
open (COORDS,"<$coords");
open (AGP,"<$agp");
open (NTFILE,"<$ntfile");
open (CLL,"<$cfile");

my %coords;
my %agp;
my %ntl;
my %cll;

while (<CLL>) {
    chomp;
    my (@tabs) = split(/\t/);
    $cll{$tabs[0]}=$tabs[1];
}

while (<COORDS>) {
    chomp;
    my (@tabs) = split(/\t/);
    $tabs[6] =~ /(\S+)\./;
    my $cl=$1;
    my $key = $tabs[0]."-".$cl;
    push (@{$coords{$key}},@tabs);
}
my %checknt;
my %seen;
my %seen2;
while (<AGP>) {
    chomp;
    my (@tabs) = split(/\t/);
    $tabs[6] =~ /(\S+)\./;
    my $cl=$1;
    my $key=$tabs[0]."-".$cl;
    my $frag_length=$tabs[3]-$tabs[2];
    $checknt{$tabs[0]}+= ($frag_length+1);
    if ($seen{$key}) {
	print STDERR "Duplicate line $key in agp file!\n";
    }
    $seen{$key}=1;
    if ($seen2{$tabs[6]}) {
	print STDERR "Clone $tabs[6] trying to be in two lines!\n";
    }
    $seen2{$tabs[6]}=1;

    if ($tabs[4] !~ /\+|\-/) {
	print STDERR "Expecting + or - in field 5, not $tabs[4] in line $_\n";
    }

    if ($tabs[5] !~ /F/) {
	print STDERR "Expecting F in field 6, not $tabs[5] in line $_\n";
    } 
    push (@{$agp{$key}},@tabs);
}

while (<NTFILE>) {
    chomp;
    my (@tabs) = split(/\t/);
    $ntl{$tabs[0]}=$tabs[1];
    my $check=$checknt{$tabs[0]};
    if ($check != $tabs[1]) {
	print STDERR "Sum of fragments for ".$tabs[0]." is $check while length should be ".$tabs[1]."\n";
    }
}

foreach my $nt_cl (keys(%agp)) {

    $nt_cl =~ /(\S+)\-(\S+)/;
    my $nt=$1;
    my $cl=$2;
    my $cl_length=$cll{$cl};
    
    my (@agp)=@{$agp{$nt_cl}};
    my $cl_start=$agp[7];
    my $cl_end=$agp[8];
    my $nt_start=$agp[2];
    my $nt_end=$agp[3];
    my $agp_strand=$agp[4];
    my $nt_length=$ntl{$nt};

    my (@coords)=@{$coords{$nt_cl}};
    my $coords_strand=$coords[4];

    #Checks!
    if ($nt_start > $nt_length) {
	print STDERR "NT start $nt_start for $nt_cl line is greater than nt length $nt_length!\n";
    }
    if ($nt_end > $nt_length) {
	print STDERR "NT end $nt_end for $nt_cl line is greater than nt length $nt_length!\n";
    }
    if ($cl_start > $cl_length) {
	print STDERR "CL start $cl_start for $nt_cl line is greater than cl length $cl_length!\n";
    }
    if ($cl_end > $cl_length) {
	print STDERR "CL end $cl_end for $nt_cl line is greater than cl length $cl_length!\n";
    }
    if ($agp_strand != $coords_strand) {
	print STDERR "Strand for line $nt_cl is different in coords and in agp\n";
    }
    

}


