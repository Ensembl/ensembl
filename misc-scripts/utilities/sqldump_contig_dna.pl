#!/usr/local/bin/perl
use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my $db=Bio::EnsEMBL::DBSQL::DBAdaptor->new(-dbname=>'mouse',-host=>'ecs1a',-user=>'ensro');

my $sth = $db->prepare("select id from contig");
$sth->execute;
my (@cids) = $sth->fetchrow_array;
foreach my $contig (@cids) {
    $sth = $db->prepare("select concat('>',c.id,'\tlength: ',c.length,'\n',d.sequence) from dna d, contig c, static_golden_path sgp where sgp.raw_id = c.internal_id and c.dna = d.id and c.internal_id = $contig");
    $sth->execute;
    my ($string) = $sth->fetchrow_array();
    my $width = 80; 
    $string =~ s/(.{1,60})/$1\n/g;
    print "$string\n";
}

    
	
