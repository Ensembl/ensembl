#!/usr/local/bin/perl
use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my $db=Bio::EnsEMBL::DBSQL::DBAdaptor->new(-dbname=>'mouse',-host=>'ecs1a',-user=>'ensro');

my $sth = $db->prepare("select internal_id from contig c, static_golden_path sgp where raw_id = c.internal_id");
$sth->execute;
my @cids; 
while (my ($contig) = $sth->fetchrow_array ) {
    push (@cids,$contig);
}
foreach my $contig (@cids) {
    $sth = $db->prepare("select concat('>',c.id,'\tlength: ',c.length,'\n',d.sequence) from dna d, contig c where c.dna = d.id and c.internal_id = $contig");
    #print STDERR "select concat('>',c.id,'\tlength: ',c.length,'\n',d.sequence) from dna d, contig c where c.dna = d.id and c.internal_id = $contig\n";
    $sth->execute;
    my ($string) = $sth->fetchrow_array();
    my $width = 60; 
    $string =~ s/(.{1,$width})/$1\n/g;
    $string =~ s/\n//;
    print "$string";
}
$db->DESTROY;
    
	
