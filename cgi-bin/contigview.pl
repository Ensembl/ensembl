#!/usr/local/bin/perl 

# makes GFF stuff for a contig.

BEGIN {
    push(@INC,"../modules");
    push(@INC,"../../bioperl-live");
}

use CGI;
use Bio::EnsEMBL::DBSQL::Obj;
use strict;

my $q = new CGI;
print $q->header();
#print "content-type: text\n\n";

my $contigid = $q->param('contig');
my @features;


eval {
    my $db = new Bio::EnsEMBL::DBSQL::Obj( -user => 'root', -db => 'ensdev' , -host => 'croc.sanger.ac.uk');
    my $contig = $db->get_Contig($contigid);
};

if( $@ ) {
    print "<p>Warning! Exception<p>\n<pre>\n$@\n</pre>\n";
} else {
    print "<h2>Contig for $contigid</h2><p><i>No graphic yet</i>\n\n<p>[<a href=\"/cgi-test/ensembl/cgi-bin/contig2gff.pl?contig=$contigid\">GFF</a>]</p>\n";
}


    
