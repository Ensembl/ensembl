#!/usr/local/bin/perl 

# makes GFF stuff for a contig.

BEGIN {
    push(@INC,"../modules");
    push(@INC,"../../bioperl-live");
}

use CGI;
use Bio::EnsEMBL::DB::Obj;
use strict;

my $q = new CGI;
print $q->header();
#print "content-type: text\n\n";

my $contigid = $q->param('contig');
my @features;

eval {
    my $db = new Bio::EnsEMBL::DB::Obj( -user => 'root', -db => 'pog' , -host => 'caldy.sanger.ac.uk');
    my $contig = $db->get_Contig($contigid);
    @features = $contig->get_all_SeqFeatures;
};

if( $@ ) {
    print "<p>Warning! Exception<p>\n<pre>\n$@\n</pre>\n";
} else {

    print "<pre>\n";
    foreach my $sf ( @features ) {
	# $sf is Bio::SeqFeature::Generic object.
	print $sf->gff_string, "\n";
    }
    print "</pre>\n";
}


    
