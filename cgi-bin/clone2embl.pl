#!/usr/local/bin/perl



BEGIN {
    unshift(@INC,"../modules");
    unshift(@INC,"../../bioperl-live");
}

use CGI;
use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::EMBL_Dump;
use Bio::AnnSeqIO::EMBL; # compile time checking of EMBL dumper
use Bio::AnnSeqIO;
use strict;

my $q = new CGI;
print $q->header();
#print "content-type: text\n\n";

my $clone = $q->param('clone');
my @contigs;


eval {
    my $db = new Bio::EnsEMBL::DBSQL::Obj( -user => 'root', -db => 'ensdev' , -host => 'croc.sanger.ac.uk');
    my $clone = $db->get_Clone($clone);

    my $as = $clone->get_AnnSeq();

    &Bio::EnsEMBL::EMBL_Dump::add_ensembl_comments($as);
    my $emblout = Bio::AnnSeqIO->new( -format => 'EMBL', -fh => \*STDOUT);
    &Bio::EnsEMBL::EMBL_Dump::ensembl_annseq_output($emblout);

    print "<pre>\n";
    $emblout->write_annseq($as);
    print "</pre>\n";


};

if( $@ ) {
    print "Apologies. An exception occured!<p><pre>$@</pre>\n";
}
