#!/usr/local/bin/perl

use Bio::EnsEMBL::DBLoader;

$dblocator = "Bio::EnsEMBL::DBSQL::Obj/host=ensrv4.sanger.ac.uk;user=ensro;dbname=ensembl_freeze17_michele";
print STDERR "About to call new..\n";
my $db = Bio::EnsEMBL::DBLoader->new($dblocator);
print STDERR "got db object\n";

$db->static_golden_path_type('UCSC');

@vclist = $db->get_StaticGoldenPathAdaptor->fetch_VirtualContig_list_sized('ctg12824',2000000,100000,4000000,100);


foreach $vc ( @vclist ) {
    print "length is ",$vc->length,"\n";
}

