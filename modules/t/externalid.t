
## We start with some black magic to print on failure.
BEGIN { $| = 1; print "1..6\n"; 
	use vars qw($loaded); }
END {print "not ok 1\n" unless $loaded;}


use lib 't';
use EnsTestDB;
$loaded = 1;
print "ok 1\n";    # 1st test passes.
    
my $ens_test = EnsTestDB->new();

$gene_obj = $ens_test->get_DBSQL_Obj->gene_Obj();

@id = $gene_obj->get_New_external_id('gene','ENSG',2);


if( $id[0] ne 'ENSG00000000001' || $id[1] ne 'ENSG00000000002' ) {
	print "not ok 2\n";
} else {
	print "ok 2\n";
}


@id = $gene_obj->get_New_external_id('gene','ENSG',2);


if( $id[0] ne 'ENSG00000000003' || $id[1] ne 'ENSG00000000004' ) {
	print "not ok 3\n";
} else {
	print "ok 3\n";
}

@id = $gene_obj->get_New_external_id('transcript','ENST',2);


if( $id[0] ne 'ENST00000000001' || $id[1] ne 'ENST00000000002' ) {
	print "not ok 4\n";
} else {
	print "ok 4\n";
}

@id = $gene_obj->get_New_external_id('exon','ENSE',2);


if( $id[0] ne 'ENSE00000000001' || $id[1] ne 'ENSE00000000002' ) {
	print "not ok 5\n";
} else {
	print "ok 5\n";
}


@id = $gene_obj->get_New_external_id('translation','ENSP',5);


if( $id[0] ne 'ENSP00000000001' || $id[4] ne 'ENSP00000000005' ) {
	print "not ok 6\n";
} else {
	print "ok 6\n";
}




