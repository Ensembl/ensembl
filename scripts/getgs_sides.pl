#!/usr/local/bin/perl


use Getopt::Long;
use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DBLoader;


# global defaults
my $host = 'localhost';
my $dbuser = 'root';
my $dbname = 'ensembl_freeze24';
my $dbpass = undef;
my $module = 'Bio::EnsEMBL::DBSQL::Obj';
my $port = '410000';
my $usefile = 0;

&GetOptions(
             'dbuser:s'  => \$dbuser,
             'dbpass:s'  => \$dbpass,
             'host:s'    => \$host,
             'dbname:s'  => \$dbname,
             'port:n'    => \$port,
             'module:s'  => \$module,
             'usefile'   => \$usefile
             );
my @est;

my $locator = "$module/host=$host;port=$port;dbname=$dbname;user=$dbuser;pass=$dbpass;perlonlyfeatures=1";
print STDERR "Using $locator for todb\n";
my $db =  Bio::EnsEMBL::DBLoader->new($locator);


$raw_id = shift;

$contig = $db->get_Contig($raw_id);

$vc = Bio::EnsEMBL::DB::VirtualContig->new( -focuscontig => $contig,
					    -focusposition => $contig->golden_start + 10,
					    -ori => 1,
					    -left => 20000000,
					    -right => 20000000,
					    );


#$usec = $contig;
$usec  = $vc;

$vc->_dump_map(\*STDOUT);

foreach $gp ( $usec->get_all_PredictionFeatures ) {
    foreach $sub ( $gp->sub_SeqFeature ) {
	if( $sub->strand == 1 ) {
	    $sub->attach_seq($usec->primary_seq);

	    print "For ",$sub->entire_seq->subseq($sub->start-2,$sub->start-1),
	    ":",$sub->entire_seq->subseq($sub->end+1,$sub->end+2),"\n";
	} else {
	    $sub->attach_seq($usec->primary_seq);

	    print "Rev ",$sub->entire_seq->trunc($sub->end+1,$sub->end+2)->revcom->seq,
	    ":",$sub->entire_seq->trunc($sub->start-2,$sub->start-1)->revcom->seq,"\n";
	}
    }
}

#$vc->dump_agp('chr1',\*STDOUT);
