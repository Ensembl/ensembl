## Bioperl Test Harness Script for Modules
##
# CVS Version
# $Id$


# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

#-----------------------------------------------------------------------
## perl test harness expects the following output syntax only!
## 1..3
## ok 1  [not ok 1 (if test fails)]
## 2..3
## ok 2  [not ok 2 (if test fails)]
## 3..3
## ok 3  [not ok 3 (if test fails)]
##
## etc. etc. etc. (continue on for each tested function in the .t file)
#-----------------------------------------------------------------------


## We start with some black magic to print on failure.
BEGIN { $| = 1; print "1..9\n"; 
	use vars qw($loaded); }
END {print "not ok 1\n" unless $loaded;}

use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::DB::WriteableVirtualContig;
$loaded=1;
print "ok \n";    # 1st test passed, loaded needed modules

#Creating test overlap database

$conf{'overlap'}      = 'testoverlap';
$conf{'mysqladmin'} = '/mysql/current/bin/mysqladmin';
$conf{'mysql'}      = '/mysql/current/bin/mysql';
$conf{'user'}       = 'root';
$conf{'perl'}       = 'perl';

if ( -e 't/overlap.conf' ) {
  print STDERR "Reading configuration from overlap.conf\n";
  open(C,"t/overlap.conf");
  while(<C>) {
    my ($key,$value) = split;
    $conf{$key} = $value;
  }
} else {
  print STDERR "Using default values\n";
  foreach $key ( keys %conf ) {
    print STDERR " $key $conf{$key}\n";
  }
  print STDERR "\nPlease use a file t/overlap.conf to alter these values if the test fails\nFile is written <key> <value> syntax\n\n";
}

$nuser = $conf{user};

my $create_overlap        = "$conf{mysqladmin} -u ".$nuser." create $conf{overlap}";

system($create_overlap)   == 0 or die "$0\nError running '$create_overlap' : $!";

print "ok 2\n";    #Databases created successfuly

#Initialising databases
my $init_overlap        = "$conf{mysql} -u ".$nuser." $conf{overlap} < ../sql/table.sql";

system($init_overlap)     == 0 or die "$0\nError running '$init_overlap' : $!";

print "ok 3\n";

#Suck test data into db
print STDERR "Inserting test data in test overlap db...\n";
my $suck_data      = "$conf{mysql} -u ".$nuser." $conf{overlap} < t/overlap.dump";
system($suck_data) == 0 or die "$0\nError running '$suck_data' : $!";

print "ok 4\n";

# Connect to test db

my $db             = new Bio::EnsEMBL::DBSQL::Obj(-host   => 'localhost',
						  -user   => $conf{user},
						  -dbname => $conf{overlap}
						 );

die "$0\nError connecting to database : $!" unless defined($db);

print "ok 5\n";

my $contig = $db->get_Contig('contig1');

die "$0\nError fetching contig1 : $!" unless defined ($contig);

print "ok 6\n";

my $wvc     = new Bio::EnsEMBL::DB::WriteableVirtualContig(-focuscontig   => $contig,
						 -focusposition => 1,
						 -ori           => 1,
						 -left          => 20,
						 -right         => 20,
                                                 -gene_obj      => $db->gene_Obj);

die ("$0\nCan't create virtual contig :$!") unless defined ($wvc);

$wvc->_dump_map(\*STDERR);

print "ok 7\n";

$gene = Bio::EnsEMBL::Gene->new();
$gene->id('gene-id-1');
$gene->version(1);
$trans = Bio::EnsEMBL::Transcript->new();
$trans->id('trans-id-1');
$trans->version(1);
$trl = Bio::EnsEMBL::Translation->new();
$trl->start_exon_id('exon-1');
$trl->end_exon_id('exon-2');
$trl->start(2);
$trl->end(13);
$trl->id('trl-id');
$trl->version(1);
$trans->translation($trl);
$trans->created(1);
$trans->modified(1);
$gene->add_Transcript($trans);
$gene->created(1);
$gene->modified(1);

$exon = Bio::EnsEMBL::Exon->new();
$exon->id('exon-1');
$exon->start(1);
$exon->end(5);
$exon->strand(1);
$exon->version(1);
$exon->phase(0);
$exon->created(1);
$exon->modified(1);
$exon->contig_id($wvc->id);
$trans->add_Exon($exon);
$exon{'exon-1'} = $exon;
$exon->attach_seq($wvc->primary_seq);

$exon = Bio::EnsEMBL::Exon->new();
$exon->id('exon-2');
$exon->start(8);
$exon->end(15);
$exon->strand(1);
$exon->version(1);
$exon->phase(0);
$exon->created(1);
$exon->modified(1);
$exon->contig_id($wvc->id);
$trans->add_Exon($exon);
$exon{'exon-2'} = $exon;
$exon->attach_seq($wvc->primary_seq);

$wvc->write_Gene($gene);

print "ok 8\n";

$newgene = $db->get_Gene('gene-id-1');

$error = 0;
foreach $exon ( $newgene->each_unique_Exon ) {
    if( exists $exon{$exon->id} ) {
	if( $exon->seq != $exon{$exon->id}->seq ) {
	    print STDERR "For exon ",$exon->id," sequences do not agree!\n";
	    print STDERR "Exon in RC is $exon [",$exon->seq->seq,"] VC is [",$exon{$exon->id}->seq->seq,"]\n";
	    $error =1;
	} else {
	    print STDERR "Got exon match for",$exon->id,"\n";
	}
    } else {
	print STDERR "Weird exon in output...\n";
    }
}

if( $error == 1 ) {
    print "not ok 9\n";
} else {
    print "ok 9\n";
}

END {
   my $drop_overlap        = "echo \"y\" | $conf{mysqladmin} -u ".$nuser." drop $conf{overlap}";
  system($drop_overlap)     == 0 or die "$0\nError running '$drop_overlap' : $!";
}




