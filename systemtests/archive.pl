#!usr/local/bin/perl

use strict;

use Bio::EnsEMBL::DBLoader;
use Getopt::Long;
use Bio::Seq;

my $dbtype = 'rdb';
my $host   = 'localhost';
my $port   = '410000';
my $dbname = 'archive';
my $dbuser = 'root';
my $dbpass = '';
my $module = 'Bio::EnsEMBL::DBArchive::Obj';
my $help;

&GetOptions( 
	     'dbtype:s'   => \$dbtype,
	     'host:s'     => \$host,
	     'port:n'     => \$port,
	     'dbname:s'   => \$dbname,
	     'dbuser:s'   => \$dbuser,
	     'dbpass:s'   => \$dbpass,
	     'module:s'   => \$module,
	     'h|help'     => \$help
	     );

my $locator = "$module/host=$host;port=$port;dbname=$dbname;user=$dbuser;pass=$dbpass";
my $db =  Bio::EnsEMBL::DBLoader->new($locator);

print "\nTesting get_seq... ";
my $seq = $db->get_seq ('ENSP0000dummy','1');
print "\n         id: ",$seq->id,"\ndescription: ",$seq->desc,"\n       type: ",$seq->type,"\n   sequence: ",$seq->seq,"\n\n";

print "Now testing get_seq_by_id...\n";
foreach $seq ($db->get_seq_by_id('ENSP0000dummy')) {
    print  "\n         id: ",$seq->id,"\ndescription: ",$seq->desc,"\n       type: ",$seq->type,"\n   sequence: ",$seq->seq,"\n\n";
}

print "Now testing get_seq_by_clone...\n";
foreach $seq ($db->get_seq_by_clone('dummy_clone_id','transcript')) {
    print  "\n         id: ",$seq->id,"\ndescription: ",$seq->desc,"\n       type: ",$seq->type,"\n   sequence: ",$seq->seq,"\n\n";
}

print "Now testing get_seq_by_clone_version...\n";
foreach $seq ($db->get_seq_by_clone_version('dummy_clone_id','1','transcript')) {
    print  "\n         id: ",$seq->id,"\ndescription: ",$seq->desc,"\n       type: ",$seq->type,"\n   sequence: ",$seq->seq,"\n\n";
}

print "Now testing get_seq_by_gene...\n";
foreach $seq ($db->get_seq_by_gene('dummy_gene_id2')) {
    print  "\n         id: ",$seq->id,"\ndescription: ",$seq->desc,"\n       type: ",$seq->type,"\n   sequence: ",$seq->seq,"\n\n";
}

print "Now testing get_seq_by_gene_version...\n";
foreach $seq ($db->get_seq_by_gene_version('dummy_gene_id2','1')) {
    print  "\n         id: ",$seq->id,"\ndescription: ",$seq->desc,"\n       type: ",$seq->type,"\n   sequence: ",$seq->seq,"\n\n";
}

print "Now testing write_seq... ";
my $newseq = Bio::Seq->new(
			   -seq=>'atcgcga',
			   -id=>'from_archive.pl',
			   -desc=>'Sequence from the EnsEMBL Archive database test system',
			   );
$db->write_seq ($newseq,'1','exon','gene_from_test','1','clone_from_test','1');
print "written!\n";
print "Testing delete_seq... ";
$db->delete_seq ('from_archive.pl','1');
print "deleted!\n";

 


