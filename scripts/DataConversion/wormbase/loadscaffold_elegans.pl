#!/usr/local/bin/perl -w

use strict;
use vars qw($USER $PASS $DB $HOST $DSN);
use Bio::EnsEMBL::DBLoader;
use Bio::SeqIO;
use Bio::EnsEMBL::RawContig;
use Bio::EnsEMBL::Clone;
use Getopt::Long;

my $dbtype = 'rdb';
my $host   = 'ecs1d';
my $port   = '';
my $dbname = 'elegans_newschema';
my $dbuser = 'ecs1dadmin';
my $dbpass = 'TyhRv';
my $module = 'Bio::EnsEMBL::DBSQL::DBAdaptor';
$| = 1;
my $help = 0;
my $pipe = 1;
my $verbose = 0;

&GetOptions( 
	     'dbtype:s'   => \$dbtype,
	     'host:s'     => \$host,
	     'port:n'     => \$port,
	     'dbname:s'   => \$dbname,
	     'dbuser:s'   => \$dbuser,
	     'dbpass:s'   => \$dbpass,
	     'module:s'   => \$module,
             'pipe'       => \$pipe,
	     'verbose'    => \$verbose,
	     'h|help'     => \$help
	     );
if ($help) {
    exec('perldoc', $0);
}


$SIG{INT} = sub {my $sig=shift;die "exited after SIG$sig";};
my $locator;
if ($dbpass) {
   $locator = "$module/host=$host;port=;dbname=$dbname;user=$dbuser;pass=$dbpass";
} else {
   $locator = "$module/host=$host;port=;dbname=$dbname;user=$dbuser;pass=";
}
my $db =  Bio::EnsEMBL::DBLoader->new($locator);

my ($seqfile) = shift;
my ($contig_id_file) = shift;

open(FH, $contig_id_file);

my %seq_got;

if( !defined $seqfile ) { die 'cannot load because sequence file to load sequences from was not passed in as argument';}

my $seqio = new Bio::SeqIO(-format => 'Fasta', 
			   -file => $seqfile);
my $count = 0;

while( my $seq = $seqio->next_seq ) {
  my @clone_ids = split / /, $seq->desc;
  #print $seq->id."\n";
  #print $seq->desc."\n";
  if($seq->seq =~ /match/i){
    die("seq = ".$seq->seq."\n");
  }
  if($seq->seq =~ /abort/i){
    die("seq = ".$seq->seq."\n");
  }  
  my $clone     = new Bio::EnsEMBL::Clone;
  my $contig    = new Bio::EnsEMBL::RawContig;
  my ($version) = $clone_ids[0] =~ /\S+\.(\d+)/;
  $seq_got{$clone_ids[0]} = 1;
  $clone->htg_phase(4);
  $clone->id($clone_ids[0]); 
  $clone->embl_id($clone_ids[0]);
  $clone->version(1);
  $clone->embl_version($version);
  $clone->adaptor($db->get_CloneAdaptor);
  my $contig_id = $clone_ids[0].".1.".$seq->length;
  $contig->name($contig_id);
  $contig->dbID($count++);
  $contig->seq($seq->seq);
  $contig->length($seq->length);
  $contig->adaptor($db->get_RawContigAdaptor);
  $clone->add_Contig($contig);
  my $time = time;
  $clone->created($time);
  $clone->modified($time);
  #print "have clone with ".$clone->id." contig ".$contig->id."\n";
  eval { 
    $db->get_CloneAdaptor->store($clone);
    $verbose && print STDERR "Written ".$clone->id." scaffold into db\n";
  };
  if( $@ ) {
    print STDERR "Could not write clone into database, error was $@\n";
  }
}


while(<FH>){
  chomp;
  my $name = $_;
  if(!$seq_got{$name}){
    print "don't seem to have sequence for ".$name."\n";
  }
}

close(FH);
