# Load data from a file of regulatory regions into a database

use strict;

use DBI;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Analysis;
use Getopt::Long;

my ($host, $user, $pass, $port, $dbname, $file, $del, $motif);

GetOptions( "host=s",   \$host,
	    "user=s",   \$user,
	    "pass=s",   \$pass,
	    "port=i",   \$port,
	    "dbname=s", \$dbname,
	    "file=s",   \$file,
	    "motif=i",  \$motif,
	    "delete",   \$del,
	    "help",     \&usage
	  );


if (!$motif) {
  print "No motif ID specified, exiting\n";
  exit(1);
}

my $dbi = DBI->connect("dbi:mysql:host=$host;port=$port;database=$dbname",
		       $user,
		       $pass,
		       {'RaiseError' => 1}) || die "Can't connect to database";

delete_existing($dbi) if ($del);

my $rr_sth = $dbi->prepare("INSERT INTO regulatory_feature (name, seq_region_id, seq_region_start, seq_region_end, seq_region_strand, analysis_id, regulatory_motif_id, influence) VALUES(?,?,?,?,?,?,?,?)");
my $rr_obj_sth = $dbi->prepare("INSERT INTO regulatory_feature_object VALUES(?,?,?)");

my $db_adaptor = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $host,
						    -user => $user,
						    -port => $port,
						    -pass => $pass,
						    -dbname => $dbname);

my $analysis_adaptor    = $db_adaptor->get_AnalysisAdaptor;
my $slice_adaptor       = $db_adaptor->get_SliceAdaptor();
my $transcript_adaptor  = $db_adaptor->get_TranscriptAdaptor();
my $translation_adaptor = $db_adaptor->get_TranslationAdaptor();
my $gene_adaptor        = $db_adaptor->get_GeneAdaptor();

my $count = 0;

open (FILE, "<$file") || die "Can't open $file";

while (<FILE>) {

  next if ($_ =~ /^\s*\#/ || $_ =~ /^\s*$/);

  my ($group, $seq, $method, $feature, $chr, $start, $end, $str, $phase, $score, $type, $id_ignore, $id) = split;

  my $strand = ($str =~ /\+/ ? 1 : -1);

  # get analysis ID for $method
  my $analysis = $analysis_adaptor->fetch_by_logic_name($method);
  if (!$analysis) {
    print STDERR "Can't get analysis for $method, skipping\n";
    next;
  }

  # get seq_region id for $chr
  my $chr_slice = $slice_adaptor->fetch_by_region('chromosome', $chr, $start, $end, $strand);

  # get internal ID from $type and $id
  my ($ensembl_type) = $type =~ /(gene|transcript|translation)/;
  if (!$ensembl_type) {
    print STDERR "Can't get ensembl type from $type, skipping\n";
    next;
  }
  $ensembl_type = ucfirst(lc($ensembl_type));
  $id =~ s/[\"\']//g;  # strip quotes

  my $ensembl_object;

  if ($ensembl_type eq "Gene") {

    $ensembl_object = $gene_adaptor->fetch_by_stable_id($id);

  } elsif ($ensembl_type eq "Transcript") {

    $ensembl_object = $transcript_adaptor->fetch_by_stable_id($id);

  } elsif ($ensembl_type eq "Translation") {

    $ensembl_object = $translation_adaptor->fetch_by_stable_id($id);

  } else {

    print STDERR "Don't know how to deal with type $ensembl_type, skipping\n";
    next;

  }

  if (!$ensembl_object) {
    print STDERR "Can't find ensembl $ensembl_type corresponding to $id\n";
    next;
  }

  # populate tables
  # this should be done via the store methods of the API, once we have some ...
  $rr_sth->execute($seq,
		   $slice_adaptor->get_seq_region_id($chr_slice),
		   $chr_slice->start(),
		   $chr_slice->end(),
		   $chr_slice->strand(),
		   $analysis->dbID(),
		   $motif,
		   'negative');

  my $rr_id = $rr_sth->{'mysql_insertid'};

  $rr_obj_sth->execute($rr_id,
		      $ensembl_type,
		      $ensembl_object->dbID());

  $count++;
  # TODO - motifs, peptides

}

close FILE;

print "Uploaded $count regulatory regions to database\n";



sub delete_existing {

  my $dbi = shift;

  $dbi->do("DELETE FROM regulatory_feature");
  #$dbi->do("DELETE FROM regulatory_motif");
  $dbi->do("DELETE FROM regulatory_feature_object");
  $dbi->do("DELETE FROM peptide_regulatory_feature");

  print "Deleted existing data in regulatory region tables\n";

}


sub usage {

  print <<EOF;
Usage: perl load_regulatory.pl
         -host   : database host
         -port   : database port
         -user   : user name
         -pass   : password
         -dbname : database name
         -motif  : ID of motif to link to. Must already exist in database.
	 -file   : file to load from
	 -delete : delete existing contents of tables first

File format expected:

GROUP SEQ METHOD FEATURE CHR START END STRAND PHASE SCORE TYPE ID

e.g.

Similarity      hsa-miR-108     miRanda miRNA_target    3	49546906        49546927        +       .       107     transcript id   "ENST00000308775"
Similarity      hsa-miR-153     miRanda miRNA_target    3	49546683        49546702        +       .       106     transcript id   "ENST00000308775"

Lines beginning with # are ignored.

EOF

  exit(0);

}
