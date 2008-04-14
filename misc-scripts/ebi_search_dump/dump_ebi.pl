# Dump gene and xref information to an XML file for indexing by the EBI's
# search engine.
#
# To copy files to the EBI so that they can be picked up:
# scp *.xml.gz glenn@puffin.ebi.ac.uk:xml/
#
# Email eb-eye@ebi.ac.uk after copying so the files can be indexed.

use strict;
use DBI;

use Getopt::Long;
use IO::Zlib;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;

use HTML::Entities;

my ($host, $user, $pass, $port, $dbpattern, $max_genes, $nogzip, $no_variation, $parallel, $dir, $no_vega);

GetOptions( "host=s",              \$host,
	    "user=s",              \$user,
	    "pass=s",              \$pass,
	    "port=i",              \$port,
	    "dbpattern|pattern=s", \$dbpattern,
	    "nogzip!",             \$nogzip,
            "max_genes=i",         \$max_genes,
	    "no_variation",        \$no_variation,
	    "no_vega",             \$no_vega,
	    "parallel",            \$parallel,
	    "dir=s",               \$dir,
	    "help" ,               \&usage
	  );

$user    = $user || "ensro";
$host    = $host || "ens-staging";
$port    = $port ||  "3306";
$dbpattern = $dbpattern || "_core_";
$dir     = $dir || "/lustre/scratch1/ensembl/gp1/xml";

my $entry_count = 0;

my $fh;

run();

sub run() {

  Bio::EnsEMBL::Registry->load_registry_from_db(-host => $host,
						-port => $port,
						-user => $user,
						-pass => $pass);

  # loop over databases

  my $dsn = "DBI:mysql:host=$host";
  $dsn .= ";port=$port" if ($port);

  my $db = DBI->connect( $dsn, $user, $pass );

  my @dbnames = map {$_->[0] } @{$db->selectall_arrayref("show databases")};

  for my $dbname (@dbnames) {

    next if ($dbname !~ /$dbpattern/);

    my $file = $dir . "/" . $dbname . ".xml";
    $file .= ".gz" unless ($nogzip);

    if ($parallel) {

      submit($dbname, $file);

    } else {

      dump_single($dbname, $file);

    }

  }

}

# -------------------------------------------------------------------------------

sub submit {

  my ($dbname, $file) = @_;

  print "Submitting job for $dbname\n";

  my $o = $dir . "/" . ${dbname} . ".out";
  my $e = $dir . "/" . ${dbname} . ".err";

  my $p = ($pass) ? "-pass $pass" : '';

  my $n = substr($dbname, 0, 10);

  my $g = ($nogzip) ? "-nogzip" : '';

  system "bsub -o $o -e $e -J $n perl dump_ebi.pl -user $user -host $host $p -port $port -dbpattern $dbname $g";

}

# -------------------------------------------------------------------------------

sub dump_single {

  my ($dbname, $file) = @_;

  unless ($nogzip) {

    $fh = new IO::Zlib;
    $fh->open("$file", "wb9") || die ("Can't open compressed stream to $file");

  } else  {

    open(FILE, ">$file") || die "Can't open $file";

  }

  print "Dumping $dbname to $file\n";

  my $start_time = time;

  my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor('-host' => $host,
					       '-port' => $port,
					       '-user' => $user,
					       '-pass' => $pass,
					       '-dbname' => $dbname,
					       '-species' => $dbname);

  header($dba, $dbname);

  content($dba);

  footer();

  print_time($start_time);

}

# -------------------------------------------------------------------------------

sub header {

  my ($dba, $dbname) = @_;

  p ("<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>");
  p ("<!DOCTYPE database [ <!ENTITY auml \"&#228;\">]>");
  p ("<database>");
  p ("<name>$dbname</name>");

  my $meta_container = $dba->get_MetaContainer();
  my $species = @{$meta_container->list_value_by_key('species.common_name')}[0];

  p ("<description>Ensembl $species database</description>");

  my $release = @{$meta_container->list_value_by_key('schema_version')}[0];

  p ("<release>$release</release>");

  p ("");
  p ("<entries>");


}

# -------------------------------------------------------------------------------

sub content {

  my ($dba) = @_;

  dump_genes($dba);

  unless ($no_vega) {

    my $db_vega = vega_attach($dba);

    dump_genes($db_vega) if ($db_vega);

  }

}

# -------------------------------------------------------------------------------

sub dump_genes {

  my ($dba) = @_;

  my $gene_adaptor = $dba->get_GeneAdaptor();

  my $meta_container = $dba->get_MetaContainer();
  my $species = $meta_container->get_Species()->common_name();

  my $db_variation = variation_attach($dba) unless $no_variation;

  my $trv_adaptor;
  if ($db_variation) { # not all species have variation databases
    $trv_adaptor = $db_variation->get_TranscriptVariationAdaptor();
  }

  my ($db_adaptor) = @_;

  my $genes = $gene_adaptor->fetch_all();

  while (my $gene = shift(@$genes)) {

    last if ($max_genes && $entry_count >= $max_genes);

    # general gene info
    p("");
    p ("<entry id=\"" . $gene->stable_id() . "\">");

    p ("<name>" . $gene->display_id() . "</name>");

    my $description = encode_entities($gene->description()); # do any other fields need encoding?

    p ("<description>" . $description . "</description>");

    my $created = $gene->created_date();
    my $modified = $gene->modified_date();

    if ($created || $modified) {
      p ("<dates>");
      p ("<date type=\"creation\" value=\"" . format_date($created) . "\"/>") if ($created);
      p ("<date type=\"last_modification\" value=\"" . format_date($modified) . "\"/>") if ($modified);
      p ("</dates>");
    }

    # xrefs - deal with protein_ids separately, as additional fields
    my @protein_ids;

    p ("<cross_references>");

    my $tax = $meta_container->get_taxonomy_id();
    p ("<ref dbname=\"taxonomy\" dbkey=\"$tax\"/>");

    my $xrefs = $gene->get_all_DBLinks();

    while (my $xref = shift(@$xrefs)) {

      if ($xref->dbname() !~ /protein_id/) {
	my $display_id = $xref->display_id();
	$display_id =~ s/<.+>//g;
	p ("<ref dbname=\"" . $xref->dbname() ."\" dbkey=\"" . $display_id. "\"/>");
      } else {
	push @protein_ids, $xref->display_id();
      }
    }

    p ("</cross_references>");

    # additional fields - transcript, translation, species etc
    p ("<additional_fields>");

    my $transcripts = $gene->get_all_Transcripts();

    foreach my $transcript (@$transcripts) { # can't use while/shift here as array needed later

      p ("<field name=\"transcript\">" . $transcript->stable_id() . "</field>");

      my $translation = $transcript->translation();
      p ("<field name=\"translation\">" . $translation->stable_id() . "</field>") if ($translation);

    }

    while (my $protein_id = shift(@protein_ids)) {
      p ("<field name=\"protein_id\">" . $protein_id . "</field>");
    }

    p ("<field name=\"species\">" . $species . "</field>");

    # SNP IDs
    if ($db_variation) {
      my $transcript_variations = $trv_adaptor->fetch_all_by_Transcripts($transcripts);
      while (my $tv = shift(@$transcript_variations)) {
	p ("<field name=\"variation_id\">" . $tv->variation_feature()->variation_name() . "</field>");
      }
    }

    p ("</additional_fields>");

    # close tag
    p ("</entry>");

    $entry_count++;

  }

}

# -------------------------------------------------------------------------------

sub footer {


  p ("</entries>");

  p ("<entry_count>$entry_count</entry_count>");

  p ("</database>");

  print "Dumped $entry_count entries\n";

  if ($nogzip) {

    close(FILE);

  } else {

    $fh->close();

  }

}


# -------------------------------------------------------------------------------

sub p {

  my $str = shift;

  # TODO - encoding

  $str .= "\n";

  if ($nogzip) {

    print FILE $str;

  } else {

    print $fh $str;

  }

}

# -------------------------------------------------------------------------------

sub format_date {

  my $t = shift;

  my ($y, $m, $d, $ss, $mm, $hh) = (localtime($t))[5,4,3,0,1,2];
  $y += 1900;
  $d = "0" . $d if ($d < 10);
  my $mm = text_month($m);

  return "$d-$mm-$y";

}

# -------------------------------------------------------------------------------

sub text_month {

  my $m = shift;

  my @months = qw[JAN FEB MAR APR MAY JUN JUL AUG SEP OCT NOV DEC];

  return $months[$m];

}

# -------------------------------------------------------------------------------

sub print_time {

  my $start = shift;

  my $t = time - $start;
  my $s = $t % 60;
  $t = ($t - $s) / 60;
  my $m = $t % 60;
  $t = ($t - $m) /60;
  my $h = $t % 60;

  print "Time taken: " . $h . "h " . $m . "m " .$s . "s\n";

}

# -------------------------------------------------------------------------------

#
# Figure out the name of a variation database from the core database name
#

sub variation_attach {

  my $db = shift;

  my $core_db_name;
  $core_db_name = $db->dbc->dbname();
  return undef if ($core_db_name !~ /_core_/);

  my $dbc = $db->dbc();
  my $sth = $dbc->prepare("show databases");
  $sth->execute();
  my $all_db_names = $sth->fetchall_arrayref();
  my %all_db_names = map {( $_->[0] , 1)} @$all_db_names;
  my $variation_db_name = $core_db_name;
  $variation_db_name =~ s/_core_/_variation_/;

  return undef if (! exists $all_db_names{$variation_db_name});

  # register the dbadaptor with the Registry
  return Bio::EnsEMBL::Variation::DBSQL::DBAdaptor->new(-host => $dbc->host(),
							-user => $dbc->username(),
							-pass => $dbc->password(),
							-port => $dbc->port(),
							-dbname => $variation_db_name);

}

# -------------------------------------------------------------------------------

#
# Figure out the name of a vega database from the core database name
#

sub vega_attach {

  my $db = shift;

  my $core_db_name;
  $core_db_name = $db->dbc->dbname();
  return undef if ($core_db_name !~ /_core_/);

  my $dbc = $db->dbc();
  my $sth = $dbc->prepare("show databases");
  $sth->execute();
  my $all_db_names = $sth->fetchall_arrayref();
  my %all_db_names = map {( $_->[0] , 1)} @$all_db_names;
  my $vega_db_name = $core_db_name;
  $vega_db_name =~ s/_core_/_vega_/;

  return undef if (! exists $all_db_names{$vega_db_name});

  # register the dbadaptor with the Registry
  return Bio::EnsEMBL::DBSQL::DBAdaptor->new(-host => $dbc->host(),
					     -user => $dbc->username(),
					     -pass => $dbc->password(),
					     -port => $dbc->port(),
					     -dbname => $vega_db_name);

}

# -------------------------------------------------------------------------------

sub usage {
  print <<EOF; exit(0);

Usage: perl $0 <options>

  -host         Database host to connect to. Defaults to ens-staging.

  -port         Database port to connect to. Defaults to 3306.

  -dbpattern    Database name regexp. Defaults to _core_

  -user         Database username. Defaults to ensro.

  -pass         Password for user.

  -dir          Directory to write output to. Defaults to /lustre/scratch1/ensembl/gp1/xml.

  -nogzip       Don't compress output as it's written.

  -max_genes    Only dump this many genes for testing.

  -no_variation Don't dump variation IDs.

  -no_vega      Don't attach to Vega databases.

  -parallel     Submit jobs in parallel.

  -help         This message.

EOF

}
