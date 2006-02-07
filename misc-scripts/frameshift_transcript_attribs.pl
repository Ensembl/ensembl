use strict;
use warnings;

# Finds all potential frameshifts (exons 1, 2 4 or 5 bp apart)
# in a database and adds transcript attributes for them.
# Attribute value is intron length.

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Attribute;

use Getopt::Long;

my ($host, $port, $user, $pass, $dbname, $nostore, $delete);

GetOptions('host=s'   => \$host,
           'user=s'   => \$user,
           'port=i'   => \$port,
           'pass=s'   => \$pass,
           'dbname=s' => \$dbname,
	   'nostore'  => \$nostore,
	   'delete'   => \$delete,
           'help'     => sub { usage(); exit(0); });

$port ||= 3306;

usage() if(!$user || !$dbname || !$host);

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(-host   => $host,
					     -user   => $user,
					     -pass   => $pass,
					     -dbname => $dbname,
					     -port   => $port);

if ($delete) {

  print STDERR "Deleting existing 'Frameshift' transcript attributes\n";
  my $sth = $db->dbc()->prepare("DELETE ta, at FROM transcript_attrib ta, attrib_type at WHERE at.attrib_type_id=ta.attrib_type_id AND at.code='Frameshift'");
  $sth->execute();

}

my $attribute_adaptor = $db->get_AttributeAdaptor();
my $transcript_adaptor = $db->get_TranscriptAdaptor();

my %biotypes = ();

print STDERR "Finding frameshifts in $dbname, creating transcript attributes ...\n";
print STDERR "Attributes will not be stored in database\n" if ($nostore);

my $sth = $db->dbc()->prepare
  (qq{SELECT t.transcript_id, g.biotype,
      MIN(IF(e1.seq_region_strand = 1,
	     e2.seq_region_start - e1.seq_region_end - 1,
	     e1.seq_region_start - e2.seq_region_end - 1)) AS intron_length
      FROM exon e1, exon e2, exon_transcript et1, exon_transcript et2,
      transcript t, gene g
      WHERE et1.exon_id = e1.exon_id
      AND et2.exon_id = e2.exon_id
      AND et1.transcript_id = et2.transcript_id
      AND et1.rank = et2.rank - 1
      AND et1.transcript_id = t.transcript_id
      AND t.gene_id = g.gene_id
      GROUP BY t.transcript_id
      HAVING intron_length IN (1,2,4,5)});

$sth->execute();

my ($transcript_id, $biotype, $intron_length, $count);
$sth->bind_columns(\$transcript_id, \$biotype, \$intron_length);

while($sth->fetch()){

  my $attribute = Bio::EnsEMBL::Attribute->new(-CODE => 'Frameshift',
					       -NAME => 'Frameshift',
					       -DESCRIPTION => 'Frameshift modelled as intron',
					       -VALUE => $intron_length);

  my @attribs = ($attribute);

  my $transcript = $transcript_adaptor->fetch_by_dbID($transcript_id);

  $attribute_adaptor->store_on_Transcript($transcript, \@attribs) if (!$nostore);

  $biotypes{$biotype}++;
  $count++;

}

print "$count short intron attributes\n";
print "Attributes not stored in database\n" if ($nostore);

print "Biotypes of affected genes:\n";
foreach my $biotype(keys %biotypes) {
  print $biotype . "\t" . $biotypes{$biotype} . "\n";
}

# ----------------------------------------------------------------------

sub usage {

  print << "EOF";

  Finds all potential frameshifts (exons 1, 2 4 or 5 bp apart) in a database 
  and adds transcript  attributes for them. Attribute value is intron length.

  perl introns_to_transcript_attribs.pl {options}

 Options ([..] indicates optional):

   --host       The database server to connect to.

   [--port]     The port to use. Defaults to 3306.

   --user       Database username. Must allow writing.

   --pass       Password for user.

   --dbname     The database to use.

  [--nostore]   Don't store the attributes, just print results.

  [--delete]    Delete any existing "Frameshift" attributes before creating new ones.

  [--help]      This text.


EOF

}

# ----------------------------------------------------------------------
