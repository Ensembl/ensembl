#!/usr/local/ensembl/bin/perl -w
=head1 NAME

  delete_genes.pl

=head1 SYNOPSIS
 
  delete_genes.pl
  deletes genes from given database whose ids are passed in through STDIN

=head1 DESCRIPTION


=head1 OPTIONS

    -host      host name for database (gets put as host= in locator)

    -port      For RDBs, what port to connect to (port= in locator)

    -dbname    For RDBs, what name to connect to (dbname= in locator)

    -dbuser    For RDBs, what username to connect as (dbuser= in locator)

    -dbpass    For RDBs, what password to use (dbpass= in locator)

    -genefile  File containing list of genes to be deleted

=cut

use strict;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;


my $host   = 'ecs2d';
my $port;
my $dbname = 'anopheles_gambiae_core_13_2';
my $dbuser = 'ensro';
my $pass   = undef;

#IPR002104,IPR001584,IPR001037,IPR003308,IPR000477,IPR000123

&GetOptions( 
	     'host:s'      => \$host,
	     'port:n'      => \$port,
	     'dbname:s'    => \$dbname,
	     'dbuser:s'    => \$dbuser,
	     'pass:s'      => \$pass,
	     );

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $host,
					    -user => $dbuser,
					    -dbname => $dbname,
					    -pass  => $pass
					   );

my $transadaptor = $db->get_TranscriptAdaptor();

my $query = "select distinct(f.translation_id),x.description from protein_feature f, interpro i, xref x where f.hit_id = i.id and i.interpro_ac = x.dbprimary_acc and (x.description like '%transcriptase%' or x.description like '%integrase%' or x.description like '%bact%')";

my $sth = $db->prepare($query);
$sth->execute();

while (my ($ac,$desc) = $sth->fetchrow()) {
    my $query1 = "select transcript_id from transcript where translation_id = $ac";
    my $sth1 = $db->prepare($query1);
    $sth1->execute();
    my $translac = $sth1->fetchrow();
    my $obj = $transadaptor->fetch_by_dbID($translac);
    my $seq = $obj->seq->seq;
    print ">$translac#transposon\t$desc\n$seq\n";
}






