#!/usr/local/bin/perl -w

=head1 NAME

Another data integrity test

=head1 SYNOPSIS

This script has 2 purposes:
- check that the length of each protein is not shorter than cut-off
given by user (default 2).
- check that protein sequence does not contain stops ('*').

How it all happens:
- get the list of fpcs
- build virtual contigs for each fpc
- get genes on each virtual contig
- get transcripts for each gene
- get translation of each transcript
- check the length of the transcript and presence of '*' characters.

Output is divided to STDERR and STDOUT

STDERR gets messages like:
Now checking genes on fpc ctg12323
22 genes with 30 transcripts 


Information about proteins, transcripts, genes shorter than
cut-off or containing stops id printed to STDOUT:

ENSP00000228176 encoded by transcript ENST00000228176 from gene ENSG00000110736 on fpc ctg12269 contains '*' character(s).
VKRAYLVHSAYDQSYNFIYKSFRIASII*X  

...

ENSP00000227580 encoded by transcript ENST00000227580 from gene ENSG00000110163 on fpc ctg12475 has length:     1  


=head1 OPTIONS

-host       db server name (default localhost)

-port       port to connect to (default 3306)

-dbname     name of the database to use (default ensembl_test)

-driver     database driver (defauly mysql)

-user       username for database access (default ensro)

-pass       password for database access (default undef)

-minlength  minimal length of the protein to be considered OK

-h|help     print out help (this text)


=head1 WARNING

This script takes quite a long time to run...

=cut

use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long; 

my $minlength = 2;
my $host = 'localhost';
my $port = 3306;
my $dbname = 'ensembl_test';
my $pass = undef;
my $user = 'ensro';
my $driver = 'mysql';
my $help;

&GetOptions
(
             'host:s'      => \$host,
             'port:n'      => \$port,
             'dbname:s'    => \$dbname,
             'user:s'    => \$user,
             'pass:s'      => \$pass,
             'driver:s'    => \$driver,
             'minlength:n' => \$minlength,
             'h|help'      => \$help,
);    

$help && exec('perldoc', $0);                                                                                                                         

# Get db adaptor
my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new
(
    -user   => $user,
    -dbname => $dbname,
    -host   => $host,
    -driver => $driver,
    -port   => $port,
    -pass   => $pass,
);


# Get list of all current cromosomes
my @fpcs;
if(@ARGV)
{
    @fpcs = @ARGV;
}
else
{
    my $sth = $db->prepare
    (
        "select distinct(fpcctg_name) from static_golden_path"
    );
    my $rv = $sth->execute();
    my $fpc;
    $sth->bind_columns(undef, \$fpc);
    while($sth->fetch)
    {
	push @fpcs, $fpc;
    }
}


# Get static golden path adaptor
my $stadp = $db->get_StaticGoldenPathAdaptor();


# Iterate over fpcs and pull out all genes
my($genecount, $transcriptcount) = (0, 0);
foreach my $fpc (@fpcs)
{
    my($fpc_genecount, $fpc_transcriptcount) = (0, 0);
    print STDERR "Now checking genes on fpc $fpc\n";
    my $vc = $stadp->fetch_VirtualContig_by_fpc_name($fpc);
    my @genes = $vc->get_all_Genes;
    foreach my $gene (@genes)
    {
	$genecount++; $fpc_genecount++; 
        foreach my $transcript ($gene->each_Transcript)
        {
	    $transcriptcount++; $fpc_transcriptcount++;
            my $pep = $transcript->translate;
            my $length = length($pep->seq);
            if($length < $minlength)
            {
		print $pep->id,
                      " encoded by transcript ",
                      $transcript->id,
                      " from gene ",
                      $gene->id,
                      " on fpc ",
                      $fpc,
                      " has length:\t",
                      $length,
                      "\n";
            }
            if($pep->seq =~ /\*/)
            {
		print $pep->id,
                      " encoded by transcript ",
                      $transcript->id,
                      " from gene ",
                      $gene->id,
                      " on fpc ",
                      $fpc,
                      " contains '*' character(s).\n",
                      $pep->seq,
                      "\n";
            }
        } #  foreach my $transcript ($gene->each_Transcript)
    } # foreach my $gene (@genes)
    print STDERR "$fpc_genecount genes with $fpc_transcriptcount transcripts\n";
} # foreach my $fpc (@fpcs)

print STDERR "Checked ", scalar(@fpcs), " fpcs with $genecount genes with $transcriptcount transcripts\n";






