#!/usr/local/bin/perl -w

=head1 NAME

Another data integrity test

=head1 SYNOPSIS

This script has 2 purposes:
- check that the length of each protein is not shorter than cut-off
given by user (default 2).
- check that protein sequence does not contain stops ('*').

How it all happens:
- get the list of "current" chromosomes (at the present we still have
things like chrUL_random) unless a list is given on the commandline
- build virtual contigs for each chromosome
- get genes on each virtual contig
- get transcripts for each gene
- get translation of each transcript
- check the length of the transcript and presence of '*' characters.

Output is divided to STDERR and STDOUT

STDERR gets messages like:
Now checking genes on chromosome chr1

...

159 genes with 172 transcripts

...

Checked 41 chromosomes with 335 genes with 351 transcripts 

Information about proteins, transcripts, genes shorter or equal to
cut-off or containing stops id printed to STDOUT:

ENSP00000211067 encoded by transcript ENST00000211067 from gene ENSG00000095908 on chromosome chrNA_random has length:  4

...

ENSP00000215395 encoded by transcript ENST00000215395 from gene ENSG00000099644 on chromosome chrUL_random contains '*' character(s).
TPSQSEDLRACFEQNKFQGIATRDGLALAIGFLEPIVQNWFQNERSRQVRQHCRESRPRPGRHGPQEGR*KRTAVTGSQTALLLRAFEKDRFPGIAAREDLAR*TGLPGSRIQIRFQNRRARHLGEAGRAPAKAGSRYNAAP  


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

=head1 WARNING

On big chromosomes (like chr1) one is likely to run out of memory
(that's what happend at least on ecs1c)

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
my @chromosomes;
if(@ARGV)
{
    @chromosomes = @ARGV;
}
else
{
    my $sth = $db->prepare
    (
        "select distinct(chr_name) from static_golden_path"
    );
    my $rv = $sth->execute();
    my $chr;
    $sth->bind_columns(undef, \$chr);
    while($sth->fetch)
    {
	push @chromosomes, $chr;
    }
}


# Get static golden path adaptor
my $stadp = $db->get_StaticGoldenPathAdaptor();


# Iterate over chromosomes and pull out all genes
my($genecount, $transcriptcount);
foreach my $chromosome (@chromosomes)
{
    my($chr_genecount, $chr_transcriptcount);
    print STDERR "Now checking genes on chromosome $chromosome\n";
    my $vc = $stadp->fetch_VirtualContig_by_chr_name($chromosome);
    my @genes = $vc->get_all_Genes;
    foreach my $gene (@genes)
    {
	$genecount++; $chr_genecount++; 
        foreach my $transcript ($gene->each_Transcript)
        {
	    $transcriptcount++; $chr_transcriptcount++;
            my $pep = $transcript->translate;
            my $length = length($pep->seq);
            if($length < $minlength)
            {
		print $pep->id,
                      " encoded by transcript ",
                      $transcript->id,
                      " from gene ",
                      $gene->id,
                      " on chromosome ",
                      $chromosome,
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
                      " on chromosome ",
                      $chromosome,
                      " contains '*' character(s).\n",
                      $pep->seq,
                      "\n";
            }
        } #  foreach my $transcript ($gene->each_Transcript)
    } # foreach my $gene (@genes)
    print STDERR "$chr_genecount genes with $chr_transcriptcount transcripts\n";
} # foreach my $chromosome (@chromosomes)

print STDERR "Checked ", scalar(@chromosomes), " chromosomes with $genecount genes with $transcriptcount transcripts\n";
