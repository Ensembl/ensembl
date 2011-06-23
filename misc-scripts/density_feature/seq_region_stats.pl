#!/usr/bin/env perl

use strict;
use warnings;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Getopt::Long;

my ( $host, $user, $pass, $port, $dbname, $pattern, $stats );

GetOptions( "host|h=s", \$host,
	    "user|u=s", \$user,
	    "pass|p=s", \$pass,
	    "port=i", \$port,
	    "dbname|d=s", \$dbname,
	    "pattern=s", \$pattern,
	    "stats|s=s", \$stats,
	    "help" ,               \&usage
	  );

usage() if (!$host || !$user || !$pass || (!$dbname && !$pattern) || !$stats || $stats !~ /^(gene|snp)$/ );


my %attrib_codes = ( 'miRNA'                => 'miRNA',
		     'snRNA'                => 'snRNA',
		     'snoRNA'               => 'snoRNA',
		     'rRNA'                 => 'rRNA',
		     'tRNA'                 => 'tRNA',
		     'snlRNA'               => 'snlRNA',
		     'known protein_coding' => 'knwCod',
		     'misc_RNA'             => 'mscRNA',
		     'novel protein_coding' => 'novCod',
		     'pseudogene'           => 'pseudo',
		     'scRNA'                => 'scRNA',
		     'Mt_tRNA'              => 'MTtRNA',
		     'Mt_rRNA'              => 'MTrRNA',
                     'ncRNA'                => 'ncRNA',
		     '3prime_overlapping_ncrna' => 'ncRNA', # (v63)
		     'polymorphic_pseudogene' => 'pseudo',
                     'havana_pseudogene'    => 'pseudo', 
                     'processed_pseudogene' => 'pseudo', 
                     'unprocessed_pseudogene' => 'pseudo', 
                     'Pseudogene'           => 'pseudo', 
                     'transcribed_pseudogene' => 'pseudo', # added for v54, else is flagged as "Unspecified biotype"
		     'scRNA_pseudogene'     => 'RNA_pseu',
		     'tRNA_pseudogene'      => 'RNA_pseu',
		     'rRNA_pseudogene'      => 'RNA_pseu',
		     'snoRNA_pseudogene'    => 'RNA_pseu',
		     'snRNA_pseudogene'     => 'RNA_pseu',
		     'misc_RNA_pseudogene'  => 'RNA_pseu',
		     'miRNA_pseudogene'     => 'RNA_pseu',
		     'Mt_tRNA_pseudogene'   => 'RNA_pseu',
		     'IG_V_gene'            => 'Ig',
		     'IG_J_gene'            => 'Ig',
		     'IG_D_gene'            => 'Ig',
		     'IG_C_gene'            => 'Ig',
		     'IG_Z_gene'            => 'Ig',
		     'IG_M_gene'            => 'Ig',
		     'TR_C_gene'            => 'Ig', #Actually it is TR but was put with Igs (v62)
		     'TR_J_gene'            => 'Ig', #Actually it is TR but was put with Igs (v62)
		     'TR_V_gene'            => 'Ig', #Actually it is TR but was put with Igs (v62)
		     'C_segment'            => 'Ig',
		     'D_segment'            => 'Ig',
		     'J_segment'            => 'Ig',
		     'V_segment'            => 'Ig',
                     'IG_pseudogene'        => 'pseudo',
		     'IG_V_pseudogene'      => 'pseudo',
		     'IG_C_pseudogene'      => 'pseudo',
		     'TR_V_pseudogene'      => 'pseudo',
		     'IG_J_pseudogene'      => 'pseudo',
		     'retrotransposed'      => 'rettran',
		     'processed_transcript' => 'proc_tr',
		     'lincRNA'              => 'lincRNA');



#get biotypes from the production database when new field attr_code is added to the biotype table

# Master database location:
# my ( $mhost, $mport ) = ( 'ens-staging1', '3306' );
# my ( $muser, $mpass ) = ( 'ensro',        undef );
# my $mdbname = 'ensembl_production';


# my $prod_dsn = sprintf( 'DBI:mysql:host=%s;port=%d;database=%s',
  #                   $mhost, $mport, $mdbname );
# my $prod_dbh = DBI->connect( $prod_dsn, $muser, $mpass,
  #                        { 'PrintError' => 1, 'RaiseError' => 1 } );

#my @attrib_codes  = map { @_[0]->@_[1] }  @{ $prod_dbh->selectall_arrayref('select distinct name, attr_code from biotype where is_current = 1 order by name') };

#my %attrib_codes = map { $_=>$_} @attrib_codes;

my @dbnames;
if (! $dbname) {
  my $dsn = sprintf( 'dbi:mysql:host=%s;port=%d', $host, $port );
  my $dbh = DBI->connect( $dsn, $user, $pass );
  @dbnames =
    map { $_->[0] } @{ $dbh->selectall_arrayref('SHOW DATABASES') };
}
else {
  @dbnames = ( $dbname )
}

my $genestats = 1 if($stats eq 'gene');
my $snpstats = 1 if($stats eq 'snp');

foreach my $name (@dbnames) {
  if ( $pattern && ($name !~ /$pattern/) ) { next }

  printf( "\nConnecting to '%s'\n", $name );

  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $host,
					      -user => $user,
					      -port => $port,
					      -pass => $pass,
					      -dbname => $name);

  my $total_count = 0;
  # delete old attributes before starting
  if ($genestats) {
      foreach my $code (values %attrib_codes) {
	  my $sth = $db->dbc()->prepare( "DELETE sa FROM seq_region_attrib sa, attrib_type at WHERE at.attrib_type_id=sa.attrib_type_id AND at.code=?" );
	  $sth->execute("GeneNo_$code");
      }
  }
  
  if ($snpstats) {
      my $sth = $db->dbc()->prepare( "DELETE sa FROM seq_region_attrib sa, attrib_type at WHERE at.attrib_type_id=sa.attrib_type_id AND at.code=?" );
      $sth->execute("SNPCount");    
  }

#
# Only run on database with genes
#

  my $genes_present;

  if($genestats) {
    my $sth = $db->dbc()->prepare( "select count(*) from gene" );
    $sth->execute();
    
    my ( $gene_count )  = $sth->fetchrow_array();
    
    $genes_present = ($gene_count) ? 1 : 0;
  } else {
    $genes_present = 0;
  }
  
#
# and seq_regions
#
  my $sth = $db->dbc()->prepare( "select count(*) from seq_region" );
  $sth->execute();
  my ( $seq_region_count ) = $sth->fetchrow_array();
  if( ! $seq_region_count ) {
    print STDERR "No seq_regions for $dbname.\n";
    exit();
  }
  
  my $snps_present;
  my $snp_db;

  if ($snpstats) {
      $snp_db = variation_attach( $db );
      if (defined $snp_db) {$snps_present = 1;}
  }

  my $slice_adaptor = $db->get_SliceAdaptor();
  my $attrib_adaptor = $db->get_AttributeAdaptor();

# Do not include non-reference sequences ie. haplotypes for human
#my $top_slices = $slice_adaptor->fetch_all( "toplevel" , undef, 1);
  my $top_slices = $slice_adaptor->fetch_all( "toplevel" );

  while (my $slice = shift(@{$top_slices})) {
#    print STDERR "Processing seq_region ", $slice->seq_region_name(), "\n";
      
    my @attribs;
   
    if($genes_present) {
      
      my %attrib_counts;
      my %counts;
      
      my $genes = $slice->get_all_Genes();
    
      while (my $gene = shift(@{$genes})) {

	my $biotype = $gene->biotype();
	if( $biotype =~ /coding/i ) {
	  if($gene->is_known()) {
	    $biotype = "known ".$biotype;
	  } else {
	    $biotype = "novel ".$biotype;
	  }
	}

	$counts{$biotype}++;

      }

      for my $biotype ( keys %counts ) {
	my $attrib_code = $attrib_codes{$biotype};
	if( !$attrib_code ) {
	  print STDERR "Unspecified biotype \"$biotype\" in database $name.\n";
	  next;
	}

	$attrib_counts{$attrib_code} += $counts{$biotype};
	
      }

      foreach my $attrib_code (keys %attrib_counts) {
	push @attribs, Bio::EnsEMBL::Attribute->new
	  (-NAME => $attrib_code.' Gene Count',
	   -CODE => 'GeneNo_'.$attrib_code,
	   -VALUE => $attrib_counts{$attrib_code},
	   -DESCRIPTION => 'Number of '.$attrib_code.' Genes');

      }

    }

    if( $snps_present ) {
	  my $sth = $snp_db->dbc->prepare("SELECT COUNT(*) FROM variation_feature WHERE seq_region_id = ?");
	  $sth->execute($slice->get_seq_region_id);
	  my $count;
	  $sth->bind_columns(undef,\$count);
	  $sth->fetch;
	  
      push @attribs, Bio::EnsEMBL::Attribute->new
	(-NAME => 'SNP Count',
	 -CODE => 'SNPCount',
	 -VALUE => $count,
	 -DESCRIPTION => 'Total Number of SNPs');
	
	  $sth->finish;
    }

    $attrib_adaptor->store_on_Slice($slice, \@attribs);
    my $slice_attrib_count = @attribs;
    $total_count += $slice_attrib_count;
#  print_chromo_stats([$slice]);
  }

  print STDOUT "Written $total_count seq reqion attributes to database $name on server $host.\n";

}



sub print_chromo_stats {
  my $chromosomes = shift;

  foreach my $chr (@$chromosomes) {
    print "\nchromosome: ",$chr->seq_region_name(),"\n";
    foreach my $attrib (@{$chr->get_all_Attributes()}) {
      print "  ", $attrib->name(), ": ", $attrib->value(), "\n";
    }
  }
}


#
# tries to attach variation database.
#

sub variation_attach {
  my $db = shift;

  my $core_db_name;
  $core_db_name = $db->dbc->dbname();
  if( $core_db_name !~ /_core_/ ) {
    return 0;
  }
  #
  # get a lost of all databases on that server
  #
  my $sth = $db->dbc->prepare( "show databases" );
  $sth->execute();
  my $all_db_names = $sth->fetchall_arrayref();
  my %all_db_names = map {( $_->[0] , 1)} @$all_db_names;
  my $snp_db_name = $core_db_name;
  $snp_db_name =~ s/_core_/_variation_/;
 


if( ! exists $all_db_names{ $snp_db_name } ) {
   return 0;
 }

 # this should register the dbadaptor with the Registry
 my $snp_db = Bio::EnsEMBL::Variation::DBSQL::DBAdaptor->new
   ( -host => $db->dbc()->host(),
     -user => $db->dbc()->username(),
     -pass => $db->dbc()->password(),
     -port => $db->dbc()->port(),
     -dbname => $snp_db_name,
     -group => "variation",
     -species => "DEFAULT"
   );

  return $snp_db;
}


sub usage {
  my $indent = ' ' x length($0);
  print <<EOF; exit(0);


For each toplevel slice, count the number of genes for each biotype
(gene stats) or count the number of SNPs (snp stats).

gene stats
What does it do?

Deletes all seq_region_attrib that have attrib_type code 
with a prefix 'GeneNo_'. All toplevel slices are fetched.

Input data: dna seqence, genes, xrefs, xref projections 
Output tables: seq_region_attrib (attrib_type code with prefix 'GeneNo')


When to run it in the release cycle?

After core have finished xref projections


Which databases to run it on?

Run on all core databases (including otherfeatures, cdna etc) for each release.


How long does it take?

It takes about 10 mins to run for a database in normal queue,


snp stats

What does it do?

Deletes out all seq_region_attrib that have attrib_type code of 'SNPCount'. 
Attach variation db if exists. All toplevel slices are fetched. 
For each slice, count the number of SNPs.

This option requires ensembl-variation in perl5lib.

Input data: top level seq regions, variation db
Output tables: seq_region_attrib (attrib_type code with prefix 'SNPCount')


When to run it in the release cycle?

When variation dbs have been handed over


Which databases to run it on?

Run on core databases only for new species or if the assembly changed, 
or if the variation positions have changed in the corresponding variation db.


How long does it take?

It takes about 20 mins to run for a database in normal queue.



Usage: 

  $0 -h host [-port port] -u user -p password \\
  $indent -d database | -pattern pattern \\
  $indent -s gene | snp  \\
  $indent [-help]  \\

  -h|host             Database host to connect to

  -port               Database port to connect to (default 3306)

  -u|user             Database username

  -p|pass             Password for user

  -d|dbname           Database name

  -pattern            Database name regexp

  -s|stats            'gene' or 'snp'

  -help               This message


EOF
 
}
