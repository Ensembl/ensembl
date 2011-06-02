#!/usr/bin/env perl

use strict;
use warnings;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Getopt::Long;

my ( $host, $user, $pass, $port, $dbname, $pattern, $genestats, $snpstats  );

GetOptions( "host=s", \$host,
	    "user=s", \$user,
	    "pass=s", \$pass,
	    "port=i", \$port,
	    "dbname=s", \$dbname,
	    "pattern=s", \$pattern,
	    "genestats", \$genestats,
	    "snpstats", \$snpstats
	  );


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
		     'lincRNA'              => 'lincRNA',);

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

foreach my $name (@dbnames) {
  if ( $pattern && ($name !~ /$pattern/) ) { next }

  printf( "\nConnecting to '%s'\n", $name );

  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $host,
					      -user => $user,
					      -port => $port,
					      -pass => $pass,
					      -dbname => $name);

  
  # do both genestats and snpstats by default
  $genestats = $snpstats = 1 if(!$genestats && !$snpstats);

 # $snpstats = 0 if ($name =~ /homo|mus/);


  # delete old attributes before starting
  foreach my $code (values %attrib_codes) {
    if ($genestats) {
      my $sth = $db->dbc()->prepare( "DELETE sa FROM seq_region_attrib sa, attrib_type at WHERE at.attrib_type_id=sa.attrib_type_id AND at.code=?" );
      $sth->execute("GeneNo_$code");
    }
    if ($snpstats) {
      my $sth = $db->dbc()->prepare( "DELETE sa FROM seq_region_attrib sa, attrib_type at WHERE at.attrib_type_id=sa.attrib_type_id AND at.code=?" );
      $sth->execute("SNPCount");
    }
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
  
  my $snp_db = variation_attach( $db );

  my $snps_present = $snpstats && $snp_db;

  my $slice_adaptor = $db->get_SliceAdaptor();
  my $attrib_adaptor = $db->get_AttributeAdaptor();

# Do not include non-reference sequences ie. haplotypes for human
#my $top_slices = $slice_adaptor->fetch_all( "toplevel" , undef, 1);
  my $top_slices = $slice_adaptor->fetch_all( "toplevel" );

  while (my $slice = shift(@{$top_slices})) {
#    print STDERR "Processing seq_region ", $slice->seq_region_name(), "\n";
    
    my @attribs;

    if($genes_present) {
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

	# not used:
	# my $no_space = $biotype;
	# $no_space =~ s/ /_/g;
	
	push @attribs, Bio::EnsEMBL::Attribute->new
	  (-NAME => $biotype.' Gene Count',
	   -CODE => 'GeneNo_'.$attrib_code,
	   -VALUE => $counts{$biotype},
	   -DESCRIPTION => 'Number of '.$biotype.' Genes');
      }
    }

    if( $snps_present ) {
	  my $sth = $snp_db->dbc->prepare("SELECT COUNT(*) FROM variation_feature WHERE seq_region_id = ?");
	  $sth->execute($slice->get_seq_region_id);
	  my $count;
	  $sth->bind_columns($count);
	  $sth->fetch;
	  
      push @attribs, Bio::EnsEMBL::Attribute->new
	(-NAME => 'SNP Count',
	 -CODE => 'SNPCount',
	 -VALUE => $count,
	 -DESCRIPTION => 'Total Number of SNPs');
	
	  $sth->finish;
    }

    $attrib_adaptor->store_on_Slice($slice, \@attribs);
#  print_chromo_stats([$slice]);
  }
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
