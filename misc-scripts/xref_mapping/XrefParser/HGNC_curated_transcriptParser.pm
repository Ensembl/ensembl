package XrefParser::HGNC_curated_transcriptParser;

use strict;
use File::Basename;

use base qw( XrefParser::BaseParser );

use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;


#my $dbi2;

if (!defined(caller())) {

  if (scalar(@ARGV) != 1) {
    print STDERR "\nUsage: HGNC_curated_transcriptParser.pm file <source_id> <species_id>\n\n";
    exit(1);
  }

  run(@ARGV);
}

sub run_script {
  my $self = shift if (defined(caller(1)));

  my $file = shift;
  my $source_id = shift;
  my $species_id = shift;
  my $verbose = shift;

  my ($type, $my_args) = split(/:/,$file);
  
  my $user  ="ensro";
  my $host;
  my $port;
  my $dbname;
  my $pass;

  if($my_args =~ /host[=][>](\S+?)[,]/){
    $host = $1;
  }
  if($my_args =~ /port[=][>](\S+?)[,]/){
    $port =  $1;
  }
  if($my_args =~ /dbname[=][>](\S+?)[,]/){
    $dbname = $1;
  }
  if($my_args =~ /pass[=][>](\S+?)[,]/){
    $pass = $1;
  }
  if($my_args =~ /user[=][>](\S+?)[,]/){
    $user = $1;
  }

  print "Using $host $dbname for Vega\n" if($verbose);


  my $clone_source_id =
    $self->get_source_id_for_source_name('Clone_based_vega_transcript');

  my $hgnc_source_id =
    $self->get_source_id_for_source_name('HGNC','havana');
 
  my $sql = 'select tsi.stable_id, x.display_label from xref x, object_xref ox , transcript_stable_id tsi, external_db e where e.external_db_id = x.external_db_id and x.xref_id = ox.xref_id and tsi.transcript_id = ox.ensembl_id and e.db_name like ?';


  my %ott_to_vega_name;
  my %ott_to_enst;
  
  my $dbi2 = $self->dbi2($host, $port, $user, $dbname, $pass);
  if(!defined($dbi2)){
    return 1;
  }
  

  my $sth = $dbi2->prepare($sql);   # funny number instead of stable id ?????
  $sth->execute("Vega_transcript") or croak( $dbi2->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
    $ott_to_vega_name{$row[0]} = $row[1];
  }
  $sth->finish;

  $sth = $dbi2->prepare($sql);   # funny number instead of stable id ?????
  $sth->execute("ENST_CDS") or croak( $dbi2->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
    $ott_to_enst{$row[0]} = $row[1];
  }
  $sth->finish;

  $sth = $dbi2->prepare($sql);
  $sth->execute("ENST_ident") or croak( $dbi2->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
    $ott_to_enst{$row[0]} = $row[1];
  }
  $sth->finish;


  my $xref_count = 0;



  
  my $dbi = $self->dbi();

  my %synonym;
  my $dbname = "HGNC";
  my $syn;
  my $name;

  $sth = $dbi->prepare('select es.synonym, x.label from synonym es, xref x, source s where x.xref_id = es.xref_id and x.source_id = s.source_id and s.name = "EntrezGene"' );
  $sth->execute();
  $sth->bind_columns(\$syn,\$name);
  while($sth->fetch){
    $synonym{$syn} = $name;
  }
  $sth->finish;


  $sth = $dbi->prepare('select es.synonym, x.label from synonym es, xref x, source s where x.xref_id = es.xref_id and x.source_id = s.source_id and s.name = "'.$dbname.'" and s.priority_description like "desc_only"');
  $sth->execute();
  $sth->bind_columns(\$syn,\$name);
  while($sth->fetch){
    $synonym{$syn} = $name;
  }
  $sth->finish;



  #get the source ids for HGNC sources

  my (%accession, %version, %description);

  $sql  = 'select source_id from source where name like "HGNC" ';
  $sql .= 'and priority_description like "desc_only" ';
  $sth = $dbi->prepare($sql);
  
  $sth->execute();


  my ($hgnc_source_id);
  $sth->bind_columns(\$hgnc_source_id);
  my @arr;
  while($sth->fetch()){
    push @arr, $hgnc_source_id;
  }
  $sth->finish;

  $sql = "select accession, label, version,  description from xref where source_id in (".join(", ",@arr).")";
  $sth = $dbi->prepare($sql);
  $sth->execute();
  my ($acc, $lab, $ver, $desc);
  my $hgnc_loaded_count = 0;
  $sth->bind_columns(\$acc, \$lab, \$ver, \$desc);
  while (my @row = $sth->fetchrow_array()) {
    $accession{$lab} = $acc;
    $version{$lab} = $ver;
    $description{$lab} = $desc;
    $hgnc_loaded_count++;
  }
  $sth->finish;

  if($hgnc_loaded_count == 0){
    die "No point continuing no hgncs there\n";
  }










  my $not_in_hgnc = 0;
  foreach my $ott (keys %ott_to_enst){
    if(defined($ott_to_vega_name{$ott})){
      my $id = $hgnc_source_id;
      my $name  = $ott_to_vega_name{$ott};
      my $acc = undef;
      my $xref_id ;
      if($name =~ /[.]/){
	$id = $clone_source_id;
        $name =~ s/[.]\d+//;    #remove .number
	$xref_id = $self->add_xref($name, "" , $name , $description{$name}, $id, $species_id, "DIRECT");
      }
      else{ 
	my $copy = $name;
	$name =~ s/-\d+$//;    #remove -number
	if(defined($accession{$name})){
	}
	elsif(defined($synonym{$name})){
	  $name = $synonym{$name};
	  if(!defined($accession{$name})){
	    print "Havana name $copy which has a synonym of $name cannot be found in the HGNC data???\n";
            $not_in_hgnc++;
	    next;
	  }
          print "Havana uses old name $copy instead of $name\n";
	}
	else{
	  print "Havana name ($copy)  $name cannot be found in the HGNC data???\n";
          $not_in_hgnc++;
	  next;
	}
	$xref_id = $self->add_xref($accession{$name}, "" , $name , $description{$name}, $id, $species_id, "DIRECT");	
      }
      $xref_count++;
      
      $self->add_direct_xref($xref_id, $ott_to_enst{$ott}, "transcript", "");
    }
  }
  
  print "$xref_count direct xrefs succesfully parsed\n" if($verbose);
  print "$not_in_hgnc xrefs could not be loaded as they were not in HGNC\n)"  if($verbose and $not_in_hgnc);
  return 0;
}





1;

