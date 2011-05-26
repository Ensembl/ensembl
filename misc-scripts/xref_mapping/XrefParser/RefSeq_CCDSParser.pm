package XrefParser::RefSeq_CCDSParser;

use strict;

use DBI;

use base qw( XrefParser::BaseParser );

# Parse file of Refseq records and assign direct xrefs

sub run_script {

  my $self = shift if (defined(caller(1)));
  my $file = shift;
  my $source_id  = shift;
  my $species_id = shift;
  my $verbose    = shift;

  my $user = "ensro";
  my $host;
  my $port;
  my $dbname;
  my $pass;

  if($file =~ /host[=][>](\S+?)[,]/){
    $host = $1;
  }
  if($file =~ /port[=][>](\S+?)[,]/){
    $port =  $1;
  }
  if($file =~ /dbname[=][>](\S+?)[,]/){
    $dbname = $1;
  }
  if($file =~ /pass[=][>](\S+?)[,]/){
    $pass = $1;
  }

  my $dna_pred = XrefParser::BaseParser->get_source_id_for_source_name("RefSeq_dna_predicted");

  # becouse the direct mapping have no descriptions etc
  # we have to steal these from the previous Refseq parser.

  my %label;
  my %version;
  my %description;

  my $dbi = $self->dbi();  
  my $sql = "select xref.accession, xref.label, xref.version,  xref.description from xref, source where xref.source_id = source.source_id and source.name = 'RefSeq_dna'";
  my $sth = $dbi->prepare($sql);
  $sth->execute();
  my ($acc, $lab, $ver, $desc);
  $sth->bind_columns(\$acc, \$lab, \$ver, \$desc);
  while (my @row = $sth->fetchrow_array()) {
    $label{$acc} = $lab;
    $version{$acc} = $ver;
    $description{$acc} = $desc;
  }
  $sth->finish;
 



  $sql = 'select x.accession, x.xref_id, d.ensembl_stable_id, "Transcript"
            from xref x, transcript_direct_xref d, source s 
             where s.source_id = x.source_id and 
                   x.xref_id = d.general_xref_id and s.name like "CCDS"'; 
 
  $sth = $dbi->prepare($sql);
  $sth->execute();
  my ($access, $old_xref_id, $stable_id, $type);
  $sth->bind_columns(\$access, \$old_xref_id, \$stable_id, \$type);
  my %ensembl_stable_id;
  my %ensembl_type;
  my %old_xref;
  while (my @row = $sth->fetchrow_array()) {
      push @{$ensembl_stable_id{$access}}, $stable_id;
      $ensembl_type{$access} = $type;
      $old_xref{$access} = $old_xref_id; 
  }
  $sth->finish;
  
 

  my $line_count = 0;
  my $xref_count = 0;
  my $direct_count = 0;
  my %seen;
  my %old_to_new;

#
# dbi2 is the ccds database
#
  my $dbi2 = $self->dbi2($host, $port, $user, $dbname, $pass);
  if(!defined($dbi2)){
    return 1;
  }


##############################NEW#########################################

  # get ccds -> xref transcript_id                 ensembl_stable_id{CCDS1} = ENST00001
  # get ccds -> internal transcript_id             ccds_to_internal_id(CCDS1} = 12345

  $sql = 'select x.dbprimary_acc, ox.ensembl_id from xref x, object_xref ox, external_db e where x.xref_id = ox.xref_id and x.external_db_id = e.external_db_id and e.db_name like ? order by x.version';

# order by version added so that the hash gets overwritten with the latest version.


  # calculate internal_id -> xref transcript_id
  my %internal_to_stable_id;
  my ($acc, $internal_id);

  my $sth = $dbi2->prepare($sql); 
  $sth->execute("CCDS") or croak( $dbi2->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
    my $acc = $row[0];
    my $internal_id = $row[1];
    if(defined($ensembl_stable_id{$acc})){
      $internal_to_stable_id{$internal_id} =  $ensembl_stable_id{$acc};
    }
    else{
      print "$acc not found in ccds database????\n";
    }
  }    

  # for each object_xref for refseq_dna change internal_id to xref transcript_id
  $sth->execute("Refseq_dna") or croak( $dbi2->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
    my $refseq = $row[0];
    my $internal_id = $row[1];

    if(defined($internal_to_stable_id{$internal_id})){
    }
    else{
      print "Problem no internal_to_stable_id for $internal_id\n"; 
      next;
    }
  
    $line_count++;
    if(!defined($seen{$refseq})){
      $seen{$refseq} = 1;
      my $new_source_id = $source_id;
      if($refseq =~ /^XM/){
	$new_source_id = $dna_pred;
      }
      my $xref_id = $self->add_xref($refseq, $version{$refseq} , $label{$refseq}||$refseq , 
				    $description{$refseq}, $new_source_id, $species_id, "DIRECT");


      foreach my $stable_id (@{$internal_to_stable_id{$internal_id}}){
	$self->add_direct_xref($xref_id, $stable_id, "Transcript", "");
	$direct_count++;
      }

      $old_to_new{$old_xref{$refseq}} = $xref_id;
      $xref_count++;
    }
  }
############################END NEW######################################

  #for each one seen get all its dependent xrefs and load them fro the new one too;

  my $add_dependent_xref_sth = $dbi->prepare("INSERT INTO dependent_xref VALUES(?,?,?,?)");
  my $get_dependent_xref_sth = $dbi->prepare("SELECT dependent_xref_id, linkage_annotation "
					    .  "FROM  dependent_xref where master_xref_id = ?");

  foreach my $old_xref (keys %old_to_new){
      my $linkage;
      my $dependent_id;
      $get_dependent_xref_sth->execute($old_xref);
      $get_dependent_xref_sth->bind_columns(\$dependent_id, \$linkage);
      while(my @row = $get_dependent_xref_sth->fetchrow_array()){
	  $add_dependent_xref_sth->execute($old_to_new{$old_xref}, $dependent_id, $linkage, $source_id); 
      }   
  }


  print "Parsed $line_count RefSeq_dna identifiers from $file, added $xref_count xrefs and $direct_count direct_xrefs  from $line_count lines.\n" if ($verbose);


  return 0;

}

1;
