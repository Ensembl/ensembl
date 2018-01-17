=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

package XrefMapper::BasicMapper;

use strict;
use warnings;
use Carp;
use Cwd;
use File::Basename;
use IPC::Open3;

use XrefMapper::db;
use Bio::EnsEMBL::DBSQL::AltAlleleGroupAdaptor;

=head2 new

  Description: Constructor for BasicMapper.
  Returntype : BasicMapper
  Exceptions : none
  Caller     : general

=cut

sub new{
  my($class, @args) = @_;

  my $self ={};
  bless $self,$class;
  return $self;
}



=head2 xref

  Arg [1]    : (optional)
  Example    : $mapper->core($new_core);
  Description: Getter / Setter for the core.
               info for the xref database.
  Returntype : XrefMapper::db
  Exceptions : none

=cut

sub xref{
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_xref} = $arg );
  return $self->{_xref};
}

=head2 farm_queue

  Arg [1]    : (optional)
  Example    : $mapper->farm_queue("long");
  Description: Getter / Setter for the farm queue.
  Returntype : string
  Exceptions : none

=cut

sub farm_queue{
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_queue} = $arg );
  return $self->{_queue};
}

=head2 exonerate

  Arg [1]    : (optional)
  Example    : $mapper->exonerate("/usr/local/exonerate1.1.1");
  Description: Getter / Setter for the exonerate executable with full path.
  Returntype : string
  Exceptions : none

=cut

sub exonerate{
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_exonerate} = $arg );
  return $self->{_exonerate};
}

=head2 core

  Arg [1]    : (optional)
  Example    : $mapper->core($new_core);
  Description: Getter / Setter for the core.
               info for the ensembl core database.
  Returntype : XrefMapper::db
  Exceptions : none

=cut

sub core{
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_core} = $arg );
  return $self->{_core};
}

=head2 previous_core

  Arg [1]    : (optional)
  Example    : $mapper->previous_core($old_core);
  Description: Getter / Setter for the previous release of the core db.
  Returntype : XrefMapper::db
  Exceptions : none

=cut

sub previous_core{
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_previous_core} = $arg );
  return $self->{_previous_core};
}


=head2 add_meta_pair

  Arg [1]    : key
  Arg [2]    : value
  Example    : $mapper->add_meta_pair("head_directory","/lustre/src/");
  Description: Adds key value pairs to the database
  Returntype : none
  Exceptions : none

=cut

sub add_meta_pair {

  my ($self, $key, $value) = @_;

  my $sth = $self->xref->dbc->prepare('insert into meta (meta_key, meta_value, date) values("'.$key.'", "'.$value.'", now())');
  $sth->execute;
  $sth->finish;
  return;
}

sub update_process_status{
  my ($self, $value) = @_;

  my $sth_stat = $self->xref->dbc->prepare("insert into process_status (status, date) values('".$value."',now())");
  $sth_stat->execute();
  $sth_stat->finish;
  return;
}

sub xref_latest_status { 
  my $self = shift;
  my $verbose = shift || 0;

  my $sth = $self->xref->dbc->prepare("select id, status, date from process_status order by id");

  $sth->execute();
  my ($id, $status, $date);
  $sth->bind_columns(\$id, \$status,\$date);
  while($sth->fetch){
    print "$status\t$date\n" if($verbose and $self->verbose);
  }
  return $status;

}

sub get_meta_value {
  my ($self, $key) = @_;

  my $sth = $self->xref->dbc->prepare('select meta_value from meta where meta_key like "'.$key.'" order by meta_id');

  $sth->execute();
  my $value;
  $sth->bind_columns(\$value);
  while($sth->fetch){   # get the last one
  }
  $sth->finish;

  return $value;
}

sub process_file {
  my $self = shift;
  my $file = shift;
  my $verbose = shift;
  

  my $xref=undef;
  my $ensembl=undef;
  my $type;
  
  my %xref_hash=();
  my %species_hash=();
  my %farm_hash=();
  
  open my $fh, "<", $file or croak ("\nCannot open input file '$file':\n $!\n");
  while( my $line = <$fh> ) {

    chomp($line);
    next if $line =~ /^#/;
    next if !$line;

    my ($key, $value) = split("=",$line);
    $value =~ s/^\s*// if defined $value;
    $value =~ s/\s*$// if defined $value;
    $key =~ s/^\s*// if defined $key;
    $key =~ s/\s*$// if defined $key;
    if($key eq "species"){
      $type = "species";
      $species_hash{'species'} = $value;
    }
    elsif($key eq "xref"){
      $type = "xref";
    }
    elsif($key eq "farm"){
      $type = "farm";
    }
    elsif($type eq "species"){ # processing species data
      $species_hash{lc($key)} = $value;
    }
    elsif($type eq "xref"){    # processing xref data
      $xref_hash{lc($key)} = $value;
    }
    elsif($type eq "farm"){
      $farm_hash{lc($key)} = $value;
    }
  }
  close $fh or croak "Can't close file";

  my $value = $species_hash{'species'};
  my $taxon = $species_hash{'taxon'};

  if ($value !~ /_/) {
    print STDERR "\'$value\' is not a recognised species - please use full species name (e.g. homo_sapiens) in $file\n";
    exit(1);
  }

  my $use_basic = 0;
  my $mapper;
  my $module;
  my $class = "XrefMapper/$value.pm";
  my $eval_test = eval {
    require $class;
  };
  if($@ or $eval_test != 1) {
    if ($@ =~ /Can\'t locate $class/) {
      if (defined $taxon) {
      	$class = "XrefMapper/$taxon.pm";
      	eval {
	  require $class;
      	};
      	if($@) {
	  if ($@ =~ /Can\'t locate $class/)  {
	    $use_basic = 1;
	  } else { die "$@"; }
       	} else {
	  $module = $taxon; 
       	}
      }
      else {
	$use_basic = 1;
      }
    }
    else {
      die "$@";
    }

  } else{
    $module = $value;
  }

  if ($use_basic or !defined $module) {
	if(defined($verbose) and $verbose) {
		my $warning_msg = "Did not find a specific mapping module XrefMapper::$value ";
		if (defined $taxon) {
			$warning_msg .= "or XrefMapper::$taxon "; 
		}
		$warning_msg .= "- using XrefMapper::BasicMapper instead\n";
 		carp($warning_msg);
	}
        require XrefMapper::BasicMapper;
        $module = "BasicMapper";
  }

  $mapper = "XrefMapper::$module"->new();

  if(defined($farm_hash{'queue'})){
    $mapper->farm_queue($farm_hash{'queue'});
  }
  if(defined($farm_hash{'exonerate'})){
    $mapper->exonerate($farm_hash{'exonerate'});
  }


  if(defined($xref_hash{host}) ){
    my ($host, $user, $dbname, $pass, $port);
    $host = $xref_hash{'host'};
    $user = $xref_hash{'user'};
    $dbname = $xref_hash{'dbname'};
    if(defined($xref_hash{'password'})){
      $pass = $xref_hash{'password'};
    }
    else{
      $pass = '';
    }
    if(defined($xref_hash{'port'})){
      $port = $xref_hash{'port'};
    }
    else{
      $port = 3306;
    }

    $xref = new XrefMapper::db(-host => $host,
			       -port => $port,
			       -user => $user,
			       -pass => $pass,
			       -group   => 'core',
			       -dbname => $dbname);

    $mapper->xref($xref);
    $mapper->add_meta_pair("xref", $host.":".$dbname);
    if(defined($xref_hash{'dir'})){
      $xref->dir($xref_hash{'dir'});
      if(!-d $xref_hash{'dir'}){
	croak "directory ".$xref_hash{'dir'}." does not exist please create this\n";
      }
    }
    else{
      croak "No directory specified for the xref fasta files\n";
    }

  }
  else {
    croak "No host name given for xref database\n";
  }
 
  if(defined($species_hash{'species'})){

    my ($host, $port, $user, $dbname, $pass);
    $host = $species_hash{'host'};
    $user = $species_hash{'user'};
    $dbname = $species_hash{'dbname'};
    if(defined($species_hash{'password'})){
      $pass = $species_hash{'password'};
    }
    else{
      $pass = '';
    }
    if(defined($species_hash{'port'})){
      $port = $species_hash{'port'};
    }
    else{
      $port = '';
    }

    my $core = new XrefMapper::db(-host => $host,
				  -port => $port,
				  -user => $user,
				  -pass => $pass,
				  -group   => 'core',
				  -dbname => $dbname);
    
    $mapper->core($core);
    
    $mapper->add_meta_pair("species", $host.":".$dbname);

    if(defined($species_hash{'dir'})){
       $core->dir($species_hash{'dir'});
       if(!-d $species_hash{'dir'}){
	  croak "directory ".$species_hash{'dir'}." does not exist please create this\n";
       }
    }    
    else{
       croak "No directory specified for the ensembl fasta files\n";
    }
    
    $core->species($value);

    #connect to previous release of core db if connection details specified in xref_input (pr_host, pr_port, pr_dbname, pr_user) 
    if (defined( $species_hash{'pr_host'}) && defined( $species_hash{'pr_user'}) && defined( $species_hash{'pr_dbname'}) ) {
	my ($pr_host, $pr_port, $pr_user, $pr_dbname);
	$pr_host = $species_hash{'pr_host'};
	$pr_user = $species_hash{'pr_user'};
	$pr_dbname = $species_hash{'pr_dbname'};
	if(defined($species_hash{'pr_port'})){
	    $pr_port = $species_hash{'pr_port'};
	}

	my $previous_core = new XrefMapper::db(-host => $pr_host,
				  -port => $pr_port,
				  -user => $pr_user,
				  -pass => '',
				  -group   => 'core',
				  -dbname => $pr_dbname);
    
	$mapper->previous_core($previous_core);

	$mapper->add_meta_pair("species", $pr_host.":".$pr_dbname);    

    }
  }

  return $mapper;
}


=head2 dumpcheck

  Arg [1]    : (optional)
  Example    : $mapper->dumpcheck("yes");
  Description: Getter / Setter for dumpcheck.
               If set the mapper will not dump fasta files
               if they exist already.
  Returntype : scalar
  Exceptions : none

=cut

sub dumpcheck {
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_dumpcheck} = $arg );
  return $self->{_dumpcheck};
}

sub nofarm {
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_nofarm} = $arg );
  return $self->{_nofarm};
}

sub verbose {
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_verbose} = $arg );
  return $self->{_verbose};
}

sub species_id {
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_species_id} = $arg );
  return $self->{_species_id};
}

sub get_id_from_species_name {
 my ($self, $species_name) = @_;

 my $sql = "select species_id from species where name = '".$species_name."'";
 my $sth = $self->xref->dbc->prepare($sql);
 $sth->execute();
 my @row = $sth->fetchrow_array();
 my $species_id;
 if (@row) {
   $species_id = $row[0];
 } else {
   print STDERR "Couldn't get ID for species ".$species_name."\n";
   print STDERR "It must be one of :-\n";
   $sql = "select name from species";
   $sth = $self->xref->dbc->prepare($sql);
   $sth->execute();
   while(my @row2 = $sth->fetchrow_array()){
     print STDERR $row2[0]."\n";
   }
 }
 $sth->finish();
 
 return $species_id;
 

}

#
# Alt alleles
#

sub get_alt_alleles {
  my $self =  shift;
  
  my $dba = $self->core->dba;
  my $aaga = Bio::EnsEMBL::DBSQL::AltAlleleGroupAdaptor->new($dba);
  
  my $aa_list = $aaga->fetch_all();
  
  my $count = scalar(@$aa_list);
  my %alt_id_to_gene_id;
  my %gene_id_to_alt_id;
  my $max_alt_id = 0;
  my %is_reference;
  my $sth;
  my $insert_sth = $self->xref->dbc->prepare("insert into alt_allele (alt_allele_id, gene_id, is_reference) values (?, ?,?)");

  if($count){
    $sth = $self->xref->dbc->prepare("delete from alt_allele");
    $sth->execute;
    my $alt_added = 0;
    my $num_of_genes = 0;
    
    # Iterate through all alt-allele groups, pushing unique alleles into the xref alt allele table.
    # Track the reference gene IDs.
    
    foreach my $aag (@$aa_list) {
        my $ref_gene = $aag->rep_Gene_id();
        # Representative gene not guaranteed, try to find an alternative best fit
        if (!$ref_gene) {
            my $genes = $aag->get_all_Genes;
            foreach my $gene (@$genes) {
                if ($gene->slice->is_reference) {
                    $ref_gene = $gene->dbID;
                }
            }
        }
        if (!$ref_gene) { 
            warn('Tried very hard but failed to select a representative gene for alt-allele-group '.$aag->dbID);
            next;
        }
        $is_reference{$ref_gene} = 1;
        my $others = $aag->get_all_Gene_ids('no rep');
        # Extra step in place to handle non-ref situations
        my @cleaned_others = grep {!/$ref_gene/} @$others;
        
        $insert_sth->execute($aag->dbID,$ref_gene,1);
        $num_of_genes++;
        $alt_added++;
        foreach my $aa (@cleaned_others) {
            $insert_sth->execute($aag->dbID,$aa,0);
            $num_of_genes++;
        }
        
        if ($aag->dbID > $max_alt_id) { $max_alt_id = $aag->dbID }
    }
    
    print "$alt_added alleles found containing $num_of_genes genes\n";
  }
  else{
    print "No alt_alleles found for this species.\n" ;
  }


  ### LRGs added as alt_alleles in the XREF system but never added to core.

  #
  # Use $max_alt_id for new ones.
  #
  
  my $sql =(<<'LRG');
SELECT  ox.ensembl_id, g.gene_id
  FROM xref x, object_xref ox, external_db e, gene g
    WHERE x.xref_id = ox.xref_id AND
          e.external_db_id = x.external_db_id AND
          e.db_name like "Ens_Hs_gene" AND
          ox.ensembl_object_type = "Gene" AND
           x.display_label = g.stable_id
LRG
  
  $sth = $self->core->dbc->prepare($sql);
  my ($core_gene_id, $lrg_gene_id);
  $sth->execute();
  $sth->bind_columns(\$lrg_gene_id, \$core_gene_id);
  
  $count =0;
  
  my $old_count = 0;
  my $new_count = 0;
  my $lrg_count = 0;
  #
  # If the core gene is already in an alt_allele set then use that alt_id for the LRG gene only.
  # Else use a new one and add both core and LRG.
  #
  

  while ($sth->fetch()){
    my $aag = $aaga->fetch_by_gene_id($core_gene_id);
    if ($aag) {
        $insert_sth->execute($aag->dbID, $lrg_gene_id, 0);
        $old_count++;
    } else {
        $aag = $aaga->fetch_by_gene_id($lrg_gene_id);
        if ($aag) {
            $insert_sth->execute($aag->dbID, $lrg_gene_id, 1);
            print "LRG perculiarity\t$core_gene_id\t$lrg_gene_id\n";
            $lrg_count++;
        } else {
            $max_alt_id++;
            $insert_sth->execute($max_alt_id, $lrg_gene_id, 0);
            $insert_sth->execute($max_alt_id, $core_gene_id, 1);
            $new_count++;
        }
    }
    $count++;
  }

  
  if($count){
    print "Added $count alt_allels for the lrgs. $old_count added to previous alt_alleles and $new_count new ones\n";
    print "LRG problem count = $lrg_count\n";
  }


  $self->update_process_status("alt_alleles_added");
  return;
  
}


#
# Default behaviour is not to do the offical naming
# Overload this method in the species file returning the
# official database name to do so. 
# (ie, human-> HGNC, mouse ->MGI, zebrafisf -> ZFIN_ID)
#
sub get_official_name {
  return;
}



#
# Biomart insists that a source is linked to only one ensembl
# object type (Gene, Transcript, Translation). So biomart_fix
# will move $dbnmae entry for type1 to type 2
# i.e. move all HGNC from transcripts to Genes.
#
sub biomart_fix{
  my ($self, $db_name, $type1, $type2, $verbose, $xref_dbc) = @_;
  $xref_dbc = $self->xref->dbc unless defined $xref_dbc;

  print "$db_name is associated with both $type1 and $type2 object types\n" if(defined($verbose));
  print "$db_name moved to Gene level.\n" if(!defined($verbose));

  my $to;
  my $from;
  my $to_id;
  my $from_id;
  if($type1 eq "Gene" or $type2 eq "Gene"){
    $to = "Gene";
    $to_id = "gene_id";
    if($type1 eq "Translation" or $type2 eq "Translation"){
      $from = "Translation";
      $from_id = "translation_id"
    }
    else{
      $from = "Transcript";
      $from_id = "transcript_id";
    }
  }
  else{
    $to = "Transcript";
    $to_id = "transcript_id";
    $from = "Translation";
    $from_id = "translation_id";
  }

  if ($db_name eq 'GO' || $db_name eq 'goslim_goa') { 
    $to = 'Translation';
    $from = 'Transcript';
    $to_id = 'translation_id';
    $from_id = 'transcript_id';
  }
  
  print "Therefore moving all associations from $from to ".$to."\n" if(defined($verbose));
  

  my $sql =(<<"EOF");
  UPDATE IGNORE object_xref, gene_transcript_translation, xref, source
    SET object_xref.ensembl_object_type = "$to",
      object_xref.ensembl_id = gene_transcript_translation.$to_id 
	WHERE object_xref.ensembl_object_type = "$from" AND
	  object_xref.ensembl_id = gene_transcript_translation.$from_id AND
	    xref.xref_id = object_xref.xref_id AND
	      xref.source_id = source.source_id AND
                object_xref.ox_status = "DUMP_OUT"  AND
		  source.name = "$db_name";
EOF
  my $result =  $xref_dbc->do($sql) ;

  if($db_name eq "GO" || $db_name eq 'goslim_goa'){
    $sql =(<<"EOF2");
  DELETE object_xref, identity_xref, go_xref
    FROM object_xref, xref, source, identity_xref, go_xref
      WHERE object_xref.ensembl_object_type = "$from" AND
        identity_xref.object_xref_id = object_xref.object_xref_id AND
	xref.xref_id = object_xref.xref_id AND
          go_xref.object_xref_id = object_xref.object_xref_id AND
	  xref.source_id = source.source_id AND
            object_xref.ox_status = "DUMP_OUT"  AND
	      source.name = "$db_name";
EOF2
    
  $result = $xref_dbc->do($sql);  

# Special tidying up for transcripts without translation
# The resulting object_xref does not have an ensembl_id to map to

    $sql=(<<"EOF4");
  DELETE object_xref, identity_xref, go_xref
    FROM object_xref, xref, source, identity_xref, go_xref
      WHERE object_xref.ensembl_object_type = "$to" AND
        identity_xref.object_xref_id = object_xref.object_xref_id AND
        xref.xref_id = object_xref.xref_id AND
          go_xref.object_xref_id = object_xref.object_xref_id AND
          xref.source_id = source.source_id AND
            object_xref.ensembl_id = 0 AND
              object_xref.ox_status = "DUMP_OUT"  AND
                source.name = "$db_name";
EOF4
  }
  else{
    $sql =(<<"EOF3");
  DELETE object_xref, identity_xref
    FROM xref, source, object_xref
      LEFT JOIN identity_xref
        ON identity_xref.object_xref_id = object_xref.object_xref_id
      WHERE object_xref.ensembl_object_type = "$from" AND
	xref.xref_id = object_xref.xref_id AND
	  xref.source_id = source.source_id AND
            object_xref.ox_status = "DUMP_OUT"  AND
	      source.name = "$db_name";
EOF3

  $result = $xref_dbc->do($sql);
  }
#  print "\n$sql\n";

  #delete dependent_xref 
  $sql =(<<'EOF4');
  DELETE FROM dependent_xref WHERE object_xref_id NOT IN 
   (SELECT object_xref_id FROM object_xref);
EOF4
  return;
}


#
# This sub finds which source lie on multiple ensembl obejct types
# and calls biomart_fix to fix this.
#
sub biomart_testing{
  my ($self) = @_;

  my $sql = 'SELECT ox.ensembl_object_type, COUNT(*), s.name  FROM xref x, object_xref ox, source s  WHERE x.xref_id = ox.xref_id AND s.source_id = x.source_id  and ox.ox_status = "DUMP_OUT" GROUP BY s.name, ox.ensembl_object_type';


  my $again = 1;
  while ($again){
    $again = 0;

    my $sth = $self->xref->dbc->prepare($sql);
    $sth->execute();
    my ($type, $count, $name);
    my ($last_type, $last_count, $last_name);
    $sth->bind_columns(\$type,\$count,\$name);
    $last_name = "DEFAULT";
    while ((!$again) and $sth->fetch){
      if($last_name eq $name){
	$again  = 1;
	$self->biomart_fix($name,$last_type, $type, 1);
      }
      $last_name = $name;
      $last_type= $type;
      $last_count = $count;
    }
    $sth->finish;  
  }

  $self->update_process_status('biomart_test_finished');
  return;
}

#
# Similar to above but just reports the problems.
# It does not fix them
#

sub biomart_test{
  my ($self) = @_;

  my $sql = 'SELECT ox.ensembl_object_type, COUNT(*), s.name  FROM xref x, object_xref ox, source s  WHERE x.xref_id = ox.xref_id AND s.source_id = x.source_id  and ox.ox_status = "DUMP_OUT" GROUP BY s.name, ox.ensembl_object_type';


  my $sth = $self->xref->dbc->prepare($sql);
  
  $sth->execute();
  my ($type, $count, $name);
  my ($last_type, $last_count, $last_name);
  $sth->bind_columns(\$type,\$count,\$name);
  $last_name = "NOTKNOWN";
  my $first = 1;
  while ($sth->fetch){
    if($last_name eq $name){
      if($first){
	print STDERR "\nProblem Biomart test fails\n";
	$first=0;
      }
      print STDERR "$last_name\t$last_count\t$last_type\n";
      print STDERR "$name\t$count\t$type\n";
    }
    $last_name = $name;
    $last_type= $type;
    $last_count = $count;
  }
  $sth->finish;
  return;
}

# remove a list of patterns from a string
sub filter_by_regexp {

  my ($self, $str, $regexps) = @_;

  foreach my $regexp (@$regexps) {
    $str =~ s/$regexp//ig;
  }

  return $str;

}


sub get_species_id_from_species_name{
  my ($self,$species) = @_;


  my $sql = "select species_id from species where name = '".$species."'";
  my $sth = $self->dbc->prepare($sql);
  $sth->execute();
  my @row = $sth->fetchrow_array();
  my $species_id;
  if (@row) {
    $species_id = $row[0];
  } else {
    print STDERR "Couldn't get ID for species ".$species."\n";
    print STDERR "It must be one of :-\n";
    $sql = "select name from species";
    $sth = $self->dbc->prepare($sql);
    $sth->execute();
    while(my @row2 = $sth->fetchrow_array()){
      print STDERR $row2[0]."\n";
    }
  }
  $sth->finish();

  return $species_id;
}


sub clean_up{
  my $self = shift;
  my $stats = shift;
  my $keep_core_data = shift;

  # remove all object_xref, identity_xref  entries

  my $sql = "TRUNCATE table object_xref";
  my $sth = $self->xref->dbc->prepare($sql);
  $sth->execute(); 

  $sql = "TRUNCATE table go_xref";
  $sth = $self->xref->dbc->prepare($sql);
  $sth->execute(); 

  $sql = "TRUNCATE table identity_xref";
  $sth = $self->xref->dbc->prepare($sql);
  $sth->execute(); 
 
  # remove all xrefs after PARSED_xref_id
  # set dumped to NULL fro all xrefs.

  my $max_xref_id = $self->get_meta_value("PARSED_xref_id");

  if($max_xref_id){
    $sql = "DELETE from xref where xref_id > $max_xref_id";
    $sth = $self->xref->dbc->prepare($sql);
    $sth->execute(); 
  }

  $sql = "UPDATE xref set dumped = null";
  $sth = $self->xref->dbc->prepare($sql);
  $sth->execute(); 

  $sql = "DELETE from display_xref_priority";
  $sth = $self->xref->dbc->prepare($sql);
  $sth->execute();
      

  $sql = "DELETE from gene_desc_priority";
  $sth = $self->xref->dbc->prepare($sql);
  $sth->execute();


  if (!$keep_core_data) {
      # remove all from core_info tables
      #        gene_transcript_translation
      #        [gene/transcript/translation]_stable_id
      #
      $sql = "DELETE from gene_transcript_translation";
      $sth = $self->xref->dbc->prepare($sql);
      $sth->execute(); 

      $sql = "DELETE from gene_stable_id";
      $sth = $self->xref->dbc->prepare($sql);
      $sth->execute(); 
 
      $sql = "DELETE from transcript_stable_id";
      $sth = $self->xref->dbc->prepare($sql);
      $sth->execute(); 
 
      $sql = "DELETE from translation_stable_id";
      $sth = $self->xref->dbc->prepare($sql);
      $sth->execute(); 
  }
  return;
}

sub remove_mapping_data{
  my $self = shift;

  my $sql = "DELETE from mapping_jobs";
  my $sth = $self->xref->dbc->prepare($sql);
  $sth->execute(); 

  $sql = "DELETE from mapping";
  $sth = $self->xref->dbc->prepare($sql);
  $sth->execute();

  $sql = "DELETE from alt_allele";
  $sth = $self->xref->dbc->prepare($sql);
  $sth->execute();

  $sql = "DELETE from source_mapping_method";
  $sth = $self->xref->dbc->prepare($sql);
  $sth->execute();

  return;
}


sub revert_to_parsing_finished{
  my $self = shift;


  $self->clean_up();
  $self->remove_mapping_data();

  $self->update_process_status('parsing_finished');

  return;
}


sub revert_to_mapping_finished{
  my $self = shift;

  
  $self->clean_up(undef,1);

  # set mapping jobs to SUBMITTED
  my $sql = 'UPDATE mapping_jobs set status = "SUBMITTED"';;
  my $sth = $self->xref->dbc->prepare($sql);
  $sth->execute(); 

  $self->update_process_status('mapping_finished');
  return;
}

#
# In case we have alt alleles with xefs, these will be direct ones
# we need to move all xrefs on to the reference
#

sub get_alt_allele_hashes{
  my $self= shift;

  my %alt_to_ref;
  my %ref_to_alts;

  my $sql = "select alt_allele_id, gene_id, is_reference from alt_allele order by alt_allele_id, is_reference DESC";

  my $sth = $self->xref->dbc->prepare($sql);
  $sth->execute();
  my ($alt_allele_id,$gene_id, $is_ref);
  $sth->bind_columns(\$alt_allele_id, \$gene_id, \$is_ref);
  my $last_alt_allele = 0;
  my $ref_gene;
  while($sth->fetch()){
      if( $alt_allele_id != $last_alt_allele) {
	  #use the first non-reference gene if there is no reference gene in an alt_allele
	  $ref_gene = $gene_id;
      } else{
	  $alt_to_ref{$gene_id} = $ref_gene;
	  push @{$ref_to_alts{$ref_gene}}, $gene_id;
      }
      $last_alt_allele = $alt_allele_id;
  }
  $sth->finish;

  return \%alt_to_ref, \%ref_to_alts;
}


sub process_alt_alleles{
  my $self = shift;
  my $dbc = shift;
  $dbc = $self->xref->dbc unless defined $dbc;

  # ALL are on the Gene level now. This may change but for now it is okay.
  my ($alt_to_ref, $ref_to_alts) = $self->get_alt_allele_hashes();

  my $tester = XrefMapper::TestMappings->new($self);
  if($tester->unlinked_entries){
    croak "Problems found before process_alt_alleles\n";
  }
  #
  # Move the xrefs on to the reference Gene.
  # NOTE: Igonore used as the xref might already be on this Gene already and we do not want it to crash
  #
  my $move_sql =(<<'MOVE');
UPDATE IGNORE object_xref ox, xref x, source s 
  SET ox.ensembl_id = ? 
    WHERE x.source_id = s.source_id AND 
          ox.xref_id = x.xref_id AND
          ox.ensembl_id = ? AND
          ox.ensembl_object_type = 'Gene' AND
          ox.ox_status = 'DUMP_OUT' AND 
          s.name in (
MOVE
$move_sql .= "'".join("', '",$self->get_gene_specific_list()) . "')";

print "MOVE SQL\n$move_sql\n";

  #
  # Now where it was already on the Gene the ignore will have stopped the move
  # so we now want to just remove those ones as they already exist.
  #
  my $del_ix_sql =(<<'DIX');
DELETE ix 
  FROM identity_xref ix, object_xref ox, xref x, source s 
    WHERE x.source_id = s.source_id AND
          ox.object_xref_id = ix.object_xref_id AND
          ox.xref_id = x.xref_id AND 
          ox.ensembl_id = ? AND 
          ox.ensembl_object_type = 'Gene' AND 
          ox.ox_status = 'DUMP_OUT' AND
           s.name in (
DIX
$del_ix_sql .= "'".join("', '",$self->get_gene_specific_list()) . "')";

  my $del_sql =(<<'DEL');
DELETE ox 
  FROM object_xref ox, xref x, source s 
    WHERE x.source_id = s.source_id AND
          ox.xref_id = x.xref_id AND 
          ox.ensembl_id = ? AND 
          ox.ensembl_object_type = 'Gene' AND 
          ox.ox_status = 'DUMP_OUT' AND
           s.name in (
DEL
$del_sql .= "'".join("', '",$self->get_gene_specific_list()) . "')";

  my $move_sth = $dbc->prepare($move_sql)  || croak "$move_sql cannot be prepared";
  my $del_ix_sth = $dbc->prepare($del_ix_sql)    || croak "$del_ix_sql cannot be prepared";
  my $del_sth = $dbc->prepare($del_sql)    || croak "$del_sql cannot be prepared";

  my $move_count = 0;
  my $del_ix_count = 0;
  my $del_ox_count = 0;
  foreach my $key (keys %$alt_to_ref){
    $move_sth->execute($alt_to_ref->{$key}, $key);
    $move_count += $move_sth->rows;

    $del_ix_sth->execute($key);
    $del_ix_count += $del_ix_sth->rows;

    $del_sth->execute($key);
    $del_ox_count += $del_sth->rows;
  }
  $move_sth->finish;
  $del_sth->finish;
  $del_ix_sth->finish;

  print "Number of rows:- moved = $move_count, identitys deleted = $del_ix_count, object_xrefs deleted = $del_ox_count\n";
  if($tester->unlinked_entries){
    croak "Problems found mid process_alt_alleles\n";
  }
  #
  # Now we have all the data on the reference Gene we want to copy all the data
  # onto the alt alleles.
  #


  my $get_data_sql=(<<'GET');
SELECT ox.object_xref_id, ox.ensembl_object_type, ox.xref_id, ox.linkage_annotation, 
       ox.linkage_type, ox.ox_status, ox.unused_priority, ox.master_xref_id,
       ix.query_identity, ix.target_identity, ix.hit_start, ix.hit_end,
       ix.translation_start, ix.translation_end, ix.cigar_line, ix.score, ix.evalue
  FROM xref x, source s, object_xref ox 
    LEFT JOIN identity_xref ix ON ox.object_xref_id =ix.object_xref_id 
      WHERE  x.source_id = s.source_id AND
             ox.xref_id = x.xref_id AND
             ox.ensembl_id = ? AND
             ox.ox_status = 'DUMP_OUT' AND
             ox.ensembl_object_type = 'Gene' AND
              s.name in (
GET

  $get_data_sql .= "'".join("', '",$self->get_gene_specific_list()) . "')";

  my $get_data_sth = $self->xref->dbc->prepare($get_data_sql) || croak "Could not prepare $get_data_sql";



  my $insert_object_xref_sql =(<<'INO');
INSERT INTO object_xref (object_xref_id, ensembl_id, ensembl_object_type, xref_id, linkage_annotation, 
            linkage_type, ox_status, unused_priority, master_xref_id) 
       VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
INO

  my $insert_ox_sth = $self->xref->dbc->prepare($insert_object_xref_sql) || croak "Could not prepare $insert_object_xref_sql";


  my $insert_identity_xref_sql = (<<'INI');
INSERT INTO identity_xref (object_xref_id, query_identity, target_identity, hit_start, hit_end,
            translation_start, translation_end, cigar_line, score, evalue ) 
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
INI

  my $insert_ix_sth = $self->xref->dbc->prepare($insert_identity_xref_sql) || croak "Could not prepare $insert_identity_xref_sql";



  my $max_object_xref_id;

  my $sth = $self->xref->dbc->prepare("SELECT MAX(object_xref_id) FROM object_xref");
  $sth->execute();
  $sth->bind_columns(\$max_object_xref_id);
  $sth->fetch;
  if((!defined($max_object_xref_id)) or (!$max_object_xref_id)){
    croak "Problem getting max object_xref_id";
  }
  $max_object_xref_id++;

  my $added_count = 0;
  my $ignored = 0;
  foreach my $key (keys %$ref_to_alts){
    $get_data_sth->execute($key);
    my ($object_xref_id, $ensembl_object_type, $xref_id, $linkage_annotation,
	$linkage_type, $ox_status, $unused_priority, $master_xref_id,
	$query_identity, $target_identity, $hit_start, $hit_end,
	$translation_start, $translation_end, $cigar_line, $score, $evalue);

    $get_data_sth->bind_columns(\$object_xref_id, \$ensembl_object_type, \$xref_id, \$linkage_annotation,
				\$linkage_type, \$ox_status, \$unused_priority, \$master_xref_id,
				\$query_identity, \$target_identity, \$hit_start, \$hit_end,
				\$translation_start, \$translation_end, \$cigar_line, \$score, \$evalue);

    while( $get_data_sth->fetch()){
      foreach my $alt (@{$ref_to_alts->{$key}}){
	$max_object_xref_id++;
        $insert_ox_sth->execute($max_object_xref_id, $alt, $ensembl_object_type, $xref_id, $linkage_annotation,
				$linkage_type, $ox_status, $unused_priority, $master_xref_id) || croak "Could not insert object_xref data";

#ONLY add identity xref if object_xref was added successfully.
	if( $insert_ox_sth->rows){
	  $added_count++;
	  $insert_ix_sth->execute($max_object_xref_id, $query_identity, $target_identity, $hit_start, $hit_end,
				$translation_start, $translation_end, $cigar_line, $score, $evalue) ||  croak "Could not insert identity_xref data";
	}
	else{
	  $ignored++;
	}
      }
    }
  }
  print "Added $added_count new mapping but ignored $ignored\n";
  
  if($tester->unlinked_entries){
    croak "Problems found after process_alt_alleles\n";
  }

  $self->update_process_status('alt_alleles_processed');
  return;
}


#
# These sources should be on the gene, even if they are mapped transcript or translation.
# We define which ones are to be moved here
#
sub get_gene_specific_list {
  my $self = shift;

  my @list = qw(DBASS3 DBASS5 EntrezGene miRBase RFAM TRNASCAN_SE RNAMMER UniGene Uniprot_gn WikiGene MIM_GENE MIM_MORBID HGNC);

  return @list;
}


#
# Here we do the moving.
#
sub source_defined_move{
  my $self = shift;
  my $dbi = shift;

  my $tester = XrefMapper::TestMappings->new($self);
  if($tester->unlinked_entries){
    croak "Problems found before source_defined_move\n";
  }
  foreach my $source ($self->get_gene_specific_list()){
    $self->biomart_fix($source,"Translation","Gene", undef, undef, $dbi);
    $self->biomart_fix($source,"Transcript","Gene", undef, undef, $dbi);
  }
  if($tester->unlinked_entries){
    croak "Problems found after source_defined_move\n";
  }
  $self->update_process_status('source_level_move_finished');
  return;
}

1;
