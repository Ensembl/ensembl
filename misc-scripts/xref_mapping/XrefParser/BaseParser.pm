package XrefParser::BaseParser;

use strict;

use DBI;
use Digest::MD5 qw(md5_hex);
use File::Path;
use File::Basename;
use POSIX qw(strftime);
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception;

my $base_dir = ".";

my $add_xref_sth = undef;
my $add_direct_xref_sth = undef;
my $add_dependent_xref_sth = undef;
my $get_xref_sth = undef;
my $add_synonym_sth = undef;

my $dbi;
my %dependent_sources;
my %taxonomy2species_id;
my %name2species_id;

my ($host, $port, $dbname, $user, $pass, $create, $release, $cleanup, $deletedownloaded);
my ($skipdownload,$drop_db,$checkdownload) ;

# --------------------------------------------------------------------------------
# Get info about files to be parsed from the database

sub run {

  ($host, $port, $dbname, $user, $pass, my $speciesr, my $sourcesr, $skipdownload, $checkdownload, 
    $create, $release, $cleanup, $drop_db, $deletedownloaded) = @_;



  my @species = @$speciesr;
  my @sources = @$sourcesr;

  my $sql_dir = dirname($0);

  create($host, $port, $user, $pass, $dbname, $sql_dir."/" , $drop_db ) if ($create);

  my $dbi = dbi();

  # validate species names
  my @species_ids = validate_species(@species);

  # validate source names
  exit(1) if (!validate_sources(@sources));

  # build SQL
  my $species_sql = "";
  if (@species_ids) {
    $species_sql .= " AND su.species_id IN (";
    for (my $i = 0; $i < @species_ids; $i++ ) {
      $species_sql .= "," if ($i ne 0);
      $species_sql .= $species_ids[$i];
    }
    $species_sql .= ") ";
  }

  my $source_sql = "";
  if (@sources) {
    $source_sql .= " AND LOWER(s.name) IN (";
    for (my $i = 0; $i < @sources; $i++ ) {
      $source_sql .= "," if ($i ne 0);
      $source_sql .= "\'" . lc($sources[$i]) . "\'";
    }
    $source_sql .= ") ";
  }

  my $sql =
    "SELECT s.source_id, su.source_url_id, s.name, su.url, su.checksum, su.parser, su.species_id " .
      "FROM source s, source_url su " .
	"WHERE s.download='Y' AND su.source_id=s.source_id " .
	  $source_sql . $species_sql .
	  "ORDER BY s.ordered";
#  print $sql . "\n";

  my $sth = $dbi->prepare($sql);
  $sth->execute();
  my ($source_id, $source_url_id, $name, $url, $checksum, $parser, $species_id);
  $sth->bind_columns(\$source_id, \$source_url_id, \$name, \$url, \$checksum, \$parser, \$species_id);
  my $last_type = "";
  my $dir;
  my %summary=();
  while (my @row = $sth->fetchrow_array()) {

    # Download each source into the appropriate directory for parsing later
    # or call the appropriate database parser if appropriate
    # Also delete previous working directory if we're starting a new source type

    # can have more than one file

    my @files = split(/\s+/,$url);

    my $parse = 0;
    my $empty = 0;
    my $file_cs=0;
    my $type = $name;
    my @new_file=();
    $dir = $base_dir . "/" . sanitise($type);
    my $dsn;
    $summary{$parser} = 0;
    ##
    # for summary purposes if 0 is returned then it is successful.
    #                         1 is returned then if failed.
    #                     undef/nothing  is returned the we do not know
    ##
    foreach my $urls (@files){

      # Database parsing
      if ($urls =~ /^mysql:/i) {
	$dsn = $urls;
	print "Parsing $dsn with $parser\n";
        eval "require XrefParser::$parser";
        my $new = "XrefParser::$parser"->new();
        if($new->run($dsn, $source_id, $species_id)){
	  $summary{$parser}++;
	}
	next;
      }

      # Local files need to be dealt with specially; assume they are specified as
      # LOCAL:location/of/file/relative/to/xref_mapper
      my ($file) = $urls =~ /.*\/(.*)/;
      if ($urls =~ /^LOCAL:(.*)/i) {
	my $local_file = $1;
	$file_cs = md5sum("$dir/$file");
	if(defined($file_cs)){
	  if (!defined $checksum || $checksum ne $file_cs) {
	    print "Checksum for $file does not match, parsing\n";
	    print "Parsing local file $local_file with $parser\n";
	    eval "require XrefParser::$parser";
	    my $new = "XrefParser::$parser"->new();
	    if($new->run($local_file, $source_id, $species_id)){
	      $summary{$parser}++;
	    }
	    else{
	      update_source($dbi, $source_url_id, $file_cs, $file);
	    }
	  }
	  else{
	    print "Ignoring $file as checksums match\n";
	  }
	}
	else{
	  $summary{$parser}++;
	}
	next;
      }
      # Download files

      if ( $checkdownload ) {   
         
         my $check_file = "$dir/$file" ;   
         $check_file =~s/\.gz//; 
         $check_file =~s/\.Z//; 
         print "checking for file $check_file\n" ; 
         if (-e "$check_file" ) { 
           print "SKIPPING $file because file $check_file already exists\n" ; 
           $skipdownload = 1 ;
         } else { 
          print "File $check_file does not exist. \n scheduling $dir/$file for download ..........\n" ; 
        }
      }


      $file =~ s/[&=]//g;
      if (length($file) > 100) {
	#	  $file = time;
	$file = md5_hex($file);
	print"URL is longer than 100 charcters; renamed to $file\n";
      }

      # File parsing
      if (!$skipdownload) {

	rmtree $dir if ($type ne $last_type && $deletedownloaded);
	mkdir $dir if (!-e $dir);

	$last_type = $type;


	print "Downloading $urls to $dir/$file\n";

	my $num_attempts = 0;
	my $missing = 1;
	while($num_attempts < 5 and $missing){
	  my $result = system("wget", "--quiet","--directory-prefix=$dir", "--output-document=$dir/$file",  $urls );
	  
	  # check that the file actually downloaded; may not (e.g. if too many anonymous users)
	  if ($result != 0) {
	    print "wget returned exit code $result; $type file $file not downloaded.\n";
	    print "waiting for 3 minutes then trying again\n";
	    sleep(180);
	  }
	  else{
	    $missing=0;
	  }
	  $num_attempts++;
	}
	if($missing){
	  die "Could not get $type file $file tried 5 times but failed\n";
	}
	# if the file is compressed, the FTP server may or may not have automatically uncompressed it
	# TODO - read .gz file directly? open (FILE, "zcat $file|") or Compress::Zlib
	if ($file =~ /(.*)\.gz$/ or $file =~ /(.*)\.Z$/) {
	  print "Uncompressing $dir/$file\n";
	  system("gunzip -f $dir/$file");
	  $file = $1;
	}

      }
      else{
	if ($file =~ /(.*)\.gz$/ or $file =~ /(.*)\.Z$/) {
	  $file = $1;
	}
      }

      push @new_file, $file;

      # compare checksums and parse/upload if necessary
      # need to check file size as some .SPC files can be of zero length
      $file_cs = md5sum("$dir/$file");
      if(defined($file_cs)){
	if (!defined $checksum || $checksum ne $file_cs) {
	  if (-s "$dir/$file") {
	    $parse =1;
	    print "Checksum for $file does not match, parsing\n";
	    
	    # Files from sources "Uniprot/SWISSPROT" and "Uniprot/SPTREMBL" are
	    # all parsed with the same parser
	    $parser = 'UniProtParser' if ($parser eq "Uniprot/SWISSPROT" || $parser eq "Uniprot/SPTREMBL");
	  }
	  else {
	    $empty = 1;
	    print $file . " has zero length, skipping\n";
	  }
	}
      }
      else{
	$summary{$parser}++;	
      }
    }
    
    if($parse and defined($new_file[0]) and defined($file_cs)){

      print "Parsing ".join(' ',@new_file)." with $parser\n";
      eval "require XrefParser::$parser";
      my $new = "XrefParser::$parser"->new();
      if($new->run("$dir/$new_file[0]", $source_id, $species_id)){
	$summary{$parser}++;
      }
      
      # update AFTER processing in case of crash.
      update_source($dbi, $source_url_id, $file_cs, $new_file[0]);
      
      # set release if specified
      set_release($release, $source_id) if ($release);
      
      unlink("$dir/$new_file[0]") if ($cleanup);
      
    }
    elsif(!$dsn && !$empty && defined($new_file[0])){
      print "Ignoring ".join(' ',@new_file)." as checksums match\n";
    }

  }
  print "---------------------------------------------------------\n";
  print "Summary of status\n";
  print "---------------------------------------------------------\n";
  foreach my $key (keys %summary){
    print $key."\t";
    if($summary{$key}){
      print "FAILED\n";
    }
    else{
      print "OKAY\n";
    }
  }



  # remove last working directory
  # TODO reinstate after debugging
  #rmtree $dir;

}

# --------------------------------------------------------------------------------

sub new {

  my $self = {};
  bless $self, "BaseParser";

  return $self;

}

# --------------------------------------------------------------------------------
# Get source ID for a particular file; matches url field

sub get_source_id_for_filename {

  my ($self, $file) = @_;
  print STDERR "FILE $file\n" ; 
  my $sql = "SELECT s.source_id FROM source s, source_url su WHERE su.source_id=s.source_id AND su.url LIKE  '%/" . $file . "%'";
  my $sth = dbi()->prepare($sql);
  $sth->execute();
  my @row = $sth->fetchrow_array();
  my $source_id;
  if (@row) {
    $source_id = $row[0];
  } 
  else {
    if($file =~ /rna.fna/ or $file =~ /gpff/){
      $source_id = 3;
    }else{ 
      warn("Couldn't get source ID for file $file\n");
      $source_id = -1;
    }
  }
  

  return $source_id;

}

sub rename_url_file{
  return undef;
}

# Get species ID for a particular file; matches url field

sub get_species_id_for_filename {

  my ($self, $file) = @_;

  my $sql = "SELECT su.species_id FROM source_url su WHERE su.url LIKE  '%/" . $file . "%'";
  my $sth = dbi()->prepare($sql);
  $sth->execute();
  my @row = $sth->fetchrow_array();
  my $source_id;
  if (@row) {
    $source_id = $row[0];
  } else {
    warn("Couldn't get species ID for file $file\n");
    $source_id = -1;
  }

  return $source_id;

}

# --------------------------------------------------------------------------------
# Get source ID for a particular source name

sub get_source_id_for_source_name {
  
  my ($self, $source_name) = @_;
  my $sql = "SELECT source_id FROM source WHERE name='" . $source_name . "'";
  my $sth = dbi()->prepare($sql);
  $sth->execute();
  my @row = $sth->fetchrow_array();
  my $source_id;
  if (@row) {
    $source_id = $row[0]; 
  } else {
    print STDERR "WARNING: There is no entity $source_name in the source-table of the \n" .
      "WARNING: xref-database. The name of the external db name ($source_name) is hardcoded\n" .
	"WARNING: in the parser\n";
    warn("WARNING: Couldn't get source ID for source name $source_name\n");

    $source_id = -1;
  }
  return $source_id;
}

sub get_source_name_for_source_id {
  my ($self, $source_id) = @_;
  my $source_name;
  
  my $sql = "SELECT name FROM source WHERE source_id= '" . $source_id. "'";
  my $sth = dbi()->prepare($sql);
  $sth->execute();
  my @row = $sth->fetchrow_array();
  if (@row) {
    $source_name = $row[0]; 
  } else {
    print STDERR "WARNING: There is no entity with source-id  $source_id  in the source-table of the \n" .
      "WARNING: xref-database. The source-id and the name of the source-id is hard-coded in populate_metadata.sql\n" .
	"WARNING: and in the parser\n";
    warn("WARNING: Couldn't get source name for source ID $source_id\n");
    $source_name = -1;
  }
  return $source_name;
}









sub get_valid_xrefs_for_dependencies{
  my ($self, $dependent_name, @reverse_ordered_source_list) = @_;

  my %dependent_2_xref;


  my $sql = "select source_id from source where name =?";
  my $sth = dbi()->prepare($sql);
  my @dependent_sources;
  $sth->execute($dependent_name);
  while(my @row = $sth->fetchrow_array()){
    push @dependent_sources,$row[0];
  }

  my @sources;
  foreach my $name (@reverse_ordered_source_list){
    $sth->execute($name);
    while(my @row = $sth->fetchrow_array()){
      push @sources,$row[0];
    }
  }
  $sth->finish;

  $sql  = "select d.master_xref_id, x2.accession ";
  $sql .= "  from dependent_xref d, xref x1, xref x2 ";
  $sql .= "    where x1.xref_id = d.master_xref_id and";
  $sql .= "          x1.source_id=? and ";
  $sql .= "          x2.xref_id = d.dependent_xref_id and";
  $sql .= "          x2.source_id=? ";
  
  my $sth = dbi()->prepare($sql);
  foreach my $d (@dependent_sources){
    foreach my $s (@sources){
       $sth->execute($s,$d);
       while(my @row = $sth->fetchrow_array()){
	 $dependent_2_xref{$row[1]} = $row[0];
       }
     }
  }

  return \%dependent_2_xref;
}

sub get_valid_xrefs_for_direct_xrefs{
  my ($self, $direct_name, @list) = @_;

  my %direct_2_xref;


  my $sql = "select source_id from source where name like ?";
  my $sth = dbi()->prepare($sql);
  my @direct_sources;
  $sth->execute($direct_name."%");
  while(my @row = $sth->fetchrow_array()){
    push @direct_sources,$row[0];
  }

  my @sources;
  foreach my $name (@list){
    $sth->execute($name);
    while(my @row = $sth->fetchrow_array()){
      push @sources,$row[0];
    }
  }
  $sth->finish;

  $sql  = "select d.general_xref_id, d.ensembl_stable_id, d.type, d.linkage_xref, x1.accession ";
  $sql .= "  from direct_xref d, xref x1 ";
  $sql .= "    where x1.xref_id = d.general_xref_id and";
  $sql .= "          x1.source_id=?";
   
  my $sth = dbi()->prepare($sql);
  foreach my $d (@direct_sources){
    $sth->execute($d);
    while(my @row = $sth->fetchrow_array()){
      $direct_2_xref{$row[4]} = $row[0]."::".$row[1]."::".$row[2]."::".$row[3];
    }
  }

  return \%direct_2_xref;
}



sub get_valid_codes{

  my ($self,$source_name,$species_id) =@_;

  # First cache synonyms so we can quickly add them later
  my %synonyms;
  my $syn_sth = dbi()->prepare("SELECT xref_id, synonym FROM synonym");
  $syn_sth->execute();

  my ($xref_id, $synonym);
  $syn_sth->bind_columns(\$xref_id, \$synonym);
  while ($syn_sth->fetch()) {

    push @{$synonyms{$xref_id}}, $synonym;

  }

  my %valid_codes;
  my @sources;

  my $sql = "select source_id from source where upper(name) like '%".uc($source_name)."%'";
  my $sth = dbi()->prepare($sql);
  $sth->execute();
  while(my @row = $sth->fetchrow_array()){
    push @sources,$row[0];
  }
  $sth->finish;

  foreach my $source (@sources){
    $sql = "select accession, xref_id from xref where species_id = $species_id and source_id = $source";
    my $sth = dbi()->prepare($sql);
    $sth->execute();
    while(my @row = $sth->fetchrow_array()){
      $valid_codes{$row[0]} =$row[1];
      # add any synonyms for this xref as well
      foreach my $syn (@{$synonyms{$row[1]}}) {
	$valid_codes{$syn} = $row[1];
      }
    }
  }
  return \%valid_codes;
}

# --------------------------------------------------------------------------------
# Upload xrefs to the database

sub upload_xref_object_graphs {
  my ($self, $rxrefs) = @_;

  my $dbi = dbi();
  print "count = ".$#$rxrefs."\n";

  if ($#$rxrefs > -1) {

    # remove all existing xrefs with same source ID(s)
#    $self->delete_by_source($rxrefs);

    # upload new ones
    print "Uploading xrefs\n";
    my $xref_sth = $dbi->prepare("INSERT INTO xref (accession,version,label,description,source_id,species_id) VALUES(?,?,?,?,?,?)");
    my $pri_insert_sth = $dbi->prepare("INSERT INTO primary_xref VALUES(?,?,?,?)");
    my $pri_update_sth = $dbi->prepare("UPDATE primary_xref SET sequence=? WHERE xref_id=?");
    my $syn_sth = $dbi->prepare("INSERT INTO synonym VALUES(?,?)");
    my $dep_sth = $dbi->prepare("INSERT INTO dependent_xref VALUES(?,?,?,?)");
    my $xref_update_label_sth = $dbi->prepare("UPDATE xref SET label=? WHERE xref_id=?");
    my $xref_update_descr_sth = $dbi->prepare("UPDATE xref SET description=? WHERE xref_id=?");
    my $pair_sth = $dbi->prepare("INSERT INTO pairs VALUES(?,?,?)");

    local $xref_sth->{RaiseError}; # disable error handling here as we'll do it ourselves
    local $xref_sth->{PrintError};

    foreach my $xref (@{$rxrefs}) {
       my $xref_id=undef;
       if(!defined($xref->{ACCESSION})){
	 print "your xref does not have an accession-number,so it can't be stored in the database\n";
	 return undef;
       }
      # Create entry in xref table and note ID
      if(! $xref_sth->execute($xref->{ACCESSION},
			 $xref->{VERSION},
			 $xref->{LABEL},
			 $xref->{DESCRIPTION},
			 $xref->{SOURCE_ID},
			 $xref->{SPECIES_ID})){
	if(!defined($xref->{SOURCE_ID})){
	  print "your xref: $xref->{ACCESSION} does not have a source-id\n";
	  return undef;
	}
	$xref_id = insert_or_select($xref_sth, $dbi->err, $xref->{ACCESSION}, $xref->{SOURCE_ID});
	$xref_update_label_sth->execute($xref->{LABEL},$xref_id) if (defined($xref->{LABEL}));
	$xref_update_descr_sth->execute($xref->{DESCRIPTION},$xref_id,) if (defined($xref->{DESCRIPTION}));
      }
      else{
	$xref_id = insert_or_select($xref_sth, $dbi->err, $xref->{ACCESSION}, $xref->{SOURCE_ID});
      }


      # create entry in primary_xref table with sequence; if this is a "cumulative"
      # entry it may already exist, and require an UPDATE rather than an INSERT
      if(!(defined($xref_id) and $xref_id)){
	print STDERR "xref_id is not set for :\n$xref->{ACCESSION}\n$xref->{LABEL}\n$xref->{DESCRIPTION}\n$xref->{SOURCE_ID}\n";
      }
      if (primary_xref_id_exists($xref_id)) {
	
	$pri_update_sth->execute($xref->{SEQUENCE}, $xref_id) || die $dbi->errstr;
	
      } else {
	
	$pri_insert_sth->execute($xref_id,
				 $xref->{SEQUENCE},
				 $xref->{SEQUENCE_TYPE},
				 $xref->{STATUS}) || die $dbi->errstr;
      }

      # if there are synonyms, add entries in the synonym table
      foreach my $syn (@{$xref->{SYNONYMS}}) {
	$syn_sth->execute($xref_id, $syn) || die "$dbi->errstr \n $xref_id\n $syn\n";
      }	# foreach syn

      # if there are dependent xrefs, add xrefs and dependent xrefs for them
      foreach my $depref (@{$xref->{DEPENDENT_XREFS}}) {

	my %dep = %$depref;

	$xref_sth->execute($dep{ACCESSION},
			   $dep{VERSION},
			   $dep{LABEL},
			   "",
			   $dep{SOURCE_ID},
			   $xref->{SPECIES_ID});

	my $dep_xref_id = insert_or_select($xref_sth, $dbi->err, $dep{ACCESSION}, $dep{SOURCE_ID});
			
	if($dbi->err){
	  print STDERR "dbi\t$dbi->err \n$dep{ACCESSION} \n $dep{SOURCE_ID} \n";
	}
	if(!defined($dep_xref_id) || $dep_xref_id ==0 ){
	  print STDERR "acc = $dep{ACCESSION} \nlink = $dep{LINKAGE_SOURCE_ID} \n".$dbi->err."\n";
	}
	$dep_sth->execute($xref_id, $dep_xref_id, $dep{LINKAGE_ANNOTATION}, $dep{LINKAGE_SOURCE_ID} ) || die $dbi->errstr;
	# TODO linkage anntation?

      }	 # foreach dep
       
       if(defined($xref_id) and defined($xref->{PAIR})){
	 $pair_sth->execute($xref->{SOURCE_ID},$xref->{ACCESSION},$xref->{PAIR});
       }				
       
              
       $xref_sth->finish() if defined $xref_sth;
       $pri_insert_sth->finish() if defined $pri_insert_sth;
       $pri_update_sth->finish() if defined $pri_update_sth;
       
     }  # foreach xref

  }
  return 1;
}

sub upload_direct_xrefs{
  my ($self, $direct_xref)  = @_;
  for my $dr(@$direct_xref) {
    # print "having now direct-XREF : $dr->{ENSEMBL_STABLE_ID} \n" ;
    my $general_xref_id = get_xref_id_by_accession_and_source($dr->{ACCESSION},$dr->{SOURCE_ID});
    if ($general_xref_id){
      # print "direct_xref:\n$general_xref_id\n$dr->{ENSEMBL_STABLE_ID}\n$dr->{ENSEMBL_TYPE}\t$dr->{LINKAGE_XREF}\n\n";
      $self->add_direct_xref($general_xref_id, $dr->{ENSEMBL_STABLE_ID},$dr->{ENSEMBL_TYPE},$dr->{LINKAGE_XREF});
    }
  }
}



# --------------------------------------------------------------------------------
# Get & cache a hash of all the source names for dependent xrefs (those that are
# in the source table but don't have an associated URL etc)

sub get_dependent_xref_sources {

  my $self = shift;

  if (!%dependent_sources) {

    my $dbi = dbi();
    my $sth = $dbi->prepare("SELECT name,source_id FROM source WHERE download='N'");
    $sth->execute() || die $dbi->errstr;
    while(my @row = $sth->fetchrow_array()) {
      my $source_name = $row[0];
      my $source_id = $row[1];
      $dependent_sources{$source_name} = $source_id;
    }
  }

  return %dependent_sources;

}

# --------------------------------------------------------------------------------
# Get & cache a hash of all the species IDs & taxonomy IDs.

sub taxonomy2species_id {

  my $self = shift;

  if (!%taxonomy2species_id) {

    my $dbi = dbi();
    my $sth = $dbi->prepare("SELECT species_id, taxonomy_id FROM species");
    $sth->execute() || die $dbi->errstr;
    while(my @row = $sth->fetchrow_array()) {
      my $species_id = $row[0];
      my $taxonomy_id = $row[1];
      $taxonomy2species_id{$taxonomy_id} = $species_id;
    }
  }

  return %taxonomy2species_id;

}

# --------------------------------------------------------------------------------
# Get & cache a hash of all the species IDs & species names.

sub name2species_id {

  my $self = shift;

  if (!%name2species_id) {

    my $dbi = dbi();
    my $sth = $dbi->prepare("SELECT species_id, name FROM species");
    $sth->execute() || die $dbi->errstr;
    while(my @row = $sth->fetchrow_array()) {
      my $species_id = $row[0];
      my $name = $row[1];
      $name2species_id{$name} = $species_id;
    }
  }

  return %name2species_id;

}

# --------------------------------------------------------------------------------
# Update a row in the source table

sub update_source {

  my ($dbi, $source_url_id, $checksum, $file) = @_;
  open(FILE, $file);
  my $file_date = POSIX::strftime('%Y%m%d%H%M%S', localtime((stat($file))[9]));
  close(FILE);

  my $sql = "UPDATE source_url SET checksum='" . $checksum . "', file_modified_date='" . $file_date . "', upload_date=NOW() WHERE source_url_id=" . $source_url_id;
  # TODO release?

  $dbi->prepare($sql)->execute() || die $dbi->errstr;

}


# --------------------------------------------------------------------------------

sub dbi {

  my $self = shift;

  if (!defined $dbi) {
    $dbi = DBI->connect("dbi:mysql:host=$host;port=$port;database=$dbname",
			"$user",
			"$pass",
			{'RaiseError' => 1}) || die "Can't connect to database";
  }

  return $dbi;

}

# --------------------------------------------------------------------------------

sub md5sum {

  my $file = shift;

  unless (-e $file) {
	print "\n\nWarning: can't find file $file - you have to download it again. \n\n" ; 
        print " SKIPPING $file\n" ; 
#        sleep(10) ; 	
	return undef; 
  }

  open(FILE, $file);
  binmode(FILE);
  my $md5 = Digest::MD5->new->addfile(*FILE)->hexdigest();
  close(FILE);

  return $md5;

}

# --------------------------------------------------------------------------------

sub get_xref_id_by_accession_and_source {

  my ($acc, $source_id) = @_;

  my $dbi = dbi();
  my $sth = $dbi->prepare("SELECT xref_id FROM xref WHERE accession=? AND source_id=?");
  $sth->execute($acc, $source_id) || die $dbi->errstr;
  my @row = $sth->fetchrow_array();
  my $xref_id = $row[0];

  return $xref_id;

}

# --------------------------------------------------------------------------------
# If there was an error, an xref with the same acc & source already exists.
# If so, find its ID, otherwise get ID of xref just inserted

sub insert_or_select {

  my ($sth, $error, $acc, $source) = @_;

  my $id;

  # TODO - check for specific error code rather than for just any error
  if ($error) {

    $id = get_xref_id_by_accession_and_source($acc, $source);
    #print STDERR "Got existing xref id " . $id . " for " . $acc . " " . $source . "\n";
	
  } else {
	
    $id = $sth->{'mysql_insertid'};
	
  }

  return $id;

}

# --------------------------------------------------------------------------------

sub primary_xref_id_exists {

  my $xref_id = shift;

  my $exists = 0;

  my $dbi = dbi();
  my $sth = $dbi->prepare("SELECT xref_id FROM primary_xref WHERE xref_id=?");
  $sth->execute($xref_id) || die $dbi->errstr;
  my @row = $sth->fetchrow_array();
  my $result = $row[0];
  $exists = 1 if (defined $result);

  return $exists;

}

# --------------------------------------------------------------------------------

# delete all xrefs & related objects

sub delete_by_source {

  my $self =shift;
  my $xrefs = shift;

  # SQL for deleting stuff
  # Note this SQL only works on MySQL version 4 and above

  #Remove direct xrefsbased on source
  my $direct_sth = $dbi->prepare("DELETE FROM direct_xref USING xref, direct_xref WHERE xref.xref_id=direct_xref.general_xref_id AND xref.source_id=?");
  
  #remove Pairs fro source
  my $pairs_sth = $dbi->prepare("DELETE FROM pairs WHERE source_id=?");

  # Remove dependent_xrefs and synonyms based on source of *xref*
  my $syn_sth = $dbi->prepare("DELETE FROM synonym USING xref, synonym WHERE xref.xref_id=synonym.xref_id AND xref.source_id=?");
  my $dep_sth = $dbi->prepare("DELETE FROM dependent_xref USING xref, dependent_xref WHERE xref.xref_id=dependent_xref.master_xref_id AND xref.source_id=?");

  # xrefs and primary_xrefs are straightforward deletes
  my $xref_sth = $dbi->prepare("DELETE FROM xref, primary_xref USING xref, primary_xref WHERE source_id=? AND primary_xref.xref_id = xref.xref_id");
#  my $p_xref_sth = $dbi->prepare("DELETE FROM primary_xref WHERE source_id=?");

  # xrefs may come from more than one source (e.g. UniProt/SP/SPtr)
  # so find all sources first
  my %source_ids;
  foreach my $xref (@$xrefs) {
    my $xref_source = $xref->{SOURCE_ID};
    $source_ids{$xref_source} = 1;
  }

  # now delete them
  foreach my $source (keys %source_ids) {
    print "Deleting pairs with source ID $source \n";
    $pairs_sth->execute($source);
    print "Deleting direct xrefs with source ID $source \n";
    $direct_sth->execute($source);
    print "Deleting synonyms of xrefs with source ID $source \n";
    $syn_sth->execute($source);
    print "Deleting dependent xrefs of xrefs with source ID $source \n";
    $dep_sth->execute($source);
    print "Deleting primary xrefs with source ID $source \n";
#    $p_xref_sth->execute($source);
    print "Deleting xrefs with source ID $source \n";
    $xref_sth->execute($source);
  }

  $syn_sth->finish() if defined $syn_sth;
  $dep_sth->finish() if defined $dep_sth;
  $xref_sth->finish() if defined $xref_sth;
#  $p_xref_sth->finish() if defined $p_xref_sth;

}

# --------------------------------------------------------------------------------

sub validate_sources {

  my @sources = @_;

  my $dbi = dbi();
  my $sth = $dbi->prepare("SELECT * FROM source WHERE LOWER(name)=?");

  foreach my $source (@sources) {

    $sth->execute(lc($source));
    if ($sth->fetchrow_array()) {
      print "Source $source is valid\n";
    } else {
      print "\nSource $source is not valid; valid sources are:\n";
      show_valid_sources();
      return 0;
    }

  }

  return 1;

}

# --------------------------------------------------------------------------------

sub show_valid_sources() {

  my $dbi = dbi();
  my $sth = $dbi->prepare("SELECT name FROM source WHERE download='Y'");

  $sth->execute();
  while (my @row = $sth->fetchrow_array()) {
    print @row[0] . "\n";
  }

}

# --------------------------------------------------------------------------------

sub validate_species {
  my @species = @_;
  my @species_ids;

  my $dbi = dbi();
  my $sth = $dbi->prepare("SELECT species_id, name FROM species WHERE LOWER(name)=? OR LOWER(aliases) LIKE ?");
  my ($species_id, $species_name);

  foreach my $sp (@species) {

    $sth->execute(lc($sp), "%" . lc($sp) . "%");
    $sth->bind_columns(\$species_id, \$species_name);
    if (my @row = $sth->fetchrow_array()) {
      print "Species $sp is valid (name = " . $species_name . ", ID = " . $species_id . ")\n";
      push @species_ids, $species_id;
    } else {
      print "Species $sp is not valid; valid species are:\n";
      show_valid_species();
      exit(1);
    }
  }
  return @species_ids;
}

# --------------------------------------------------------------------------------

sub show_valid_species() {

  my $dbi = dbi();
  my $sth = $dbi->prepare("SELECT name, aliases FROM species");

  $sth->execute();
  while (my @row = $sth->fetchrow_array()) {
    print @row[0] . " (aliases: " . $row[1] . ")\n";
  }

}

sub get_taxonomy_from_species_id{
  my ($self,$species_id) = @_;

  my $dbi = dbi();
  my $sth = $dbi->prepare("SELECT taxonomy_id FROM species WHERE species_id = $species_id");
  $sth->execute() || die $dbi->errstr;
  if(my @row = $sth->fetchrow_array()) {
    return $row[0];
  }   
  $sth->finish;
  return undef;
}
sub get_direct_xref{
  my ($self,$stable_id,$type,$link) = @_;

  my $direct_sth;
  if(!defined($direct_sth)){
    my $sql = "select general_xref_id from direct_xref d where ensembl_stable_id = ? and type = ?  and linkage_xref= ?";
    $direct_sth = $dbi->prepare($sql);  
  }
  
  $direct_sth->execute($stable_id, $type, $link) || die $dbi->errstr;
  if(my @row = $direct_sth->fetchrow_array()) {
    return $row[0];
  }   
  return undef;
}

sub get_xref{
  my ($self,$acc,$source) = @_;

  if(!defined($get_xref_sth)){
    my $sql = "select xref_id from xref where accession = ? and source_id = ?";
    $get_xref_sth = $dbi->prepare($sql);  
  }
  
  $get_xref_sth->execute($acc, $source) || die $dbi->errstr;
  if(my @row = $get_xref_sth->fetchrow_array()) {
    return $row[0];
  }   
  return undef;
}

sub add_xref {

  my ($self,$acc,$version,$label,$description,$source_id,$species_id) = @_;

  if(!defined($add_xref_sth)){
    $add_xref_sth = dbi->prepare("INSERT INTO xref (accession,version,label,description,source_id,species_id) VALUES(?,?,?,?,?,?)");
  }
 $add_xref_sth->execute($acc,$version,$label,$description,$source_id,$species_id) 
   || die "$acc\t$label\t\t$source_id\t$species_id\n";

  return $add_xref_sth->{'mysql_insertid'};

}


sub add_to_xrefs{
  my ($self,$master_xref,$acc,$version,$label,$description,$linkage,$source_id,$species_id) = @_;

  if(!defined($add_xref_sth)){
    $add_xref_sth = dbi->prepare("INSERT INTO xref (accession,version,label,description,source_id,species_id) VALUES(?,?,?,?,?,?)");
    $add_dependent_xref_sth = dbi->prepare("INSERT INTO dependent_xref VALUES(?,?,?,?)");
  }

  my $dependent_id = $self->get_xref($acc, $source_id);
  if(!defined($dependent_id)){
    $add_xref_sth->execute($acc,$version,$label,$description,$source_id,$species_id) || die "$acc\t$label\t\t$source_id\t$species_id\n";
  }
  $dependent_id = $self->get_xref($acc, $source_id);
  if(!defined($dependent_id)){
    die "$acc\t$label\t\t$source_id\t$species_id\n";
  }
  $add_dependent_xref_sth->execute($master_xref, $dependent_id,  $linkage, $source_id)|| die "$master_xref\t$dependent_id\t$linkage\t$source_id";


}

sub add_to_syn{
  my ($self, $acc, $source_id, $syn) = @_;

  if(!defined($add_synonym_sth)){
    $add_synonym_sth =  $dbi->prepare("INSERT INTO synonym VALUES(?,?)");
  }
  my $xref_id = $self->get_xref($acc, $source_id);
  if(defined($xref_id)){
    $add_synonym_sth->execute($xref_id, $syn) || die "$dbi->errstr \n $xref_id\n $syn\n";
  }
  else{
    die "Could not find acc $acc in xref table source = $source_id\n";
  }
}

# --------------------------------------------------------------------------------
# Add a single record to the direct_xref table.
# Note that an xref must already have been added to the xref table (xref_id passed as 1st arg)

sub add_direct_xref {

  my ($self, $general_xref_id, $ensembl_stable_id, $ensembl_type, $linkage_type) = @_;

  $add_direct_xref_sth = dbi->prepare("INSERT INTO direct_xref VALUES(?,?,?,?)") if (!defined($add_direct_xref_sth));

  $add_direct_xref_sth->execute($general_xref_id, $ensembl_stable_id, $ensembl_type, $linkage_type);

}

# --------------------------------------------------------------------------------

# Remove potentially problematic characters from string used for file/dir names

sub sanitise {

  my $str = shift;
  my $ret = $str;

  $ret =~ s/[\/\:]//g;
  return $ret;

}

# --------------------------------------------------------------------------------

# Create database if required. Assumes sql/table.sql and sql/populate_metadata.sql
# are present.

sub create {

  my ($host, $port, $user, $pass, $dbname, $sql_dir,$drop_db ) = @_;

  my $dbh = DBI->connect( "DBI:mysql:host=$host:port=$port", $user, $pass,
                          {'RaiseError' => 1});

  # check to see if the database already exists
  my %dbs = map {$_->[0] => 1} @{$dbh->selectall_arrayref('SHOW DATABASES')};

  if ($dbs{$dbname}) {

    if ( $drop_db ) {     
	$dbh->do( "DROP DATABASE $dbname" );
	print "Database $dbname dropped\n" ; 
    }
  
    if ( $create && !$drop_db ) {
      print "WARNING: about to drop database $dbname on $host:$port; yes to confirm, otherwise exit: ";
      $| = 1; # flush stdout
      my $p = <STDIN>;
      chomp $p;
      if ($p eq "yes") {
	$dbh->do( "DROP DATABASE $dbname" );
	print "Removed existing database $dbname\n";
      } else {
	print "$dbname NOT removed\n";
	exit(1);
      }
    } elsif ( !$create) {
      die("Database $dbname already exists. Use -create option to overwrite it.");
    }
  }

  $dbh->do( "CREATE DATABASE " . $dbname );

  print "Creating $dbname from ".$sql_dir."sql/table.sql\n";
  die "Cannot open  ".$sql_dir."sql/table.sql" if (! -e $sql_dir."sql/table.sql");
  my $cmd = "mysql -u $user -p$pass -P $port -h $host $dbname < ".$sql_dir."sql/table.sql";
  system ($cmd);

  print "Populating metadata in $dbname from ".$sql_dir."sql/populate_metadata.sql\n";
  die "Cannot open ".$sql_dir."sql/populate_metadata.sql" if (! -e $sql_dir."sql/populate_metadata.sql");
  $cmd = "mysql -u $user -p$pass -P $port -h $host $dbname < ".$sql_dir."sql/populate_metadata.sql";
  system($cmd);

}

sub get_label_to_accession{
  my ($self, $name) = @_;
  my %hash1=();

  my $dbi = dbi();
  my $sql = "select xref.accession, xref.label from xref, source where source.name like '$name%' and xref.source_id = source.source_id";
  my $sub_sth = dbi->prepare($sql);    

  $sub_sth->execute();
  while(my @row = $sub_sth->fetchrow_array()) {
    $hash1{$row[1]} = $row[0];
  }   	  
  return \%hash1;
}


sub get_accession_from_label{
  my ($self, $name) = @_;
  
  my $dbi = dbi();
  my $sql = "select xref.accession from xref where xref.label like '$name'";
  my $sub_sth = dbi->prepare($sql);    
  
  $sub_sth->execute();
  while(my @row = $sub_sth->fetchrow_array()) {
    return $row[0];
  }   	  
  return undef;
  
}

sub get_sub_list{
  my ($self, $name) = @_;
  my @list=();

  my $dbi = dbi();
  my $sql = "select xref.accession from xref where xref.accession like '$name%'";
  my $sub_sth = dbi->prepare($sql);    

  $sub_sth->execute();
  while(my @row = $sub_sth->fetchrow_array()) {
    push @list, $row[0];
  }   	  
  return @list;
}

# --------------------------------------------------------------------------------

# Set release for a source.

sub set_release {

  my ($release, $source_id) = @_;

  my $dbi = dbi();

  my $sth = $dbi->prepare("UPDATE source SET release=? WHERE source_id=?");

  $sth->execute($release, $source_id);

  print "Set release to $release for source ID $source_id\n";

}

# --------------------------------------------------------------------------------
1;

