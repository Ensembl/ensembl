package XrefParser::BaseParser;

use strict;

use DBI;
use Digest::MD5 qw(md5_hex);
use File::Path;
use POSIX qw(strftime);
use Getopt::Long;

my $base_dir = ".";

my $add_xref_sth = undef;
my $add_dependent_xref_sth = undef;
my $get_xref_sth = undef;

my $dbi;
my %dependent_sources;
my %taxonomy2species_id;
my %name2species_id;

my ($host, $port, $dbname, $user, $pass);
my $skipdownload;

# --------------------------------------------------------------------------------
# Get info about files to be parsed from the database

sub run {

  ($host, $port, $dbname, $user, $pass, my $speciesr, my $sourcesr, $skipdownload) = @_;

  my @species = @$speciesr;
  my @sources = @$sourcesr;

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
  while (my @row = $sth->fetchrow_array()) {


    # Download each source into the appropriate directory for parsing later
    # Also delete previous working directory if we're starting a new source type

    #can have more than one file

    my @files = split(/\s+/,$url);
    
    my $parse = 0;
    my $file_cs=0;
    my $type = $name;
    my @new_file=();
    $dir = $base_dir . "/" . $type;
    foreach my $urls (@files){
      my ($file) = $urls =~ /.*\/(.*)/;
      $last_type = $type;
      
      if (!$skipdownload) {
	
	rmtree $dir if ($type ne $last_type);
	mkdir $dir if (!-e $dir);
	
	print "Downloading $urls to $dir/$file\n";
	my $result = system("wget", "--quiet","--directory-prefix=$dir", $urls);
	
	# if the file is compressed, the FTP server may or may not have automatically uncompressed it
	# TODO - read .gz file directly? open (FILE, "zcat $file|") or Compress::Zlib
	if ($file =~ /(.*)\.gz$/) {
	  print "Uncompressing $dir/$file\n";
	  system("gunzip -f $dir/$file");
	  $file = $1;
	}
	
      }
      else{
	if ($file =~ /(.*)\.gz$/) {
	  $file = $1;
	}
      }
      push @new_file, $file;

      # compare checksums and parse/upload if necessary
      # need to check file size as some .SPC files can be of zero length
      print "HELLO: $dir/$file\n";
      $file_cs = md5sum("$dir/$file");
      if (!defined $checksum || $checksum ne $file_cs) {
	
	if (-s "$dir/$file") {
	  $parse =1;
	  print "Checksum for $file does not match, parsing\n";
	  
	  # Files from sources "UniProtSwissProt" and "UniProtSpTREMBL" are
	  # all parsed with the same parser
	  $parser = 'UniProtParser' if ($parser =~ /UniProt/i);
	}
	else {
	  
	  print $file . " has zero length, skipping\n";
	  
	}
      }
    }
    if($parse){

	print "Parsing ".join(' ',@new_file)." with $parser\n";
	eval "require XrefParser::$parser";
	my $new = "XrefParser::$parser"->new();
	$new->run("$dir/$new_file[0]", $source_id, $species_id);

	# update AFTER processing in case of crash.
	update_source($dbi, $source_url_id, $file_cs, $new_file[0]);

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
    warn("Couldn't get source ID for source name $source_name\n");
    $source_id = -1;
  }

  return $source_id;

}
sub get_valid_codes{
  my ($self,$source_name,$species_id) =@_;
  my %valid_codes;
  my @sources;

  my $sql = "select source_id from source where name like '%".$source_name."%'";
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
    }
  }
  return %valid_codes;
}

# --------------------------------------------------------------------------------
# Upload xrefs to the database

sub upload_xrefs {

  my ($self, @xrefs) = @_;

  my $dbi = dbi();

  if ($#xrefs > -1) {

    # remove all existing xrefs with same source ID(s)
    #delete_by_source(\@xrefs);

    # upload new ones
    print "Uploading xrefs\n";
    my $xref_sth = $dbi->prepare("INSERT INTO xref (accession,version,label,description,source_id,species_id) VALUES(?,?,?,?,?,?)");
    my $pri_insert_sth = $dbi->prepare("INSERT INTO primary_xref VALUES(?,?,?,?,?)");
    my $pri_update_sth = $dbi->prepare("UPDATE primary_xref SET sequence=? WHERE xref_id=?");
    my $syn_sth = $dbi->prepare("INSERT INTO synonym VALUES(?,?,?)");
    my $dep_sth = $dbi->prepare("INSERT INTO dependent_xref VALUES(?,?,?,?)");
    my $xref_update_label_sth = $dbi->prepare("UPDATE xref SET label=? WHERE xref_id=?");
    my $xref_update_descr_sth = $dbi->prepare("UPDATE xref SET description=? WHERE xref_id=?");

    local $xref_sth->{RaiseError}; # disable error handling here as we'll do it ourselves
    local $xref_sth->{PrintError};

    foreach my $xref (@xrefs) {
       my $xref_id;
      # Create entry in xref table and note ID
      if(! $xref_sth->execute($xref->{ACCESSION},
			 $xref->{VERSION},
			 $xref->{LABEL},
			 $xref->{DESCRIPTION},
			 $xref->{SOURCE_ID},
			 $xref->{SPECIES_ID})){
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
	print STDERR "xref_id is not set \n$xref->{ACCESSION}\n  $xref->{LABEL}\n $xref->{DESCRIPTION}\n $xref->{SOURCE_ID}\n";
      }
      if (primary_xref_id_exists($xref_id)) {
	
	$pri_update_sth->execute($xref->{SEQUENCE}, $xref_id) || die $dbi->errstr;
	
      } else {
	
	$pri_insert_sth->execute($xref_id,
				 $xref->{SEQUENCE},
				 $xref->{SEQUENCE_TYPE},
				 $xref->{STATUS},
				 $xref->{SOURCE_ID}) || die $dbi->errstr;
      }

      # if there are synonyms, create xrefs for them and entries in the synonym table
      foreach my $syn (@{$xref->{SYNONYMS}}) {

	$xref_sth->execute($syn,
			   "",
			   "",
			   "",
			   $xref->{SOURCE_ID},
			   $xref->{SPECIES_ID});

	my $syn_xref_id = insert_or_select($xref_sth, $dbi->err, $syn, $xref->{SOURCE_ID});

	$syn_sth->execute($xref_id, $syn_xref_id, $xref->{SOURCE_ID} ) || die $dbi->errstr;

      }				# foreach syn

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
	if(!defined($dep_xref_id)){
	  print STDERR "$dep{ACCESSION} \n $dep{SOURCE_ID} \n".$dbi->err."\n";
	}
	$dep_sth->execute($xref_id, $dep_xref_id, $dep{LINKAGE_ANNOTATION}, $dep{SOURCE_ID} ) || die $dbi->errstr;
	# TODO linkage anntation?

      }				# foreach dep


    }				# foreach xref


    $xref_sth->finish() if defined $xref_sth;
    $pri_insert_sth->finish() if defined $pri_insert_sth;
    $pri_update_sth->finish() if defined $pri_update_sth;

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

  if ($error) {

    $id = get_xref_id_by_accession_and_source($acc, $source);
#    print STDERR "Got existing xref id " . $id . " for " . $acc . " " . $source . "\n";

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

  # Remove dependent_xrefs and synonyms based on source of *xref*
  my $syn_sth = $dbi->prepare("DELETE FROM synonym USING xref, synonym WHERE xref.xref_id=synonym.xref_id AND xref.source_id=?");
  my $dep_sth = $dbi->prepare("DELETE FROM dependent_xref USING xref, dependent_xref WHERE xref.xref_id=dependent_xref.master_xref_id AND xref.source_id=?");

  # xrefs and primary_xrefs are straightforward deletes
  my $xref_sth = $dbi->prepare("DELETE FROM xref WHERE source_id=?");
  my $p_xref_sth = $dbi->prepare("DELETE FROM primary_xref WHERE source_id=?");

  # xrefs may come from more than one source (e.g. UniProt/SP/SPtr)
  # so find all sources first
  my %source_ids;
  foreach my $xref (@$xrefs) {
    my $xref_source = $xref->{SOURCE_ID};
    $source_ids{$xref_source} = 1;
  }

  # now delete them
  foreach my $source (keys %source_ids) {
    print "Deleting synonyms of xrefs with source ID $source \n";
    $syn_sth->execute($source);
    print "Deleting dependent xrefs of xrefs with source ID $source \n";
    $dep_sth->execute($source);
    print "Deleting primary xrefs with source ID $source \n";
    $p_xref_sth->execute($source);
    print "Deleting xrefs with source ID $source \n";
    $xref_sth->execute($source);
  }

  $syn_sth->finish() if defined $syn_sth;
  $dep_sth->finish() if defined $dep_sth;
  $xref_sth->finish() if defined $xref_sth;
  $p_xref_sth->finish() if defined $p_xref_sth;

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

sub get_xref{
  my ($acc,$source) = @_;

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

sub add_xref{
 my ($self,$acc,$version,$label,$description,$source_id,$species_id) = @_;

  if(!defined($add_xref_sth)){
    $add_xref_sth = dbi->prepare("INSERT INTO xref (accession,version,label,description,source_id,species_id) VALUES(?,?,?,?,?,?)");
  }
 $add_xref_sth->execute($acc,$version,$label,$description,$source_id,$species_id) 
   || die "$acc\t$label\t\t$source_id\t$species_id\n";
 
}


sub add_to_xrefs{
  my ($self,$master_xref,$acc,$version,$label,$description,$linkage,$source_id,$species_id) = @_;

  if(!defined($add_xref_sth)){
    $add_xref_sth = dbi->prepare("INSERT INTO xref (accession,version,label,description,source_id,species_id) VALUES(?,?,?,?,?,?)");
    $add_dependent_xref_sth = dbi->prepare("INSERT INTO dependent_xref VALUES(?,?,?,?)");
  }

  my $dependent_id = get_xref($acc, $source_id);
  if(!defined($dependent_id)){
    $add_xref_sth->execute($acc,$version,$label,$description,$source_id,$species_id) || die "$acc\t$label\t\t$source_id\t$species_id\n";
  }
  $dependent_id = get_xref($acc, $source_id);
  $add_dependent_xref_sth->execute($master_xref, $dependent_id,  $linkage, $source_id)|| die "$master_xref\t$dependent_id\t$linkage\t$source_id";


}

# --------------------------------------------------------------------------------

1;

