package XrefParser::BaseParser;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception;
use XrefParser::FetchFiles;
use XrefParser::Database;
use Carp;
use DBI;
use Getopt::Long;
#use POSIX qw(strftime);


#use File::Spec::Functions;
#use IO::File;
#use Net::FTP;
#use URI;
#use URI::file;
#use Text::Glob qw( match_glob );
#use LWP::UserAgent;


my $base_dir = File::Spec->curdir();

my $add_xref_sth = undef;
my %add_direct_xref_sth;
my $add_dependent_xref_sth = undef;
my $get_xref_sth = undef;
my $add_synonym_sth = undef;

my $dbi;
my %dependent_sources;
my %taxonomy2species_id;
my %species_id2taxonomy;
my %name2species_id;
my %species_id2name;
my %xref_dependent_mapped;


my $verbose;



###################################################
# Create new object.
#   set global $verbose
#   Store the dbi form the database for easy access
###################################################
sub new
{
  my ($proto, $database, $is_verbose) = @_;

  if((!defined $database)){# or (!$database->isa(XrefPArserDatabase)))
    croak "No database specfied";
  }
  $verbose = $is_verbose;
  $dbi = $database->dbi;

  my $class = ref $proto || $proto;
  my $self =  bless {}, $class;
  $self->dbi($dbi);
  return $self;
}


##################################
# Getter/Setter for the dbi object
##################################
sub dbi {
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_dbi} = $arg );
  return $self->{_dbi};
}


#######################################################################
# Given a file name, returns a IO::Handle object.  If the file is
# gzipped, the handle will be to an unseekable stream coming out of a
# zcat pipe.  If the given file name doesn't correspond to an existing
# file, the routine will try to add '.gz' to the file name or to remove
# any .'Z' or '.gz' and try again.  Returns undef on failure and will
# write a warning to stderr.
#######################################################################
sub get_filehandle
{
    my ($self, $file_name) = @_;

    my $io;

    my $alt_file_name = $file_name;
    $alt_file_name =~ s/\.(gz|Z)$//x;

    if ( $alt_file_name eq $file_name ) {
        $alt_file_name .= '.gz';
    }

    if ( !-e $file_name ) {
        carp(   "File '$file_name' does not exist, "
              . "will try '$alt_file_name'" );
        $file_name = $alt_file_name;
    }

    if ( $file_name =~ /\.(gz|Z)$/x ) {
        # Read from zcat pipe
        $io = IO::File->new("zcat $file_name |")
          or carp("Can not open file '$file_name' with 'zcat'");
    } else {
        # Read file normally
        $io = IO::File->new($file_name)
          or carp("Can not open file '$file_name'");
    }

    if ( !defined $io ) { return }

    if ($verbose) {
      print("Reading from '$file_name'...\n");
    }

    return $io;
}


#############################################
# Get source ID for a particular source name
#
# Arg[1] source name
# Arg[2] priority description
#
# Returns source_id or -1 if not found
#############################################
sub get_source_id_for_source_name {
  my ($self, $source_name,$priority_desc) = @_;

  my $sql = "SELECT source_id FROM source WHERE LOWER(name)='" . lc($source_name) . "'";
  if(defined($priority_desc)){
    $sql .= " AND LOWER(priority_description)='".lc($priority_desc)."'";
    $source_name .= " ($priority_desc)";
  }
  my $sth = $dbi->prepare($sql);
  $sth->execute();
  my @row = $sth->fetchrow_array();
  my $source_id;
  if (@row) {
    $source_id = $row[0]; 
  } else {
    carp "WARNING: There is no entity $source_name in the source-table of the xref database.\n";
    carp "WARNING:. The external db name ($source_name) is hardcoded in the parser\n";
    carp "WARNING: Couldn't get source ID for source name $source_name\n";

    $source_id = -1;
  }
  return $source_id;
}



############################################################
# Get a set of source IDs matching a source name pattern
#
# Adds % to each end of the source name and doe a like query
# to find all the matching source names source_ids.
#
# Returns an empty list if none found.
############################################################
sub get_source_ids_for_source_name_pattern {

  my ($self, $source_name) = @_;

  my $sql = "SELECT source_id FROM source WHERE upper(name) LIKE '%".uc($source_name)."%'";

  my $sth = $dbi->prepare($sql);
  my @sources;
  $sth->execute();
  while(my @row = $sth->fetchrow_array()){
    push @sources,$row[0];
  }
  $sth->finish;

  return @sources;

}


###############################
# From a source_id get the name
###############################
sub get_source_name_for_source_id {
  my ($self, $source_id) = @_;
  my $source_name;

  my $sql = "SELECT name FROM source WHERE source_id= '" . $source_id. "'";
  my $sth = $dbi->prepare($sql);
  $sth->execute();
  my @row = $sth->fetchrow_array();
  if (@row) {
    $source_name = $row[0]; 
  } else {
    carp "There is no entity with source-id  $source_id  in the source-table of the \n";
    carp "xref-database. The source-id and the name of the source-id is hard-coded in populate_metadata.sql\n" ;
    carp "and in the parser\n";
    carp "Couldn't get source name for source ID $source_id\n";
    $source_name = -1;
  }
  return $source_name;
}



####################################################
# Get a hash to go from accession of a dependent xref
# to master_xref_id for all of source names given
#####################################################
sub get_valid_xrefs_for_dependencies{
  my ($self, $dependent_name, @reverse_ordered_source_list) = @_;

  my %dependent_2_xref;


  my $sql = "select source_id from source where LOWER(name) =?";
  my $sth = $dbi->prepare($sql);
  my @dependent_sources;
  $sth->execute(lc($dependent_name));
  while(my @row = $sth->fetchrow_array()){
   push @dependent_sources,$row[0];
  }

  my @sources;
  foreach my $name (@reverse_ordered_source_list){
    $sth->execute(lc($name));
    while(my @row = $sth->fetchrow_array()){
      push @sources,$row[0];
    }
  }
  $sth->finish;

  my $dep_sql = (<<"DSS");
SELECT d.master_xref_id, x2.accession
  FROM dependent_xref d, xref x1, xref x2
    WHERE x1.xref_id = d.master_xref_id AND
          x1.source_id = ? AND
          x2.xref_id = d.dependent_xref_id AND
          x2.source_id = ?
DSS

  $sth = $dbi->prepare($dep_sql);
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



####################################################
# Get a hash to go from accession of a direct xref
# to master_xref_id for all of source names given
#####################################################
sub get_valid_xrefs_for_direct_xrefs{
  my ($self, $direct_name, @list) = @_;

  my %direct_2_xref;


  my $sql = "select source_id from source where name like ?";
  my $sth = $dbi->prepare($sql);
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

  my $gen_sql =(<<"GDS");
SELECT d.general_xref_id, d.ensembl_stable_id, 'TYPE', d.linkage_xref, x1.accession
  FROM TABLE_direct_xref d, xref x1
    WHERE x1.xref_id = d.general_xref_id AND
          x1.source_id=?
GDS

  my @sth;
  my $i=0;
  foreach my $type (qw(Gene Transcript Translation)){
    my $t_sql = $gen_sql;
    my $table = lc($type);
    $t_sql =~ s/TABLE/$table/x;
    $t_sql =~ s/TYPE/$type/x;

    $sth[$i++] = $dbi->prepare($t_sql);
  }

  foreach my $d (@direct_sources){
    for (my $ii =0; $i<3; $i++){
      $sth[$ii]->execute($d);
      while(my @row = $sth[$ii]->fetchrow_array()){
	$direct_2_xref{$row[4]} = $row[0]."::".$row[1]."::".$row[2]."::".$row[3];
      }
    }
  }

  return \%direct_2_xref;
}


#############################################
# Get a hash of label to acc for a particular
# source name and species_id
#############################################
sub label_to_acc{

  my ($self,$source_name,$species_id) =@_;

  # First cache synonyms so we can quickly add them later
  my %synonyms;
  my $syn_sth = $dbi->prepare("SELECT xref_id, synonym FROM synonym");
  $syn_sth->execute();

  my ($xref_id, $synonym);
  $syn_sth->bind_columns(\$xref_id, \$synonym);
  while ($syn_sth->fetch()) {

    push @{$synonyms{$xref_id}}, $synonym;

  }

  my %valid_codes;
  my @sources;

  my $sql = "select source_id from source where upper(name) like '%".uc($source_name)."%'";
  my $sth = $dbi->prepare($sql);
  $sth->execute();
  while(my @row = $sth->fetchrow_array()){
    push @sources,$row[0];
  }
  $sth->finish;

  foreach my $source (@sources){
    $sql = "select label, xref_id from xref where species_id = $species_id and source_id = $source";
    $sth = $dbi->prepare($sql);
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


####################################################
# get_valid_codes
#
# hash of accession to array of xrefs.
# This is an array becouse more than one entry can
# exist. i.e. for uniprot and refseq we have direct
# and sequence match sets and we need to give both.
####################################################
sub get_valid_codes{

  my ($self,$source_name,$species_id) =@_;

  # First cache synonyms so we can quickly add them later
  my %synonyms;
  my $syn_sth = $dbi->prepare("SELECT xref_id, synonym FROM synonym");
  $syn_sth->execute();

  my ($xref_id, $synonym);
  $syn_sth->bind_columns(\$xref_id, \$synonym);
  while ($syn_sth->fetch()) {

    push @{$synonyms{$xref_id}}, $synonym;

  }

  my %valid_codes;
  my @sources;

  my $sql = "select source_id from source where upper(name) like '%".uc($source_name)."%'";
  my $sth = $dbi->prepare($sql);
  $sth->execute();
  while(my @row = $sth->fetchrow_array()){
    push @sources,$row[0];
  }
  $sth->finish;

  foreach my $source (@sources){
    $sql = "select accession, xref_id from xref where species_id = $species_id and source_id = $source";
    $sth = $dbi->prepare($sql);
    $sth->execute();
    while(my @row = $sth->fetchrow_array()){
      push @{$valid_codes{$row[0]}}, $row[1];
      # add any synonyms for this xref as well
      foreach my $syn (@{$synonyms{$row[1]}}) {
	push @{$valid_codes{$syn}}, $row[1];
      }
    }
  }
  return \%valid_codes;
}

##############################
# Upload xrefs to the database
##############################
sub upload_xref_object_graphs {
  my ($self, $rxrefs) = @_;

  print "count = ".$#$rxrefs."\n" if($verbose);

  if ($#$rxrefs > -1) {

    # upload new ones
    print "Uploading xrefs\n" if($verbose);
    my $xref_sth = $dbi->prepare("INSERT INTO xref (accession,version,label,description,source_id,species_id, info_type) VALUES(?,?,?,?,?,?,?)");
    my $pri_insert_sth = $dbi->prepare("INSERT INTO primary_xref VALUES(?,?,?,?)");
    my $pri_update_sth = $dbi->prepare("UPDATE primary_xref SET sequence=? WHERE xref_id=?");
    my $syn_sth = $dbi->prepare("INSERT INTO synonym VALUES(?,?)");
    my $dep_sth = $dbi->prepare("INSERT INTO dependent_xref (master_xref_id, dependent_xref_id, linkage_annotation, linkage_source_id) VALUES(?,?,?,?)");
    my $xref_update_label_sth = $dbi->prepare("UPDATE xref SET label=? WHERE xref_id=?");
    my $xref_update_descr_sth = $dbi->prepare("UPDATE xref SET description=? WHERE xref_id=?");
    my $pair_sth = $dbi->prepare("INSERT INTO pairs VALUES(?,?,?)");

    local $xref_sth->{RaiseError} = 0; # disable error handling here as we'll do it ourselves
    local $xref_sth->{PrintError} = 0;

    foreach my $xref (@{$rxrefs}) {
       my $xref_id=undef;
       if(!defined($xref->{ACCESSION})){
	 print "your xref does not have an accession-number,so it can't be stored in the database\n";
	 return;
       }
      # Create entry in xref table and note ID
       if(! $xref_sth->execute($xref->{ACCESSION},
			 $xref->{VERSION} || 0,
			 $xref->{LABEL}|| $xref->{ACCESSION},
			 $xref->{DESCRIPTION},
			 $xref->{SOURCE_ID},
			 $xref->{SPECIES_ID},
			 $xref->{INFO_TYPE} || "MISC")){
	 if(!defined($xref->{SOURCE_ID})){
	   print "your xref: $xref->{ACCESSION} does not have a source-id\n";
	   return;
	 }
	 $xref_id = $self->insert_or_select($xref_sth, $dbi->err, $xref->{ACCESSION}, $xref->{SOURCE_ID}, $xref->{SPECIES_ID});
	 $xref_update_label_sth->execute($xref->{LABEL},$xref_id) if (defined($xref->{LABEL}));
	 $xref_update_descr_sth->execute($xref->{DESCRIPTION},$xref_id,) if (defined($xref->{DESCRIPTION}));
       }
       else{
	 $xref_id = $self->insert_or_select($xref_sth, $dbi->err, $xref->{ACCESSION}, $xref->{SOURCE_ID}, $xref->{SPECIES_ID});
       }

       # create entry in primary_xref table with sequence; if this is a "cumulative"
       # entry it may already exist, and require an UPDATE rather than an INSERT
       if(defined($xref->{SEQUENCE})){
	 if(!(defined($xref_id) and $xref_id)){
	   print STDERR "xref_id is not set for :\n$xref->{ACCESSION}\n$xref->{LABEL}\n$xref->{DESCRIPTION}\n$xref->{SOURCE_ID}\n$xref->{SPECIES_ID}\n";
	 }
	 if ( primary_xref_id_exists($xref_id) ) {
	   $pri_update_sth->execute( $xref->{SEQUENCE}, $xref_id )
	     or croak( $dbi->errstr() );
	 } else {

	   $pri_insert_sth->execute( $xref_id, $xref->{SEQUENCE},
				     $xref->{SEQUENCE_TYPE},
				     $xref->{STATUS} )
	     or croak( $dbi->errstr() );
	 }
       }
       # if there are synonyms, add entries in the synonym table
       foreach my $syn ( @{ $xref->{SYNONYMS} } ) {
	 $syn_sth->execute( $xref_id, $syn )
            or croak( $dbi->errstr() . "\n $xref_id\n $syn\n" );
       } # foreach syn

      # if there are dependent xrefs, add xrefs and dependent xrefs for them
      foreach my $depref (@{$xref->{DEPENDENT_XREFS}}) {

	my %dep = %$depref;

	$xref_sth->execute($dep{ACCESSION},
			   $dep{VERSION} || 0,
			   $dep{LABEL} || $dep{ACCESSION},
			   $dep{DESCRIPTION} || "",
			   $dep{SOURCE_ID},
			   $xref->{SPECIES_ID},
                           "DEPENDENT");

	my $dep_xref_id = $self->insert_or_select($xref_sth, $dbi->err, $dep{ACCESSION}, $dep{SOURCE_ID}, $xref->{SPECIES_ID});

	if(!defined($dep_xref_id) || $dep_xref_id ==0 ){
	  print STDERR "acc = $dep{ACCESSION} \nlink = $dep{LINKAGE_SOURCE_ID} \n".$dbi->err."\n";
	  print STDERR "source = $dep{SOURCE_ID}\n";
	}
        $dep_sth->execute( $xref_id, $dep_xref_id,
            $dep{LINKAGE_ANNOTATION},
            $dep{LINKAGE_SOURCE_ID} )
          or croak( $dbi->errstr() );
	# if there are synonyms, add entries in the synonym table
	foreach my $syn ( @{ $dep{SYNONYMS} } ) {
	  $syn_sth->execute( $dep_xref_id, $syn )
            or croak( $dbi->errstr() . "\n $xref_id\n $syn\n" );
        } # foreach syn
       } # foreach dep

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

    my $general_xref_id = get_xref($dr->{ACCESSION},$dr->{SOURCE_ID},$dr->{SPECIES_ID});
    if ($general_xref_id){
      $self->add_direct_xref($general_xref_id, $dr->{ENSEMBL_STABLE_ID},$dr->{ENSEMBL_TYPE},$dr->{LINKAGE_XREF});
    }
  }
  return;
}

sub add_meta_pair {

  my ($self, $key, $value) = @_;

  my $sth = $dbi->prepare('insert into meta (meta_key, meta_value, date) values("'.$key.'", "'.$value.'", now())');
  $sth->execute;
  $sth->finish;
  return;
}



# --------------------------------------------------------------------------------
# Get & cache a hash of all the source names for dependent xrefs (those that are
# in the source table but don't have an associated URL etc)

sub get_dependent_xref_sources {

  my $self = shift;

  if (!%dependent_sources) {

    my $sth = $dbi->prepare("SELECT name,source_id FROM source");
    $sth->execute() or croak( $dbi->errstr() );
    while(my @row = $sth->fetchrow_array()) {
      my $source_name = $row[0];
      my $source_id = $row[1];
      $dependent_sources{$source_name} = $source_id;
    }
    $sth->finish;
  }

  return %dependent_sources;

}

# --------------------------------------------------------------------------------
# Get & cache a hash of all the species IDs & taxonomy IDs.

sub taxonomy2species_id {
  carp ( "[DEPRECATED] taxonomy2species_id is a deprecated (unsafe) method. ".
        "Please use species_id2taxonomy instead. Called by ".
        join( ', ', (caller(0))[1..2] ) );

  my $self = shift;

  if (!%taxonomy2species_id) {

    my $sth = $dbi->prepare("SELECT species_id, taxonomy_id FROM species");
    $sth->execute() or croak( $dbi->errstr() );
    while(my @row = $sth->fetchrow_array()) {
      my $species_id = $row[0];
      my $taxonomy_id = $row[1];
      if( my $ori =$taxonomy2species_id{$taxonomy_id} ){
        croak ( "Taxon $taxonomy_id already used for species $ori. ".
               "Cannot assign to species $species_id as well. ".
               "Consider using the species_id2taxonomy call instead. ".
               "Called by ". join( ', ', (caller(0))[1..2] ) );
      }
      $taxonomy2species_id{$taxonomy_id} = $species_id;
    }
  }

  return %taxonomy2species_id;

}


sub species_id2taxonomy {

  my $self = shift;

  if (!%species_id2taxonomy) {

    my $sth = $dbi->prepare("SELECT species_id, taxonomy_id FROM species");
    $sth->execute() or croak( $dbi->errstr() );
    while(my @row = $sth->fetchrow_array()) {
      my $species_id = $row[0];
      my $taxonomy_id = $row[1];
      if(defined($species_id2taxonomy{$species_id})){
	push @{$species_id2taxonomy{$species_id}}, $taxonomy_id;
      }
      else{
	$species_id2taxonomy{$species_id} = [$taxonomy_id];
      }
    }
  }
  return %species_id2taxonomy;

}

# --------------------------------------------------------------------------------
# Get & cache a hash of all the species IDs & species names.

sub name2species_id {
  carp ( "[DEPRECATED] name2species_id is a deprecated (unsafe) method. ".
        "Please use species_id2name instead. Called by ".
        join( ', ', (caller(0))[1..2] ) );

    my $self = shift;

    if ( !%name2species_id ) {

        my $sth = $dbi->prepare("SELECT species_id, name FROM species");
        $sth->execute() or croak( $dbi->errstr() );
        while ( my @row = $sth->fetchrow_array() ) {
            my $species_id = $row[0];
            my $name       = $row[1];
            $name2species_id{$name} = $species_id;
        }

        # Also populate the hash with all the aliases.
        $sth = $dbi->prepare("SELECT species_id, aliases FROM species");
        $sth->execute() or croak( $dbi->errstr() );
        while ( my @row = $sth->fetchrow_array() ) {
            my $species_id = $row[0];
            foreach my $name ( split /,\s*/x, $row[1] ) {
                if ( my $ori = $name2species_id{$name} ) {
                  croak ( "Name $name already used for species $ori. ".
                       "Cannot assign to species $species_id as well. ".
                       "Consider using the species_id2name call instead. ".
                       "Called by ". join( ', ', (caller(0))[1..2] ) );
                } else {
                    $name2species_id{$name} = $species_id;
                }
            }
        }

    } ## end if ( !%name2species_id)

    return %name2species_id;
} ## end sub name2species_id

sub species_id2name {
  my $self = shift;

  if ( !%species_id2name ) {

    my $sth = $dbi->prepare("SELECT species_id, name FROM species");
    $sth->execute() or croak( $dbi->errstr() );
    while ( my @row = $sth->fetchrow_array() ) {
      my $species_id = $row[0];
      my $name       = $row[1];
      $species_id2name{$species_id} = [ $name ];
    }
    
    # Also populate the hash with all the aliases.
    $sth = $dbi->prepare("SELECT species_id, aliases FROM species");
    $sth->execute() or croak( $dbi->errstr() );
    while ( my @row = $sth->fetchrow_array() ) {
      my $species_id = $row[0];
      foreach my $name ( split /,\s*/x, $row[1] ) {
        $species_id2name{$species_id} ||= [];
        push @{$species_id2name{$species_id}}, $name;
      }
    }  
  } ## end if ( !%species_id2name)

  return %species_id2name;
} ## end sub species_id2name

# --------------------------------------------------------------------------------
# Update a row in the source table

# --------------------------------------------------------------------------------
sub dbi2{

    my ($self, $host2, $port2, $user2, $dbname2, $pass2) = @_;
    my $dbi2 = undef;

    if ( !defined $dbi2 || !$dbi2->ping() ) {
        my $connect_string =
          sprintf( "dbi:mysql:host=%s;port=%s;database=%s",
            $host2, $port2||3306, $dbname2 );

        $dbi2 =
          DBI->connect( $connect_string, $user2, $pass2)
          or carp( "Can't connect to database: " . $DBI::errstr ) and return;
        $dbi2->{'mysql_auto_reconnect'} = 1; # Reconnect on timeout
    }

    return $dbi2;
}

# --------------------------------------------------------------------------------


# --------------------------------------------------------------------------------

# --------------------------------------------------------------------------------
# If there was an error, an xref with the same acc & source already exists.
# If so, find its ID, otherwise get ID of xref just inserted

sub insert_or_select {

  my ($self, $sth, $error, $acc, $source, $species) = @_;

  my $id;

  if ($error and ($error == 1062)) {  # duplicate (okay so get the original) 

    $id = $self->get_xref($acc, $source, $species);
	
  } 
  elsif ($error){
    die "Error $error";
  }
  else {
	
    $id = $sth->{'mysql_insertid'};
	
  }

  return $id;

}

# --------------------------------------------------------------------------------

sub primary_xref_id_exists {

  my $xref_id = shift;

  my $exists = 0;

  my $sth = $dbi->prepare("SELECT xref_id FROM primary_xref WHERE xref_id=?");
  $sth->execute($xref_id) or croak( $dbi->errstr() );
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
    print "Deleting pairs with source ID $source \n" if($verbose);
    $pairs_sth->execute($source);
    print "Deleting direct xrefs with source ID $source \n" if($verbose);
    $direct_sth->execute($source);
    print "Deleting synonyms of xrefs with source ID $source \n" if($verbose);
    $syn_sth->execute($source);
    print "Deleting dependent xrefs of xrefs with source ID $source \n" if($verbose);
    $dep_sth->execute($source);
    print "Deleting primary xrefs with source ID $source \n" if($verbose);
#    $p_xref_sth->execute($source);
    print "Deleting xrefs with source ID $source \n" if($verbose);
    $xref_sth->execute($source);
  }

  $syn_sth->finish() if defined $syn_sth;
  $dep_sth->finish() if defined $dep_sth;
  $xref_sth->finish() if defined $xref_sth;
#  $p_xref_sth->finish() if defined $p_xref_sth;
  return;
}

sub get_taxonomy_from_species_id{
  my ($self,$species_id) = @_;
  my %hash;

  my $sth = $dbi->prepare("SELECT taxonomy_id FROM species WHERE species_id = $species_id");
  $sth->execute() or croak( $dbi->errstr() );
  while(my @row = $sth->fetchrow_array()) {
    $hash{$row[0]} = 1;
  }   
  $sth->finish;
  return \%hash;
}
# --------------------------------------------------------------------------------


sub get_direct_xref{
 my ($self,$stable_id,$type,$link) = @_;

 $type = lc($type);

 my $sql = "select general_xref_id from ${type}_direct_xref d where ensembl_stable_id = ?  and linkage_xref= ?";
 my  $direct_sth = $dbi->prepare($sql);

 $direct_sth->execute( $stable_id, $link ) or croak( $dbi->errstr() );
 if(my @row = $direct_sth->fetchrow_array()) {
   return $row[0];
 }
 return;
}


sub get_xref{
  my ($self,$acc,$source, $species_id) = @_;

  if(!defined($get_xref_sth)){
    my $sql = "select xref_id from xref where accession = ? and source_id = ? and species_id = ?";
    $get_xref_sth = $dbi->prepare($sql);  
  }
  
  $get_xref_sth->execute( $acc, $source, $species_id ) or croak( $dbi->errstr() );
  if(my @row = $get_xref_sth->fetchrow_array()) {
    return $row[0];
  }   
  return;
}

sub add_xref {
  my ( $self, $acc, $version, $label, $description, $source_id,
       $species_id, $info_type )
    = @_;

  my $xref_id = $self->get_xref($acc,$source_id, $species_id);
  if(defined($xref_id)){
    return $xref_id;
  }
  if ( !defined($add_xref_sth) ) {
    $add_xref_sth =
      $dbi->prepare( "INSERT INTO xref "
         . "(accession,version,label,description,source_id,species_id, info_type) "
         . "VALUES(?,?,?,?,?,?,?)" );
  }

  # If the description is more than 255 characters, chop it off and add
  # an indication that it has been truncated to the end of it.

  if (defined($description) && (length($description) > 255 ) ) {
    my $truncmsg = ' /.../';
    substr( $description, 255 - length($truncmsg),
            length($truncmsg), $truncmsg );
  }

  $add_xref_sth->execute( $acc, $version || 0, $label,
                          $description, $source_id, $species_id, $info_type
  ) or croak("$acc\t$label\t\t$source_id\t$species_id\n");

  return $add_xref_sth->{'mysql_insertid'};
} ## end sub add_xref


sub add_to_direct_xrefs{
  my ($self,$direct_xref,$type, $acc,$version,$label,$description,$linkage,$source_id,$species_id) = @_;

  $direct_xref || croak( "Need a direct_xref on which this xref linked too" );
  $acc         || croak( "Need an accession of this dependent xref" );
  $version     ||= 0;
  $label       ||= $acc;
  $description ||= undef;
  $linkage     ||= undef;
  $source_id   || croak( "Need a source_id for this dependent xref" );
  $species_id  || croak( "Need a species_id for this dependent xref" );

  if(!defined($add_xref_sth)){
    my $sql = (<<"AXX");
INSERT INTO xref
  (accession,version,label,description,source_id,species_id, info_type)
  VALUES (?,?,?,?,?,?,?)
AXX
    $add_xref_sth = $dbi->prepare($sql);
  }


  my $direct_id = $self->get_xref($acc, $source_id, $species_id);
  if(!defined($direct_id)){
    $add_xref_sth->execute(
        $acc, $version || 0, $label,
        $description, $source_id, $species_id, "DIRECT"
    ) or croak("$acc\t$label\t\t$source_id\t$species_id\n");
  }
  $direct_id = $self->get_xref($acc, $source_id, $species_id);

  $self->add_direct_xref($direct_id, $direct_xref, $type, "");
  return;
}

sub add_to_xrefs{
  my ($self,$master_xref,$acc,$version,$label,$description,$linkage,$source_id,$species_id) = @_;

  $master_xref || croak( "Need a master_xref_id on which this xref depends" );
  $acc         || croak( "Need an accession of this dependent xref" );
  $version     ||= 0;
  $label       ||= $acc;
  $description ||= undef;
  $linkage     ||= undef;
  $source_id   || croak( "Need a source_id for this dependent xref" );
  $species_id  || croak( "Need a species_id for this dependent xref" );

  if(!defined($add_xref_sth)){
    my $sql = (<<"IXR");
INSERT INTO xref 
  (accession,version,label,description,source_id,species_id, info_type)
  VALUES (?,?,?,?,?,?,?)
IXR
    $add_xref_sth = $dbi->prepare($sql);
  }
  if(!defined($add_dependent_xref_sth)){
    my $sql = (<<"ADX");
INSERT INTO dependent_xref 
  (master_xref_id,dependent_xref_id,linkage_annotation,linkage_source_id)
  VALUES (?,?,?,?)
ADX
    $add_dependent_xref_sth = $dbi->prepare($sql);
  }

  my $dependent_id = $self->get_xref($acc, $source_id, $species_id);
  if(!defined($dependent_id)){
    $add_xref_sth->execute(
        $acc, $version || 0, $label,
        $description, $source_id, $species_id, "DEPENDENT"
    ) or croak("$acc\t$label\t\t$source_id\t$species_id\n");
  }
  $dependent_id = $self->get_xref($acc, $source_id, $species_id);
  if(!defined($dependent_id)){
    croak("$acc\t$label\t\t$source_id\t$species_id\n");
  }

  if(!defined($xref_dependent_mapped{$master_xref."|".$dependent_id})){
    $add_dependent_xref_sth->execute( $master_xref, $dependent_id, $linkage,
				      $source_id )
      or croak("$master_xref\t$dependent_id\t$linkage\t$source_id");
    $xref_dependent_mapped{$master_xref."|".$dependent_id} = 1;
  }
  return $dependent_id;
}

sub add_to_syn_for_mult_sources{
  my ($self, $acc, $sources, $syn, $species_id) = @_;

  if(!defined($add_synonym_sth)){
    $add_synonym_sth =  $dbi->prepare("INSERT INTO synonym VALUES(?,?)");
  }
  my $found =0;
  foreach my $source_id (@$sources){
    my $xref_id = $self->get_xref($acc, $source_id, $species_id);
    if(defined($xref_id)){
      $add_synonym_sth->execute( $xref_id, $syn )
        or croak( $dbi->errstr() . "\n $xref_id\n $syn\n" );
      $found = 1;
    }
  }
  return;
}


sub add_to_syn{
  my ($self, $acc, $source_id, $syn, $species_id) = @_;

  if(!defined($add_synonym_sth)){
    $add_synonym_sth =  $dbi->prepare("INSERT INTO synonym VALUES(?,?)");
  }
  my $xref_id = $self->get_xref($acc, $source_id, $species_id);
  if(defined($xref_id)){
    $add_synonym_sth->execute( $xref_id, $syn )
      or croak( $dbi->errstr() . "\n $xref_id\n $syn\n" );
  }
  else {
      carp (  "Could not find acc $acc in "
            . "xref table source = $source_id of species $species_id\n" );
  }
  return;
}


sub add_synonym{
  my ($self, $xref_id, $syn) = @_;

  if(!defined($add_synonym_sth)){
    $add_synonym_sth =  $dbi->prepare("INSERT INTO synonym VALUES(?,?)");
  }

    $add_synonym_sth->execute( $xref_id, $dbi->quote($syn) )
      or croak( $dbi->errstr() . "\n $xref_id\n " . $dbi->quote($syn) . "\n" );

  return;
}



# --------------------------------------------------------------------------------
# Add a single record to the direct_xref table.
# Note that an xref must already have been added to the xref table (xref_id passed as 1st arg)

sub add_direct_xref {

  my ($self, $general_xref_id, $ensembl_stable_id, $ensembl_type, $linkage_type) = @_;

  if (!defined($add_direct_xref_sth{$ensembl_type})){
    my $add_gene_direct_xref_sth = $dbi->prepare("INSERT INTO gene_direct_xref VALUES(?,?,?)");
    my $add_tr_direct_xref_sth = $dbi->prepare("INSERT INTO transcript_direct_xref VALUES(?,?,?)");
    my $add_tl_direct_xref_sth = $dbi->prepare("INSERT INTO translation_direct_xref VALUES(?,?,?)");
    $add_direct_xref_sth{"gene"} = $add_gene_direct_xref_sth;
    $add_direct_xref_sth{"transcript"} = $add_tr_direct_xref_sth;
    $add_direct_xref_sth{"translation"} = $add_tl_direct_xref_sth;
    $add_direct_xref_sth{"Gene"} = $add_gene_direct_xref_sth;
    $add_direct_xref_sth{"Transcript"} = $add_tr_direct_xref_sth;
    $add_direct_xref_sth{"Translation"} = $add_tl_direct_xref_sth;
  }
  
  if(!defined($add_direct_xref_sth{$ensembl_type})){
    print "ERROR add_direct_xref_sth does not exist for $ensembl_type ???\n"; 
  }
  else{
    $add_direct_xref_sth{$ensembl_type}->execute($general_xref_id, $ensembl_stable_id, $linkage_type);
  }
  return;
}

# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------

sub get_label_to_acc{
  my ($self, $name, $species_id, $prio_desc) = @_;
  my %hash1=();

  my $sql = "select xref.accession, xref.label from xref, source where source.name like '$name%' and xref.source_id = source.source_id";
  if(defined($prio_desc)){
    $sql .= " and source.priority_description like '".$prio_desc."'";
  }
  if(defined($species_id)){
    $sql .= " and xref.species_id  = $species_id";
  }
  my $sub_sth = $dbi->prepare($sql);

  $sub_sth->execute();
  while(my @row = $sub_sth->fetchrow_array()) {
    $hash1{$row[1]} = $row[0];
  }


  #
  # Remember synonyms
  #

 $sql = "select xref.accession, synonym.synonym from xref, source, synonym where synonym.xref_id = xref.xref_id and source.name like '$name%' and xref.source_id = source.source_id";
  if(defined($prio_desc)){
    $sql .= " and source.priority_description like '".$prio_desc."'";
  }
  if(defined($species_id)){
    $sql .= " and xref.species_id  = $species_id";
  }
  $sub_sth = $dbi->prepare($sql);    

  $sub_sth->execute();
  while(my @row = $sub_sth->fetchrow_array()) {
    $hash1{$row[1]} = $row[0];
  }

  return \%hash1;
}


sub get_label_to_desc{
  my ($self, $name, $species_id, $prio_desc) = @_;
  my %hash1=();

  my $sql = "select xref.description, xref.label from xref, source where source.name like '$name%' and xref.source_id = source.source_id";
  if(defined($prio_desc)){
    $sql .= " and source.priority_description like '".$prio_desc."'";
  }
  if(defined($species_id)){
    $sql .= " and xref.species_id  = $species_id";
  }
  my $sub_sth = $dbi->prepare($sql);

  $sub_sth->execute();
  while(my @row = $sub_sth->fetchrow_array()) {
    $hash1{$row[1]} = $row[0];
  }

  #
  # Also include the synonyms
  #

  $sql = "select xref.description, synonym.synonym from xref, source, synonym where synonym.xref_id = xref.xref_id and source.name like '$name%' and xref.source_id = source.source_id";
  if(defined($prio_desc)){
    $sql .= " and source.priority_description like '".$prio_desc."'";
  }
  if(defined($species_id)){
    $sql .= " and xref.species_id  = $species_id";
  }
  $sub_sth = $dbi->prepare($sql);

  $sub_sth->execute();
  while(my @row = $sub_sth->fetchrow_array()) {
    $hash1{$row[1]} = $row[0];
  }

  return \%hash1;
}



sub get_accession_from_label{
  my ($self, $name) = @_;

  my $sql = "select xref.accession from xref where xref.label like '$name'";
  my $sub_sth = $dbi->prepare($sql);

  $sub_sth->execute();
  while(my @row = $sub_sth->fetchrow_array()) {
    return $row[0];
  }
  return;
}

sub get_sub_list{
  my ($self, $name) = @_;
  my @list=();

  my $sql = "select xref.accession from xref where xref.accession like '$name%'";
  my $sub_sth = $dbi->prepare($sql);

  $sub_sth->execute();
  while(my @row = $sub_sth->fetchrow_array()) {
    push @list, $row[0];
  }
  return @list;
}

# --------------------------------------------------------------------------------

# Set release for a source.

sub set_release
{
    my ($self, $source_id, $s_release ) = @_;

    my $sth =
      $dbi->prepare(
        "UPDATE source SET source_release=? WHERE source_id=?");

    print "Setting release to '$s_release' for source ID '$source_id'\n" if($verbose);

    $sth->execute( $s_release, $source_id );
    return;
}

sub get_dependent_mappings {
  my $self = shift;
  my $source_id = shift;

  my $sth =
    $dbi->prepare(
		  "select d.master_xref_id, d.dependent_xref_id from dependent_xref d, xref x where x.xref_id = d.dependent_xref_id and x.source_id = $source_id");

  $sth->execute();
  my $master_xref;
  my $dependent_xref;
  $sth->bind_columns(\$master_xref,\$dependent_xref);
  while($sth->fetch){
    $xref_dependent_mapped{$master_xref."|".$dependent_xref}=1;
  }
  $sth->finish;
  return;
}

sub get_ext_synonyms{
  my $self = shift;
  my $source_name = shift;
  my %ext_syns;
  my %seen;          # can be in more than once fro each type of external source.

  my $sql = 'SELECT  x.accession, x.label, sy.synonym FROM xref x, source so, synonym sy WHERE x.xref_id = sy.xref_id AND so.source_id = x.source_id AND so.name like "' . $source_name . '"';

  my $sth = $dbi->prepare($sql);

  $sth->execute;
  my ($acc, $label, $syn);
  $sth->bind_columns(\$acc, \$label, \$syn);

  my $count = 0;
  while($sth->fetch){
    if(!defined($seen{$acc.":".$syn})){
      push @{$ext_syns{$acc}}, $syn;
      push @{$ext_syns{$label}}, $syn;
      $count++;
    }
    $seen{$acc.":".$syn} = 1;
  }
  $sth->finish;

  return \%ext_syns;

}


#
# Store data needed to beable to revert to same stage as after parsing
#

sub parsing_finished_store_data {
  my $self = shift;

  # Store max id for

  # gene/transcript/translation_direct_xref     general_xref_id  #Does this change??

  # xref                                        xref_id
  # dependent_xref                              object_xref_id is all null
  # go_xref                                     object_xref_id
  # object_xref                                 object_xref_id
  # identity_xref                               object_xref_id

  my %table_and_key =
    ( 'xref' => "xref_id", 'object_xref' => "object_xref_id" );

  foreach my $table ( keys %table_and_key ) {
    my $sth = $dbi->prepare(
             "select MAX(" . $table_and_key{$table} . ") from $table" );
    $sth->execute;
    my $max_val;
    $sth->bind_columns( \$max_val );
    $sth->fetch;
    $sth->finish;
    $self->add_meta_pair( "PARSED_" . $table_and_key{$table},
                          $max_val || 1 );
  }
  return;
} ## end sub parsing_finished_store_data


1;

