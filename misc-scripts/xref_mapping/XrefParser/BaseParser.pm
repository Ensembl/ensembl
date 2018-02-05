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

package XrefParser::BaseParser;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception;
use XrefParser::FetchFiles;
use XrefParser::Database;
use Carp;
use DBI;
use Getopt::Long;

my $base_dir = File::Spec->curdir();

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
    croak 'No database specfied';
  }
  $verbose = $is_verbose;
  my $dbi = $database->dbi;

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

    my $io =undef;

    if(!(defined $file_name) or $file_name eq ''){
      confess "No file name";
    }
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
      print "Reading from '$file_name'...\n" || croak 'Could not print out message';
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
  my ($self, $source_name,$priority_desc, $dbi) = @_;
  $dbi = $self->dbi unless defined $dbi;

  my $low_name = lc $source_name;
  my $sql = "SELECT source_id FROM source WHERE LOWER(name)='$low_name'";
  if(defined $priority_desc){
    $low_name = lc $priority_desc;
    $sql .= " AND LOWER(priority_description)='$low_name'";
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

    $source_id = '-1';
  }
  $sth->finish();
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

  my ($self, $source_name, $dbi) = @_;
  $dbi = $self->dbi unless defined $dbi;

  my $big_name = uc $source_name;
  my $sql = "SELECT source_id FROM source WHERE upper(name) LIKE '%${big_name}%'";

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
  my ($self, $source_id, $dbi) = @_;
  $dbi = $self->dbi unless defined $dbi;
  my $source_name;

  my $sql = "SELECT name FROM source WHERE source_id= '$source_id'";
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
    $source_name = '-1';
  }
  $sth->finish;
  return $source_name;
}



####################################################
# Get a hash to go from accession of a dependent xref
# to master_xref_id for all of source names given
#####################################################
sub get_valid_xrefs_for_dependencies{
  my ($self, $dependent_name, $dbi, @reverse_ordered_source_list) = @_;
  $dbi = $self->dbi unless defined $dbi;

  my %dependent_2_xref;

  my $sql = 'select source_id from source where LOWER(name) =?';
  my $sth = $dbi->prepare($sql);
  my @dependent_sources;
  $sth->execute(lc $dependent_name);
  while(my @row = $sth->fetchrow_array()){
   push @dependent_sources,$row[0];
  }

  my @sources;
  foreach my $name (@reverse_ordered_source_list){
    $sth->execute(lc $name);
    while(my @row = $sth->fetchrow_array()){
      push @sources,$row[0];
    }
  }
  $sth->finish;

  my $dep_sql = (<<'DSS');
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
  $sth->finish;
  return \%dependent_2_xref;
}



####################################################
# Get a hash to go from accession of a direct xref
# to master_xref_id for all of source names given
#####################################################
sub get_valid_xrefs_for_direct_xrefs{
  my ($self, $direct_name, $separator, $dbi) = @_;
  $dbi = $self->dbi unless defined $dbi;

  my %direct_2_xref;

  my $sql = 'select source_id from source where name like ?';
  my $sth = $dbi->prepare($sql);
  my @direct_sources;
  $sth->execute("${direct_name}%");
  while(my @row = $sth->fetchrow_array()){
    push @direct_sources,$row[0];
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
    my $table = lc $type;
    $t_sql =~ s/TABLE/$table/xsm;
    $t_sql =~ s/TYPE/$type/xsm;

    $sth[$i++] = $dbi->prepare($t_sql);
  }

  foreach my $d (@direct_sources){
    for my $ii (0..2) {
      $sth[$ii]->execute($d);
      while(my ($gen_xref_id, $stable_id, $type, $link, $acc) = $sth[$ii]->fetchrow_array()){
	$direct_2_xref{$acc} = $gen_xref_id.$separator.$stable_id.$separator.$type.$separator.$link;
      }
      $sth[$ii]->finish();
    }
  }

  return \%direct_2_xref;
}


#############################################
# Get a hash of label to acc for a particular
# source name and species_id
#############################################
sub label_to_acc{

  my ($self,$source_name,$species_id, $dbi) =@_;
  $dbi = $self->dbi unless defined $dbi;

  # First cache synonyms so we can quickly add them later
  my %synonyms;
  my $syn_sth = $dbi->prepare('SELECT xref_id, synonym FROM synonym');
  $syn_sth->execute();

  my ($xref_id, $synonym);
  $syn_sth->bind_columns(\$xref_id, \$synonym);
  while ($syn_sth->fetch()) {

    push @{$synonyms{$xref_id}}, $synonym;

  }
  $syn_sth->finish;

  my %valid_codes;
  my @sources;

  my $big_name = uc $source_name;
  my $sql = "select source_id from source where upper(name) like '%${big_name}%'";
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
  $sth->finish;
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

  my ($self,$source_name,$species_id, $dbi) =@_;

  my %valid_codes;
  my @sources;
  $dbi = $self->dbi unless defined $dbi;

  my $big_name = uc $source_name;
  my $sql = "select source_id from source where upper(name) like '%$big_name%'";
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
    }
  }
  $sth->finish();
  return \%valid_codes;
}

##############################
# Upload xrefs to the database
##############################
sub upload_xref_object_graphs {
  my ($self, $rxrefs, $dbi) = @_;

  my $count = scalar @{$rxrefs};
  if($verbose) {
    print "count = $count\n" || croak 'Could not print out count';
  }

  if ($count) {

    #################
    # upload new ones
    ##################
    if ($verbose) {
      print "Uploading xrefs\n"
	|| croak 'Could not print string';
    }


    #################################################################################
    # Start of sql needed to add xrefs, primary_xrefs, synonym, dependent_xrefs etc..
    #################################################################################
    $dbi = $self->dbi unless defined $dbi;
    my $xref_sth = $dbi->prepare('INSERT INTO xref (accession,version,label,description,source_id,species_id, info_type) VALUES(?,?,?,?,?,?,?)');
    my $pri_insert_sth = $dbi->prepare('INSERT INTO primary_xref VALUES(?,?,?,?)');
    my $pri_update_sth = $dbi->prepare('UPDATE primary_xref SET sequence=? WHERE xref_id=?');
    my $syn_sth = $dbi->prepare('INSERT IGNORE INTO synonym VALUES(?,?)');
    my $dep_sth = $dbi->prepare('INSERT INTO dependent_xref (master_xref_id, dependent_xref_id, linkage_annotation, linkage_source_id) VALUES(?,?,?,?)');
    my $xref_update_label_sth = $dbi->prepare('UPDATE xref SET label=? WHERE xref_id=?');
    my $xref_update_descr_sth = $dbi->prepare('UPDATE xref SET description=? WHERE xref_id=?');
    my $pair_sth = $dbi->prepare('INSERT INTO pairs VALUES(?,?,?)');
    my $xref_id_sth = $dbi->prepare("SELECT xref_id FROM xref WHERE accession = ? AND source_id = ? AND species_id = ?");
    my $primary_xref_id_sth = $dbi->prepare('SELECT xref_id FROM primary_xref WHERE xref_id=?');



    # disable error handling here as we'll do it ourselves
    # reenabled it, as errorcodes are really unhelpful
    $xref_sth->{RaiseError} = 0;
    $xref_sth->{PrintError} = 0;

    #################################################################################
    # End of sql needed to add xrefs, primary_xrefs, synonym, dependent_xrefs etc..
    #################################################################################


    foreach my $xref (@{$rxrefs}) {
       my ($xref_id, $direct_xref_id);
       if(!(defined $xref->{ACCESSION} )){
	 print "Your xref does not have an accession-number,so it can't be stored in the database\n"
	   || croak 'Could not write message';
	 return;
       }

       ########################################
       # Create entry in xref table and note ID
       ########################################
       if(! $xref_sth->execute($xref->{ACCESSION},
			 $xref->{VERSION} || 0,
			 $xref->{LABEL}|| $xref->{ACCESSION},
			 $xref->{DESCRIPTION},
			 $xref->{SOURCE_ID},
			 $xref->{SPECIES_ID},
			 $xref->{INFO_TYPE} || 'MISC')){
	 #
	 #  if we failed to add the xref it must already exist so go find the xref_id for this
	 #
	 if(!(defined $xref->{SOURCE_ID})){
	   print "your xref: $xref->{ACCESSION} does not have a source-id\n";
	   return;
	 }
         $xref_id_sth->execute(
                   $xref->{ACCESSION},
                   $xref->{SOURCE_ID},
                   $xref->{SPECIES_ID} );
         $xref_id = ($xref_id_sth->fetchrow_array())[0];
	 if(defined $xref->{LABEL} ) {
	   $xref_update_label_sth->execute($xref->{LABEL},$xref_id) ;
	 }
	 if(defined $xref->{DESCRIPTION} ){
	   $xref_update_descr_sth->execute($xref->{DESCRIPTION},$xref_id);
	 }
       }
       else{
         $xref_id_sth->execute(
                   $xref->{ACCESSION},
                   $xref->{SOURCE_ID},
                   $xref->{SPECIES_ID} );
         $xref_id = ($xref_id_sth->fetchrow_array())[0];
       }

       foreach my $direct_xref (@{$xref->{DIRECT_XREFS}}) {
         $xref_sth->execute( $xref->{ACCESSION},
                             $xref->{VERSION} || 0,
                             $xref->{LABEL} || $xref->{ACCESSION},
                             $xref->{DESCRIPTION},
                             $direct_xref->{SOURCE_ID},
                             $xref->{SPECIES_ID},
                             $direct_xref->{LINKAGE_TYPE});
         $xref_id_sth->execute(
                   $xref->{ACCESSION},
                   $direct_xref->{SOURCE_ID},
                   $xref->{SPECIES_ID} );
         $direct_xref_id = ($xref_id_sth->fetchrow_array())[0];
         $self->add_direct_xref($direct_xref_id, $direct_xref->{STABLE_ID}, $direct_xref->{ENSEMBL_TYPE},$direct_xref->{LINKAGE_TYPE}, $dbi);
       }

       ################
       # Error checking
       ################
       if(!((defined $xref_id) and $xref_id)){
	 print STDERR "xref_id is not set for :\n".
	   "$xref->{ACCESSION}\n$xref->{LABEL}\n".
	     "$xref->{DESCRIPTION}\n$xref->{SOURCE_ID}\n".
	       "$xref->{SPECIES_ID}\n";
       }


       #############################################################################
       # create entry in primary_xref table with sequence; if this is a "cumulative"
       # entry it may already exist, and require an UPDATE rather than an INSERT
       #############################################################################
       if(defined $xref->{SEQUENCE} ){
         $primary_xref_id_sth->execute($xref_id) or croak( $dbi->errstr() );
         my @row = $primary_xref_id_sth->fetchrow_array();
         my $exists = $row[0];
	 if ( $exists ) {
	   $pri_update_sth->execute( $xref->{SEQUENCE}, $xref_id )
	     or croak( $dbi->errstr() );
	 } else {
	   $pri_insert_sth->execute( $xref_id, $xref->{SEQUENCE},
				     $xref->{SEQUENCE_TYPE},
				     $xref->{STATUS} )
	     or croak( $dbi->errstr() );
	 }
       }

       ##########################################################
       # if there are synonyms, add entries in the synonym table
       ##########################################################
       foreach my $syn ( @{ $xref->{SYNONYMS} } ) {
	 $syn_sth->execute( $xref_id, $syn )
	   or croak( $dbi->errstr() . "\n $xref_id\n $syn\n" );
       }

       #######################################################################
       # if there are dependent xrefs, add xrefs and dependent xrefs for them
       #######################################################################
       foreach my $depref (@{$xref->{DEPENDENT_XREFS}}) {
	 my %dep = %{$depref};

	 #################
	 # Insert the xref
	 #################
	 # print "inserting $dep{ACCESSION},$dep{VERSION},$dep{LABEL},$dep{DESCRIPTION},$dep{SOURCE_ID},${\$xref->{SPECIES_ID}}\n";
	 $xref_sth->execute($dep{ACCESSION},
			   $dep{VERSION} || 0,
			   $dep{LABEL} || $dep{ACCESSION},
			   $dep{DESCRIPTION} || '',
			   $dep{SOURCE_ID},
			   $xref->{SPECIES_ID},
                           'DEPENDENT');

	 #####################################
	 # find the xref_id for dependent xref
	 #####################################
	 $xref_id_sth->execute(
                   $dep{ACCESSION},
                   $dep{SOURCE_ID},
                   $xref->{SPECIES_ID} );
         my $dep_xref_id = ($xref_id_sth->fetchrow_array())[0];

	 if(!(defined $dep_xref_id) || $dep_xref_id ==0 ){
	   print STDERR "acc = $dep{ACCESSION} \nlink = $dep{LINKAGE_SOURCE_ID} \n".$dbi->err."\n";
	   print STDERR "source = $dep{SOURCE_ID}\n";
	 }

	 #
	 # Add the linkage_annotation and source id it came from
	 #
	 $dep_sth->execute( $xref_id, $dep_xref_id,
			    $dep{LINKAGE_ANNOTATION},
			    $dep{LINKAGE_SOURCE_ID} )
	   or croak( $dbi->errstr() );

	 #########################################################
	 # if there are synonyms, add entries in the synonym table
	 #########################################################
	 foreach my $syn ( @{ $dep{SYNONYMS} } ) {
	   $syn_sth->execute( $dep_xref_id, $syn )
	     or croak( $dbi->errstr() . "\n $xref_id\n $syn\n" );
	 } # foreach syn

       } # foreach dep

       #################################################
       # Add the pair data. refseq dna/pep pairs usually
       #################################################
       if(defined $xref_id and defined $xref->{PAIR} ){
	 $pair_sth->execute($xref->{SOURCE_ID},$xref->{ACCESSION},$xref->{PAIR});
       }


       ###########################
       # tidy up statement handles
       ###########################
       if(defined $xref_sth) {$xref_sth->finish()};
       if(defined $pri_insert_sth) {$pri_insert_sth->finish()} ;
       if(defined $pri_update_sth) {$pri_update_sth->finish()};
       if(defined $syn_sth) { $syn_sth->finish()};
       if(defined $dep_sth) { $dep_sth->finish()};
       if(defined $xref_update_label_sth) { $xref_update_label_sth->finish()};
       if(defined $xref_update_descr_sth) { $xref_update_descr_sth->finish()};
       if(defined $pair_sth) { $pair_sth->finish()};
       if(defined $xref_id_sth) { $xref_id_sth->finish()};
       if(defined $primary_xref_id_sth) { $primary_xref_id_sth->finish()};

     }  # foreach xref

  }
  return 1;
}

######################################################################################
# Add direct xref to the table XXX_direct_xref. (XXX -> Gene.Transcript or Translation
# Xref has to exist already, this module just adds ot yo the direct_xref table.
# $direct_xref is a reference to an array of hash objects.
######################################################################################
sub upload_direct_xrefs{
  my ($self, $direct_xref, $dbi)  = @_;
  $dbi = $self->dbi unless defined $dbi;
  for my $dr(@{$direct_xref}) {

    ################################################
    # Find the xref_id for this accession and source
    ################################################
    my $general_xref_id = $self->get_xref($dr->{ACCESSION},$dr->{SOURCE_ID},$dr->{SPECIES_ID}, $dbi);

    #######################################################
    # If found add the direct xref else write error message
    #######################################################
    if ($general_xref_id){
      $self->add_direct_xref($general_xref_id, $dr->{ENSEMBL_STABLE_ID},$dr->{ENSEMBL_TYPE},$dr->{LINKAGE_XREF}, $dbi);
    }
    else{
      print {*STDERR} 'Problem Could not find accession '.$dr->{ACCESSION}.' for source '.$dr->{SOURCE}.
	' so not able to add direct xref to '.$dr->{ENSEMBL_STABLE_ID}."\n";
    }
  }
  return;
}




###############################################
# Insert into the meta table the key and value.
###############################################
sub add_meta_pair {

  my ($self, $key, $value, $dbi) = @_;
  $dbi = $self->dbi unless defined $dbi;

  my $sth = $dbi->prepare('insert into meta (meta_key, meta_value, date) values("'.$key.'", "'.$value.'", now())');
  $sth->execute;
  $sth->finish;
  return;
}



#################################################
# Create a hash of all the source names for xrefs
#################################################
sub get_xref_sources {

  my $self = shift;
  my $dbi = shift;
  $dbi = $self->dbi unless defined $dbi;
  my %sourcename_to_sourceid;

  my $sth = $dbi->prepare('SELECT name,source_id FROM source');
  $sth->execute() or croak( $dbi->errstr() );
  while(my @row = $sth->fetchrow_array()) {
    my $source_name = $row[0];
    my $source_id = $row[1];
    $sourcename_to_sourceid{$source_name} = $source_id;
  }
  $sth->finish;

  return %sourcename_to_sourceid;
}

########################################################################
# Create and return a hash that that goes from species_id to taxonomy_id
########################################################################
sub species_id2taxonomy {

  my $self = shift;
  my $dbi = shift;
  $dbi = $self->dbi unless defined $dbi;

  my %species_id2taxonomy;

  my $sth = $dbi->prepare('SELECT species_id, taxonomy_id FROM species');
  $sth->execute() or croak( $dbi->errstr() );
  while(my @row = $sth->fetchrow_array()) {
    my $species_id = $row[0];
    my $taxonomy_id = $row[1];
    if(defined $species_id2taxonomy{$species_id} ){
      push @{$species_id2taxonomy{$species_id}}, $taxonomy_id;
    }
    else{
      $species_id2taxonomy{$species_id} = [$taxonomy_id];
    }
  }
  $sth->finish();
  return %species_id2taxonomy;
}



#########################################################################
# Create and return a hash that that goes from species_id to species name
#########################################################################
sub species_id2name {
  my $self = shift;
  my $dbi = shift;
  $dbi = $self->dbi unless defined $dbi;

  my %species_id2name;

  my $sth = $dbi->prepare('SELECT species_id, name FROM species');
  $sth->execute() or croak( $dbi->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
    my $species_id = $row[0];
    my $name       = $row[1];
    $species_id2name{$species_id} = [ $name ];
  }
  $sth->finish();

  ##############################################
  # Also populate the hash with all the aliases.
  ##############################################
  $sth = $dbi->prepare('SELECT species_id, aliases FROM species');
  $sth->execute() or croak( $dbi->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
    my $species_id = $row[0];
    foreach my $name ( split /,\s*/xms, $row[1] ) {
      $species_id2name{$species_id} ||= [];
      push @{$species_id2name{$species_id}}, $name;
    }
  }
  $sth->finish();

  return %species_id2name;
} ## end sub species_id2name


###########################################################################
# If there was an error, an xref with the same acc & source already exists.
# If so, find its ID, otherwise get ID of xref just inserted
###########################################################################
sub get_xref_id {
  my ($self, $arg_ref) = @_;
  my $sth     = $arg_ref->{sth}        || croak 'Need a statement handle for get_xref_id';
  my $acc     = $arg_ref->{acc}        || croak 'Need an accession for get_xref_id';
  my $source  = $arg_ref->{source_id}  || croak 'Need an source_id for get_xref_id';
  my $species = $arg_ref->{species_id} || confess 'Need an species_id for get_xref_id';
  my $error   = $arg_ref->{error};
  my $dbi     = $arg_ref->{dbi};

  my $id = $self->get_xref($acc, $source, $species, $dbi);

  return $id;
}

##################################################################
# If primary xref already exists for a partiuclar xref_id return 1
# else return 0;
##################################################################
sub primary_xref_id_exists {

  my ($self, $xref_id, $dbi) = @_;
  $dbi = $self->dbi unless defined $dbi;

  my $exists = 0;

  my $sth = $dbi->prepare('SELECT xref_id FROM primary_xref WHERE xref_id=?');
  $sth->execute($xref_id) or croak( $dbi->errstr() );
  my @row = $sth->fetchrow_array();
  my $result = $row[0];
  if (defined $result) {$exists = 1; }
  $sth->finish();

  return $exists;

}

############################################
# Get the tax id for a particular species id
############################################
sub get_taxonomy_from_species_id{
  my ($self,$species_id, $dbi) = @_;
  my %hash;

  $dbi = $self->dbi unless defined $dbi;
  my $sth = $dbi->prepare("SELECT taxonomy_id FROM species WHERE species_id = $species_id");
  $sth->execute() or croak( $dbi->errstr() );
  while(my @row = $sth->fetchrow_array()) {
    $hash{$row[0]} = 1;
  }
  $sth->finish;
  return \%hash;
}


################################################
# xref_id for a given stable id and linkage_xref
# Only used in GOParser at the moment
################################################
sub get_direct_xref{
 my ($self,$stable_id,$type,$link, $dbi) = @_;
 $dbi = $self->dbi unless defined $dbi;

 $type = lc $type;

 my $sql = "select general_xref_id from ${type}_direct_xref d where ensembl_stable_id = ?  and linkage_xref= ?";
 my  $direct_sth = $dbi->prepare($sql);

 $direct_sth->execute( $stable_id, $link ) or croak( $dbi->errstr() );
 if(my @row = $direct_sth->fetchrow_array()) {
   return $row[0];
 }
 $direct_sth->finish();
 return;
}


###################################################################
# return the xref_id for a particular accession, source and species
# if not found return undef;
###################################################################
sub get_xref{
  my ($self,$acc,$source, $species_id, $dbi) = @_;
  $dbi = $self->dbi unless defined $dbi;

  #
  # If the statement handle does nt exist create it.
  #
  my $sql = 'select xref_id from xref where accession = ? and source_id = ? and species_id = ?';
  my $get_xref_sth = $dbi->prepare($sql);

  #
  # Find the xref_id using the sql above
  #
  $get_xref_sth->execute( $acc, $source, $species_id ) or croak( $dbi->errstr() );
  if(my @row = $get_xref_sth->fetchrow_array()) {
    return $row[0];
  }
  $get_xref_sth->finish();
  return;
}

###################################################################
# return the object_xref_id for a particular xref_id, ensembl_id and ensembl_object_type
# if not found return undef;
###################################################################
sub get_object_xref {
  my ($self, $xref_id, $ensembl_id, $object_type, $dbi) = @_;
  $dbi = $self->dbi unless defined $dbi;

  my $sql = 'select object_xref_id from object_xref where xref_id = ? and ensembl_object_type = ? and ensembl_id = ?';
  my $get_object_xref_sth = $dbi->prepare($sql);

  #
  # Find the object_xref_id using the sql above
  #
  $get_object_xref_sth->execute( $xref_id, $ensembl_id, $object_type) or croak( $dbi->errstr() );
  if(my @row = $get_object_xref_sth->fetchrow_array()) {
    return $row[0];
  }
  $get_object_xref_sth->finish();
  return;
}


###########################################################
# Create an xref..
# If it already exists it return that xrefs xref_id
# else creates it and return the new xre_id
###########################################################
sub add_xref {
  my ( $self, $arg_ref) = @_;

  my $acc         = $arg_ref->{acc}        || croak 'add_xref needs aa acc';
  my $source_id   = $arg_ref->{source_id}  || croak 'add_xref needs a source_id';
  my $species_id  = $arg_ref->{species_id} || croak 'add_xref needs a species_id';
  my $label       = $arg_ref->{label}      || $acc;
  my $description = $arg_ref->{desc}       || '';
  my $version     = $arg_ref->{version}    || 0;
  my $info_type   = $arg_ref->{info_type}  || 'MISC';
  my $info_text   = $arg_ref->{info_text}  || '';
  my $dbi         = $arg_ref->{dbi};

  $dbi = $self->dbi unless defined $dbi;

  ##################################################################
  # See if it already exists. It so return the xref_id for this one.
  ##################################################################
  my $xref_id = $self->get_xref($acc,$source_id, $species_id, $dbi);
  if(defined $xref_id){
    return $xref_id;
  }

  my $add_xref_sth =
      $dbi->prepare( 'INSERT INTO xref '
         . '(accession,version,label,description,source_id,species_id, info_type, info_text) '
         . 'VALUES(?,?,?,?,?,?,?,?)' );

  ######################################################################
  # If the description is more than 255 characters, chop it off and add
  # an indication that it has been truncated to the end of it.
  ######################################################################
  if (defined $description && ((length $description) > 255 ) ) {
    my $truncmsg = ' /.../';
    substr $description, 255 - (length $truncmsg),
            length $truncmsg, $truncmsg;
  }


  ####################################
  # Add the xref and croak if it fails
  ####################################
  $add_xref_sth->execute( $acc, $version || 0, $label,
                          $description, $source_id, $species_id, $info_type, $info_text
  ) or croak("$acc\t$label\t\t$source_id\t$species_id\n");

  $add_xref_sth->finish();
  return $add_xref_sth->{'mysql_insertid'};
} ## end sub add_xref


###########################################################
# Create an object_xref..
# If it already exists it return the object_xref_id
# else creates it and returns the new object_xref_id
###########################################################

sub add_object_xref {
  my ($self, $arg_ref) = @_;

  my $xref_id     = $arg_ref->{xref_id}     || croak 'add_object_xref needs an xref_id';
  my $ensembl_id  = $arg_ref->{ensembl_id}  || croak 'add_object_xref needs a ensembl_id';
  my $object_type = $arg_ref->{object_type} || croak 'add_object_xref needs an object_type';
  my $dbi         = $arg_ref->{dbi};

  $dbi = $self->dbi unless defined $dbi;

  ##################################################################
  # See if it already exists. It so return the xref_id for this one.
  ##################################################################

  my $object_xref_id = $self->get_object_xref($xref_id, $ensembl_id, $object_type, $dbi);
  if(defined $object_xref_id){
    return $object_xref_id;
  }

  my $add_object_xref_sth =
      $dbi->prepare( 'INSERT INTO object_xref'
         . '(ensembl_id, ensembl_object_type, xref_id) '
         . 'VALUES(?,?,?)' );

  ####################################
  # Add the object_xref and croak if it fails
  ####################################
  $add_object_xref_sth->execute($ensembl_id, $object_type, $xref_id
  ) or croak("$ensembl_id\t$object_type\t\t$xref_id\n");

  $add_object_xref_sth->finish();
  return $add_object_xref_sth->{'mysql_insertid'};
}

###########################################################
# Create an identity_xref
###########################################################

sub add_identity_xref {
  my ($self, $arg_ref) = @_;

  my $object_xref_id  = $arg_ref->{object_xref_id}  || croak 'add_identity_xref needs an object_xref_id';
  my $score           = $arg_ref->{score}           || croak 'add_identity_xref needs a score';
  my $target_identity = $arg_ref->{target_identity} || croak 'add_identity_xref needs a target_identity';
  my $query_identity  = $arg_ref->{query_identity}  || croak 'add_identity_xref needs a query_identity';
  my $dbi             = $arg_ref->{dbi};

  $dbi = $self->dbi unless defined $dbi;

  my $add_identity_xref_sth =
      $dbi->prepare( 'INSERT INTO identity_xref'
         . '(object_xref_id, score, query_identity, target_identity) '
         . 'VALUES(?,?,?,?)' );

  ####################################
  # Add the object_xref and croak if it fails
  ####################################
  $add_identity_xref_sth->execute($object_xref_id, $score, $query_identity, $target_identity
  ) or croak("$object_xref_id\t$score\t\t$query_identity\t$target_identity\n");
  $add_identity_xref_sth->finish();
  return;
}


###################################################################
# Create new xref if needed and add as a direct xref to a stable_id
###################################################################
sub add_to_direct_xrefs{
  my ($self, $arg_ref) = @_;

  my $stable_id   = $arg_ref->{stable_id}   || croak ('Need a direct_xref on which this xref linked too' );
  my $type        = $arg_ref->{type}        || croak ('Need a table type on which to add');
  my $acc         = $arg_ref->{acc}         || croak ('Need an accession of this direct xref' );
  my $source_id   = $arg_ref->{source_id}   || croak ('Need a source_id for this direct xref' );
  my $species_id  = $arg_ref->{species_id}  || croak ('Need a species_id for this direct xref' );
  my $version     = $arg_ref->{version}     || 0;
  my $label       = $arg_ref->{label}       || $acc;
  my $description = $arg_ref->{desc};
  my $linkage     = $arg_ref->{linkage};
  my $info_text   = $arg_ref->{info_text}  || '';
  my $dbi         = $arg_ref->{dbi};

  $dbi = $self->dbi unless defined $dbi;

  my $sql = (<<'AXX');
  INSERT INTO xref (accession,version,label,description,source_id,species_id, info_type, info_text)
          VALUES (?,?,?,?,?,?,?,?)
AXX
  my $add_xref_sth = $dbi->prepare($sql);

  ###############################################################
  # If the acc already has an xrefs find it else cretae a new one
  ###############################################################
  my $direct_id = $self->get_xref($acc, $source_id, $species_id, $dbi);
  if(!(defined $direct_id)){
    $add_xref_sth->execute(
        $acc, $version || 0, $label,
        $description, $source_id, $species_id, 'DIRECT', $info_text
    ) or croak("$acc\t$label\t\t$source_id\t$species_id\n");
  }
  $add_xref_sth->finish();

  $direct_id = $self->get_xref($acc, $source_id, $species_id, $dbi);

  #########################
  # Now add the direct info
  #########################
  $self->add_direct_xref($direct_id, $stable_id, $type, '', $dbi);
  return;
}


##################################################################
# Add a single record to the direct_xref table.
# Note that an xref must already have been added to the xref table
##################################################################
sub add_direct_xref {
  my ($self, $general_xref_id, $ensembl_stable_id, $ensembl_type, $linkage_type, $dbi) = @_;

  $dbi = $self->dbi unless defined $dbi;

  $ensembl_type = lc($ensembl_type);
  my $sql = "INSERT INTO " . $ensembl_type . "_direct_xref VALUES (?,?,?)";
  my $add_direct_xref_sth = $dbi->prepare($sql);

  $add_direct_xref_sth->execute($general_xref_id, $ensembl_stable_id, $linkage_type);
  $add_direct_xref_sth->finish();
  return;
}


##########################################################
# Create/Add xref and add it as a dependency of the master
##########################################################
sub add_dependent_xref{
  my ($self, $arg_ref) = @_;

  my $master_xref = $arg_ref->{master_xref_id} || croak( 'Need a master_xref_id on which this xref depends on' );
  my $acc         = $arg_ref->{acc}            || croak( 'Need an accession of this dependent xref' );
  my $source_id   = $arg_ref->{source_id}      || croak( 'Need a source_id for this dependent xref' );
  my $species_id  = $arg_ref->{species_id}     || croak( 'Need a species_id for this dependent xref' );
  my $version     = $arg_ref->{version}        ||  0;
  my $label       = $arg_ref->{label}          || $acc;
  my $description = $arg_ref->{desc};
  my $linkage     = $arg_ref->{linkage};
  my $info_text   = $arg_ref->{info_text} || '';
  my $dbi         = $arg_ref->{dbi};

  $dbi = $self->dbi unless defined $dbi;

  my $sql = (<<'IXR');
INSERT INTO xref
  (accession,version,label,description,source_id,species_id, info_type, info_text)
  VALUES (?,?,?,?,?,?,?,?)
IXR
  my $add_xref_sth = $dbi->prepare($sql);
  $sql = (<<'ADX');
INSERT INTO dependent_xref 
  (master_xref_id,dependent_xref_id,linkage_annotation,linkage_source_id)
  VALUES (?,?,?,?)
ADX
  my $add_dependent_xref_sth = $dbi->prepare($sql);

  ####################################################
  # Does the xref already exist. If so get its xref_id
  # else create it and get the new xref_id
  ####################################################
  my $dependent_id = $self->get_xref($acc, $source_id, $species_id, $dbi);
  if(!(defined $dependent_id)){
    $add_xref_sth->execute(
        $acc, $version, $label,
        $description, $source_id, $species_id, 'DEPENDENT', $info_text
    ) or croak("$acc\t$label\t\t$source_id\t$species_id\n");
  }
  $add_xref_sth->finish();
  $dependent_id = $self->get_xref($acc, $source_id, $species_id, $dbi);

  ################################################
  # Croak if we have failed to create.get the xref
  ################################################
  if(!(defined $dependent_id)){
    croak("$acc\t$label\t\t$source_id\t$species_id\n");
  }

  ########################################################################################
  # If the dependency has not already been set ( is already in hash xref_dependent_mapped)
  # then add it
  ########################################################################################
  if(!(defined $xref_dependent_mapped{"$master_xref|$dependent_id"}) || $xref_dependent_mapped{"$master_xref|$dependent_id"} ne $linkage){
    $add_dependent_xref_sth->execute( $master_xref, $dependent_id, $linkage,
				      $source_id )
      or croak("$master_xref\t$dependent_id\t$linkage\t$source_id");
    $xref_dependent_mapped{"$master_xref|$dependent_id"} = $linkage;
  }
  $add_dependent_xref_sth->finish();

  return $dependent_id;
}

##################################################################
# Add synonyms for a particular accession for one or more sources.
# This is for priority xrefs where we have more than one source
# but want to write synonyms for each with the same accession
##################################################################
sub add_to_syn_for_mult_sources{
  my ($self, $acc, $sources, $syn, $species_id, $dbi) = @_;

  $dbi = $self->dbi unless defined $dbi;
  my $add_synonym_sth =  $dbi->prepare('INSERT IGNORE INTO synonym VALUES(?,?)');

  foreach my $source_id (@{$sources}){
    my $xref_id = $self->get_xref($acc, $source_id, $species_id, $dbi);
    if(defined $xref_id){
      $add_synonym_sth->execute( $xref_id, $syn )
        or croak( $dbi->errstr() . "\n $xref_id\n $syn\n" );
    }
  }
  $add_synonym_sth->finish();
  return;
}


##########################################################
# Add synomyn for an xref given by accession and source_id
##########################################################
sub add_to_syn{
  my ($self, $acc, $source_id, $syn, $species_id, $dbi) = @_;

  $dbi = $self->dbi unless defined $dbi;
  my $add_synonym_sth =  $dbi->prepare('INSERT IGNORE INTO synonym VALUES(?,?)');
  my $xref_id = $self->get_xref($acc, $source_id, $species_id, $dbi);
  if(defined $xref_id){
    $add_synonym_sth->execute( $xref_id, $syn )
      or croak( $dbi->errstr() . "\n $xref_id\n $syn\n" );
  }
  else {
      carp (  "Could not find acc $acc in "
            . "xref table source = $source_id of species $species_id\n" );
  }
  $add_synonym_sth->finish();
  return;
}


##########################################
# Add synomyn for an xref given by xref_id
##########################################
sub add_synonym{
  my ($self, $xref_id, $syn, $dbi) = @_;

  $dbi = $self->dbi unless defined $dbi;
  my $add_synonym_sth =  $dbi->prepare('INSERT IGNORE INTO synonym VALUES(?,?)');
  $add_synonym_sth->execute( $xref_id, $syn ) 
    or croak( $dbi->errstr()."\n $xref_id\n $syn\n\n" );

  $add_synonym_sth->finish();
  return;
}

########################################################
# Create a hash that uses the label as a key
# and the acc as the value. Also add synonyms for these
# as keys.
#######################################################
sub get_label_to_acc{
  my ($self, $name, $species_id, $prio_desc, $dbi) = @_;
  my %hash1=();

  $dbi = $self->dbi unless defined $dbi;

  my $sql =(<<"GLA");
SELECT  xref.accession, xref.label
  FROM xref, source
    WHERE source.name LIKE '$name%' AND
          xref.source_id = source.source_id
GLA
  if(defined $prio_desc){
    $sql .= " and source.priority_description like '$prio_desc'";
  }
  if(defined $species_id){
    $sql .= " and xref.species_id  = $species_id";
  }
  my $sub_sth = $dbi->prepare($sql);

  $sub_sth->execute();
  while(my @row = $sub_sth->fetchrow_array()) {
    $hash1{$row[1]} = $row[0];
  }


  ####################
  # Remember synonyms
  ####################

 $sql =(<<"GLS");
SELECT  xref.accession, synonym.synonym 
  FROM xref, source, synonym 
    WHERE synonym.xref_id = xref.xref_id AND
          source.name like '$name%' AND
           xref.source_id = source.source_id
GLS

  if(defined $prio_desc){
    $sql .= " AND source.priority_description LIKE '$prio_desc'";
  }
  if(defined $species_id){
    $sql .= " AND xref.species_id  = $species_id";
  }
  $sub_sth = $dbi->prepare($sql);

  $sub_sth->execute();
  while(my @row = $sub_sth->fetchrow_array()) {
    $hash1{$row[1]} = $row[0];
  }
  $sub_sth->finish();

  return \%hash1;
}

########################################################
# Create a hash that uses the accession as a key
# and the label as the value.
#######################################################
sub get_acc_to_label{
  my ($self, $name, $species_id, $prio_desc, $dbi) = @_;
  my %hash1=();

  $dbi = $self->dbi unless defined $dbi;

  my $sql =(<<"GLA");
SELECT  xref.accession, xref.label
  FROM xref, source
    WHERE source.name LIKE '$name%' AND
          xref.source_id = source.source_id
GLA
  if(defined $prio_desc){
    $sql .= " and source.priority_description like '$prio_desc'";
  }
  if(defined $species_id){
    $sql .= " and xref.species_id  = $species_id";
  }
  my $sub_sth = $dbi->prepare($sql);

  $sub_sth->execute();
  while(my @row = $sub_sth->fetchrow_array()) {
    $hash1{$row[0]} = $row[1];
  }
  $sub_sth->finish();

  return \%hash1;
}


########################################################
# Create a hash that uses the label as a key
# and the desc as the value. Also add synonyms for these
# as keys.
#######################################################
sub get_label_to_desc{
  my ($self, $name, $species_id, $prio_desc, $dbi) = @_;
  my %hash1=();

  $dbi = $self->dbi unless defined $dbi;

  my $sql =(<<"GDH");
  SELECT xref.description, xref.label 
    FROM xref, source 
      WHERE source.name LIKE '$name%' AND 
            xref.source_id = source.source_id
GDH
  if(defined $prio_desc){
    $sql .= " and source.priority_description like '$prio_desc'";
  }
  if(defined $species_id){
    $sql .= " and xref.species_id  = $species_id";
  }
  my $sub_sth = $dbi->prepare($sql);

  $sub_sth->execute();
  while(my @row = $sub_sth->fetchrow_array()) {
    $hash1{$row[1]} = $row[0];
  }

  ###########################
  # Also include the synonyms
  ###########################

  my $syn_sql =(<<"GDS");
  SELECT xref.description, synonym.synonym
    FROM xref, source, synonym
      WHERE synonym.xref_id = xref.xref_id AND
            source.name like '$name%' AND
             xref.source_id = source.source_id
GDS

  if(defined $prio_desc){
    $syn_sql .= " AND source.priority_description LIKE '$prio_desc'";
  }
  if(defined $species_id){
    $syn_sql .= " AND xref.species_id  = $species_id";
  }
  $sub_sth = $dbi->prepare($syn_sql);

  $sub_sth->execute();
  while(my @row = $sub_sth->fetchrow_array()) {
    $hash1{$row[1]} = $row[0];
  }
  $sub_sth->finish();

  return \%hash1;
}


########################################
# Set release for a particular source_id.
########################################
sub set_release{
    my ($self, $source_id, $s_release, $dbi ) = @_;

   $dbi = $self->dbi unless defined $dbi;

    my $sth =
      $dbi->prepare('UPDATE source SET source_release=? WHERE source_id=?');

    if($verbose) { print "Setting release to '$s_release' for source ID '$source_id'\n"; }

    $sth->execute( $s_release, $source_id );
    $sth->finish();
    return;
}


#############################################################################
# create a hash of all the dependent mapping that exist for a given source_id
# Of the format {master_xref_id|dependent_xref_id}
#############################################################################
sub get_dependent_mappings {
  my $self = shift;
  my $source_id = shift;
  my $dbi = shift;

  $dbi = $self->dbi unless defined $dbi;

  my $sql =(<<"GDM");
  SELECT  d.master_xref_id, d.dependent_xref_id
    FROM dependent_xref d, xref x
      WHERE x.xref_id = d.dependent_xref_id AND
            x.source_id = $source_id
GDM
  my $sth = $dbi->prepare($sql);
  $sth->execute();
  my $master_xref;
  my $dependent_xref;
  $sth->bind_columns(\$master_xref,\$dependent_xref);
  while($sth->fetch){
    $xref_dependent_mapped{"$master_xref|$dependent_xref"}=1;
  }
  $sth->finish;
  return;
}


##########################################################
# Create a has that uses the accession and labels for keys
# and an array of the synonyms as the vaules
##########################################################
sub get_ext_synonyms{
  my $self = shift;
  my $source_name = shift;
  my $dbi = shift;
  $dbi = $self->dbi unless defined $dbi;
  my %ext_syns;
  my %seen;          # can be in more than once fro each type of external source.
  my $separator = qw{:};

  my $sql =(<<"GES");
  SELECT  x.accession, x.label, sy.synonym
    FROM xref x, source so, synonym sy
      WHERE x.xref_id = sy.xref_id AND
            so.source_id = x.source_id AND
            so.name like '$source_name'
GES
  my $sth = $dbi->prepare($sql);

  $sth->execute;
  my ($acc, $label, $syn);
  $sth->bind_columns(\$acc, \$label, \$syn);

  my $count = 0;
  while($sth->fetch){
    if(!(defined $seen{$acc.$separator.$syn})){
      push @{$ext_syns{$acc}}, $syn;
      push @{$ext_syns{$label}}, $syn;
      $count++;
    }
    $seen{$acc.$separator.$syn} = 1;
  }
  $sth->finish;

  return \%ext_syns;

}


######################################################################
# Store data needed to beable to revert to same stage as after parsing
######################################################################
sub parsing_finished_store_data {
  my $self = shift;
  my $dbi = shift;
  $dbi = $self->dbi unless defined $dbi;

  # Store max id for

  # gene/transcript/translation_direct_xref     general_xref_id  #Does this change??

  # xref                                        xref_id
  # dependent_xref                              object_xref_id is all null
  # go_xref                                     object_xref_id
  # object_xref                                 object_xref_id
  # identity_xref                               object_xref_id

  my %table_and_key =
    ( 'xref' => 'xref_id', 'object_xref' => 'object_xref_id' );

  foreach my $table ( keys %table_and_key ) {
    my $sth = $dbi->prepare(
             'select MAX(' . $table_and_key{$table} . ") from $table" );
    $sth->execute;
    my $max_val;
    $sth->bind_columns( \$max_val );
    $sth->fetch;
    $sth->finish;
    $self->add_meta_pair( 'PARSED_' . $table_and_key{$table},
                          $max_val || 1, $dbi );
  }
  return;
} ## end sub parsing_finished_store_data


sub get_meta_value {
  my ($self, $key, $dbi) = @_;
  $dbi = $self->dbi unless defined $dbi;

  my $sth = $dbi->prepare('select meta_value from meta where meta_key like "'.$key.'" order by meta_id');
  $sth->execute();
  my $value;
  $sth->bind_columns(\$value);
  while($sth->fetch){   # get the last one
  }
  $sth->finish;

  return $value;
}



1;

