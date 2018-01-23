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

package XrefParser::ProcessData;

use strict;
use warnings;
use Carp;
use XrefParser::Database;

use File::Basename;
use File::Spec::Functions;
use IO::File;
use Digest::MD5 qw(md5_hex);

sub new
{
    my ($proto) = @_;

    my $class = ref $proto || $proto;
    return bless {}, $class;
}


##############################################################
# Main routine (run).
#   1) creates xref database if new one needed
#   2) process the downloadable sources;
#      a) fetch files if needed
#      b) parse and store xrefs/object xrefs etc
#      c) set checksum so that we know they have been processed
###############################################################
sub run {
  my ($self, $ref_arg) = @_;

  my $base_dir   = $ref_arg->{dl_path};
  my $release    = $ref_arg->{release};
  my $verbose    = $ref_arg->{verbose};
  my $unzip      = $ref_arg->{unzip};
  my $stats      = $ref_arg->{stats};
  my $cleanup    = $ref_arg->{cleanup};
  my $rspecies   = $ref_arg->{speciesr};
  my $sources    = $ref_arg->{sourcesr};
  my $notsources = $ref_arg->{notsourcesr};

  my $sql_dir = dirname($0);

  my $dbc = XrefParser::Database->new({ host    => $ref_arg->{host},
					dbname  => $ref_arg->{dbname},
					port    => $ref_arg->{port},
					user    => $ref_arg->{user},
					pass    => $ref_arg->{pass},
					verbose => $ref_arg->{verbose} });

  $self->database($dbc);

  if ($ref_arg->{create}) {
    $dbc->create($sql_dir,
		 $ref_arg->{force},
		 $ref_arg->{drop_db});
  }

  my $dbi = $dbc->dbi();
  $self->dbi($dbi);

  my $sth_c = $dbi->prepare("insert into process_status (status, date) values('xref_created',now())");
  $sth_c->execute;


  # need to use name now and that download = 'Y' as well

  my $sql = (<<"DSS");
SELECT s.name, s.source_id 
  FROM source s, dependent_source ds, source_url su
    WHERE su.source_id = s.source_id AND
          ds.dependent_name = s.name AND
          ds.master_source_id = ? AND
          su.species_id = ? AND
          s.download = 'Y' AND
          s.source_id != ? AND
          su.checksum is null
DSS
  my $dep_sth = $dbi->prepare($sql);

  # validate species names
  my @species_ids = $self->validate_species($rspecies, $verbose);

  # validate source names
  exit(1) if ( !$self->validate_sources(\@species_ids,$sources, $verbose) );
  exit(1) if ( !$self->validate_sources(\@species_ids,$notsources, $verbose) );

  # build SQL
  my $species_sql = "";
  if (@species_ids) {
    $species_sql .= " AND su.species_id IN (";
    for ( my $i = 0 ; $i < @species_ids ; $i++ ) {
      $species_sql .= "," if ( $i != 0 );
      $species_sql .= $species_ids[$i];
    }
    $species_sql .= ") ";
  }

  my $source_sql = "";
  if (defined $sources and @$sources) {
    $source_sql .= " AND LOWER(s.name) IN (";
    for ( my $i = 0 ; $i < @$sources ; $i++ ) {
      $source_sql .= "," if ( $i != 0 );
      $source_sql .= "\'" . lc( $$sources[$i] ) . "\'";
    }
    $source_sql .= ") ";
  }

  if (defined $notsources and @$notsources) {
    $source_sql .= " AND LOWER(s.name) NOT IN (";
    for ( my $i = 0 ; $i < @$notsources ; $i++ ) {
      $source_sql .= "," if ( $i != 0 );
      $source_sql .= "\'" . lc( $$notsources[$i] ) . "\'";
    }
    $source_sql .= ") ";
  }
  my $sth = $dbi->prepare("insert into process_status (status, date) values('parsing_started',now())");
  $sth->execute;

  $sql =
      "SELECT DISTINCT(s.source_id), su.source_url_id, s.name, su.url, "
      . "su.release_url, su.checksum, su.parser, su.species_id "
      . "FROM source s, source_url su, species sp "
      . "WHERE s.download='Y' AND su.source_id=s.source_id "
      . "AND su.species_id=sp.species_id "
      . $source_sql
      . $species_sql
      . "ORDER BY s.ordered";
#  print {*STDERR} $sql . "\n";

  $sth = $dbi->prepare($sql);
  $sth->execute();

  my ( $source_id, $source_url_id, $name, $url, $release_url,
       $checksum, $parser, $species_id );

    $sth->bind_columns( \$source_id,   \$source_url_id,
                        \$name,        \$url,
                        \$release_url, \$checksum,
                        \$parser,      \$species_id );

  my $dir;
  my %summary = ();

  my %sum_xrefs;
  my %sum_prim;
  my %sum_dep;
  my %sum_dir;
  my %sum_coord;
  my %sum_syn;

  #
  # We may be continuing on form a previous run  so find what xreefs we have already
  # there so that when we calculate the number of new xrefs etc we start with the correct number
  #
  $self->get_stats({ xrefs   => \%sum_xrefs,
		     primary => \%sum_prim,
		     depend  => \%sum_dep,
		     direct  => \%sum_dir,
		     coord   => \%sum_coord,
		     synonym => \%sum_syn,
		     dbi     => $dbi,
		     print   => 0 });


  my $start_transaction_sth = $dbi->prepare("start transaction");

  my $end_transaction_sth = $dbi->prepare("commit");


  while ( my @row = $sth->fetchrow_array() ) {
    print '-' x 4, "{ $name }", '-' x ( 72 - length($name) ), "\n" if ($verbose);

    my $cs;
    my $file_cs = "";
    my $parse   = 0;
    my $empty   = 0;
    my $type    = $name;
    my $dsn;
    
    my @files = split( /\s+/x, $url );
    my @files_to_parse = ();

    $dir = catdir( $base_dir, sanitise($type) );

    # For summary purposes: If 0 is returned (in
    # $summary{$name}->{$parser}) then it is successful.  If 1 is
    # returned then it failed.  If undef/nothing is returned the we
    # do not know.
    $summary{$name}->{$parser} = 0;

    my $ff=  XrefParser::FetchFiles->new();
    @files = $ff->fetch_files( {dest_dir  =>  $dir,
				user_uris => \@files,
				del_down  => $ref_arg->{deletedownloaded},
				chk_down  => $ref_arg->{checkdownload},
				verbose   => $verbose
			       });
    if ( !@files ) {
      # Fetching failed.
      ++$summary{$name}->{$parser};
      next;
    }
    if ( defined($release_url) ) {
      my @rel=();
      push @rel , $release_url;
      @rel = $ff->fetch_files( {dest_dir  =>  $dir,
				user_uris => \@rel,
				del_down  => $ref_arg->{deletedownloaded},
				chk_down  => $ref_arg->{checkdownload},
				verbose   => $verbose}
			     );
      $release_url = $rel[-1];
    }
    $start_transaction_sth->execute();

    foreach my $file (@files) {
	
      # check dependencies are loaded all ready
      if(!($self->all_dependencies_loaded($source_id, $species_id, $name, $dep_sth))){
	++$summary{$name}->{$parser};
	next;
      }
      # Database parsing
      if ( $file =~ /^mysql:/ix ) {
	$dsn = $file;
	print "Parsing $dsn with $parser\n" if ($verbose);
	my $eval_test = eval "require XrefParser::$parser";
	if($@ or $eval_test != 1) {
	  croak "Could not require XrefParser::$parser\ndollar=at=$@\neval_test = $eval_test\n";
	}
	my $new = "XrefParser::$parser"->new($self->database, $verbose);
	if (
	    $new->run( { dsn        => $dsn,
			 source_id  => $source_id,
			 species_id => $species_id,
			 name       => $name,
			 verbose    => $verbose }) )
	  {
	    ++$summary{$name}->{$parser};
	  }
	next;
      }
      if ( $file =~ /^script:/ix ) {
	if(!defined($checksum) || $checksum == 0){
	  print "Parsing $file with $parser\n" if ($verbose);
	  my $eval_test = eval "require XrefParser::$parser";
	  if($@ or $eval_test != 1) {
	    croak "Could not require XrefParser::$parser\ndollar=at=$@\neval_test = $eval_test\n";
	  }
	  my $new = "XrefParser::$parser"->new($self->database, $verbose);
	  my $sqlu =
	    "UPDATE source_url SET checksum=1, upload_date=NOW() WHERE source_url_id=$source_url_id";

	  if (
	      $new->run_script( {file       => $file,
				 source_id  => $source_id, 
				 species_id => $species_id, 
				 verbose    => $verbose }) )
	    {
	      ++$summary{$name}->{$parser};
	    }
	  else{
	    # set the checksum to 1 so that we know the script has been ran successfully
	    $dbi->prepare($sqlu)->execute() || croak( $dbi->errstr() );
	  }
	}
	else{
	  print "$file has already been run with $parser and so will not be run again\n" if($verbose);
	}
	next;
      }
	
	
      if ( $unzip && ( $file =~ /\.    # anything
				 (gz|Z) # followed by gz or Z
				 $      # at the end
				 /x ) ) {
	printf( "Uncompressing '%s' using 'gunzip'\n", $file ) if ($verbose);
	system( "gunzip", "-f", $file );
      }
      # remove the gz or Z at the end of the file name
      if ($unzip) { $file =~ s/\.(gz|Z)$//x }
	
      # Compare checksums and parse/upload if necessary need to
      # check file size as some .SPC files can be of zero length
	
      if ( !defined( $cs = md5sum($file) ) ) {
	printf( "Download '%s'\n", $file ) if($verbose);
	++$summary{$name}->{$parser};
      } else {
	$file_cs = md5_hex($file_cs.$cs);
	if ( !defined $checksum
	     || index( $checksum, $file_cs ) == -1 )
	  {
	    if ( -s $file ) {
	      $parse = 1;
	      print "Checksum for '$file' does not match, "
		. "will parse...\n" if ($verbose);

	      # Files from sources "Uniprot/SWISSPROT" and
	      # "Uniprot/SPTREMBL" are all parsed with the
	      # same parser
	      if (    $parser eq "Uniprot/SWISSPROT"
		      || $parser eq "Uniprot/SPTREMBL" )
		{
		  print STDERR "No idea why this is being done here??\n";
		  print STDERR "parser was $parser now being set to UniProtParser\n";
		  $parser = 'UniProtParser';
		}
	    } else {
	      $empty = 1;
	      printf(
		     "The file '%s' has zero length, skipping\n",
		     $file ) if ($verbose);
	    }
	  }
      } ## end else [ if ( !defined( $cs = md5sum...
	
      # Push this file to the list of files to parsed.  The files
      # are *actually* parsed only if $parse == 1.
      push @files_to_parse, $file;
	
    } ## end foreach my $file (@files)

    if ( $parse and @files_to_parse and defined $file_cs ) {
      print "Parsing '"
	. join( "', '", @files_to_parse )
	  . "' with $parser\n" if ($verbose);
	
      eval "require XrefParser::$parser";
      $@ && carp( "[ERROR] Cannot require $parser: $@" );
      my $new = "XrefParser::$parser"->new($self->database, $verbose);###########

      if ( defined $release_url ) {
	# Run with $release_url.
	if (
	    $new->run( { source_id  => $source_id,
			 species_id => $species_id,
			 files      => [@files_to_parse],
			 rel_file   => $release_url,
			 verbose    => $verbose } ) )
	  {
	    ++$summary{$name}->{$parser};
	  }
      } else {
	# Run without $release_url.
	if (
	    $new->run(  { source_id  => $source_id,
			  species_id => $species_id,
			  files      => [@files_to_parse],
			  verbose    => $verbose } ))
	  {
	    ++$summary{$name}->{$parser};
	  }
      }
	
      # update AFTER processing in case of crash.
      $self->update_source( $source_url_id,
		     $file_cs, $files_to_parse[0] );
	
      # Set release if specified
      if ( defined $release ) {
	$self->set_release( $source_id, $release );
      }

    } elsif ( !$dsn && !$empty && @files_to_parse ) {
      print(   "Ignoring '"
	       . join( "', '", @files_to_parse )
	       . "' as checksums match\n" ) if ($verbose);
    }

    if ($cleanup) {
      foreach my $file (@files_to_parse) {
	printf( "Deleting '%s'\n", $file ) if($verbose);
	unlink($file);
      }
    }

    $end_transaction_sth->execute();
    if($stats){
      $self->get_stats({ xrefs   => \%sum_xrefs,
			 primary => \%sum_prim,
			 depend  => \%sum_dep,
			 direct  => \%sum_dir,
			 coord   => \%sum_coord,
			 synonym => \%sum_syn,
			 dbi     => $dbi,
			 print   => 1 })
    }

  } ## end while ( my @row = $sth->fetchrow_array...

  print "\n", '=' x 80, "\n";
  print "Summary of status\n";
  print '=' x 80, "\n";
    
    
  foreach my $source_name ( sort keys %summary ) {
    foreach my $parser_name ( keys %{ $summary{$source_name} } ) {
      printf( "%30s %-20s\t%s\n",
	      $source_name,
	      $parser_name, (
			     defined $summary{$source_name}->{$parser_name}
			     && $summary{$source_name}->{$parser_name}
			     ? 'FAILED'
			     : 'OKAY'
			    ) );
    }
  }
  
  if($stats){
    %sum_xrefs = (); # reset we now want total numbers
    %sum_prim  = ();
    %sum_dep   = ();
    %sum_dir   = ();
    %sum_coord = ();
    %sum_syn   = ();
    $self->get_stats({ xrefs   => \%sum_xrefs,
		       primary => \%sum_prim,
		       depend  => \%sum_dep,
		       direct  => \%sum_dir,
		       coord   => \%sum_coord,
		       synonym => \%sum_syn,
		       dbi     => $dbi,
		       print   => 1 })
  }

  $sth = $dbi->prepare("insert into process_status (status, date) values('parsing_finished',now())");
  $sth->execute;

  return 1;
} ## end sub run



#################################################################################
# Print the statistics for numbers of xrefs, object xref etc that have been added
#################################################################################
sub get_stats {
  my ($self, $ref_arg) = @_;

  my $sum_xrefs = $ref_arg->{xrefs};
  my $sum_prim  = $ref_arg->{primary};
  my $sum_dep   = $ref_arg->{depend};
  my $sum_dir   = $ref_arg->{direct};
  my $sum_coord = $ref_arg->{coord};
  my $sum_syn   = $ref_arg->{synonym};
  my $dbi       = $ref_arg->{dbi};
  my $print     = $ref_arg->{print};

  # produce summary of what has been added
  my %sum_line;

  # first the number of xrefs;
  my $group_sql = "SELECT count(1), s.name from source s, xref x where s.source_id = x.source_id group by s.name";

  my $sum_sth = $dbi->prepare($group_sql);
  $sum_sth->execute();

  my ($sum_count, $sum_name);
  $sum_sth->bind_columns(\$sum_count, \$sum_name);

  while($sum_sth->fetch){
    if(defined($sum_xrefs->{$sum_name})){
      if($sum_count != $sum_xrefs->{$sum_name}){
	my $diff = ($sum_count - $sum_xrefs->{$sum_name});
	$sum_line{$sum_name} = [$diff, 0, 0, 0, 0, 0, 0, 0];
      }
    }
    else{
      $sum_line{$sum_name}  = [$sum_count, 0, 0, 0, 0, 0, ,0 ,0];
    }
    $sum_xrefs->{$sum_name} = $sum_count;
  }
  $sum_sth->finish;


  # second the number of primary xrefs
  $group_sql = "SELECT count(1), s.name from source s, primary_xref px, xref x where s.source_id = x.source_id and px.xref_id = x.xref_id group by s.name";

  $sum_sth = $dbi->prepare($group_sql);
  $sum_sth->execute();

  $sum_sth->bind_columns(\$sum_count, \$sum_name);
 
  while($sum_sth->fetch){
    if ( defined($sum_prim->{$sum_name}) && ($sum_count != $sum_prim->{$sum_name}) ){
      my $diff = ($sum_count - $sum_prim->{$sum_name});
      $sum_line{$sum_name}[1] = $diff; 
    }
    elsif(!defined($sum_prim->{$sum_name})){
      $sum_line{$sum_name}[1] = $sum_count;
    }
    $sum_prim->{$sum_name} = $sum_count;
  }
  $sum_sth->finish;


  # third the number of dependent xrefs
  $group_sql = "SELECT count(1), s.name from source s, dependent_xref dx, xref x where s.source_id = x.source_id and dx.dependent_xref_id = x.xref_id group by s.name";

  $sum_sth = $dbi->prepare($group_sql);
  $sum_sth->execute();

  $sum_sth->bind_columns(\$sum_count, \$sum_name);

  while($sum_sth->fetch){
    if ( defined($sum_dep->{$sum_name}) && ($sum_count != $sum_dep->{$sum_name}) ){
      my $diff = ($sum_count - $sum_dep->{$sum_name});
      $sum_line{$sum_name}[2] = $diff;
    }
    elsif(!defined($sum_dep->{$sum_name})){
      $sum_line{$sum_name}[2] = $sum_count;
    }
    $sum_dep->{$sum_name} = $sum_count;
  }
  $sum_sth->finish;



  # fourth,fifth and sixth the number of direct xrefs

  my $type_count =0;
  foreach my $type (qw (gene transcript translation)){

    $group_sql = "SELECT count(1), s.name from source s, ".$type."_direct_xref dx, xref x where s.source_id = x.source_id and dx.general_xref_id = x.xref_id group by s.name";

    $sum_sth = $dbi->prepare($group_sql);
    $sum_sth->execute();

    $sum_sth->bind_columns(\$sum_count, \$sum_name);

    while($sum_sth->fetch){
      $sum_name .= "_$type";
      if ( defined($sum_dir->{$sum_name}) && ($sum_count != $sum_dir->{$sum_name}) ){
	my $diff = ($sum_count - $sum_dir->{$sum_name});
	$sum_line{$sum_name}[3+$type_count] = $diff;
      }
      elsif(!defined($sum_dir->{$sum_name})){
	$sum_line{$sum_name}[3+$type_count] = $sum_count;
      }
      $sum_dir->{$sum_name} = $sum_count;
    }
    $sum_sth->finish;
    $type_count++;
  }

  # seventh the number of coordinate xrefs
  $group_sql = "SELECT count(1), s.name from source s, coordinate_xref cx  where s.source_id = cx.source_id group by s.name";

  $sum_sth = $dbi->prepare($group_sql);
  $sum_sth->execute();

  $sum_sth->bind_columns(\$sum_count, \$sum_name);

  while($sum_sth->fetch){
    if ( defined($sum_coord->{$sum_name}) && ($sum_count != $sum_coord->{$sum_name}) ){
      my $diff = ($sum_count - $sum_coord->{$sum_name});
      $sum_line{$sum_name}[6] = $diff;
    }
    elsif(!defined($sum_coord->{$sum_name})){
      $sum_line{$sum_name}[6] = $sum_count;
    }
    $sum_coord->{$sum_name} = $sum_count;
  }

  $sum_sth->finish;


  # eigth the number of synonyms
  $group_sql = "select count(1), s.name from source s, xref x, synonym o where s.source_id = x.source_id and x.xref_id = o.xref_id group by s.name";

  $sum_sth = $dbi->prepare($group_sql);
  $sum_sth->execute();

  $sum_sth->bind_columns(\$sum_count, \$sum_name);

  while($sum_sth->fetch){
    if (defined($sum_syn->{$sum_name}) && ($sum_count != $sum_syn->{$sum_name}) ){
      my $diff = ($sum_count - $sum_syn->{$sum_name});
      $sum_line{$sum_name}[7] = $diff;
    }
    elsif(!defined($sum_syn->{$sum_name})) {
      $sum_line{$sum_name}[7] = $sum_count;
    }
    $sum_syn->{$sum_name} = $sum_count;
  }
  $sum_sth->finish;

  if($print){
    ###################
    # Print the header
    ###################
    my $max_name_length = 6; # (source)
    my $width = 8;
    if(scalar(keys %sum_line)){
      foreach my $sum_name (keys %sum_line){
	if(length($sum_name) > $max_name_length){
	  $max_name_length = length($sum_name);
	}
      }

      print "\nsource". " " x ($max_name_length - 3); #( 3 = length(source) - 3 spaces)
      foreach my $val (qw(xrefs prim dep gdir tdir tdir coord synonyms)){
	print $val." " x ($width - length($val) );
      }
      print "\n";
    }

    ###################
    # Print the numbers
    ###################
    $max_name_length += 3; # lets have 3 spaces after
    foreach my $sum_name (keys %sum_line){
      $sum_name ||= 0;
      print $sum_name. " " x ( $max_name_length - length($sum_name));
      foreach my $val (@{$sum_line{$sum_name}}){
	$val ||= 0;
	print $val." " x ($width - length($val));
      }
      print "\n";
    }
    print "\n";
  }
  return;
}


###################################
# Getter/Setter for database object
###################################
sub database {
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_database} = $arg );
  return $self->{_database};
}


##############################
# Getter/Setter for dbi object
##############################
sub dbi {
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_dbi} = $arg );
  return $self->{_dbi};
}


###########################################################
# Check if all the species name in a ref to list is valid
#   If they all are; returns a list of species
#   It not prints the values that are aloud and dies.
###########################################################
sub validate_species {
  my ($self, $species, $verbose) = @_;
  my @species_ids;

  my $dbi = $self->dbi();
  my $sth = $dbi->prepare("SELECT species_id, name FROM species WHERE LOWER(name)=? OR LOWER(aliases) REGEXP ?");
  my ($species_id, $species_name);

  foreach my $sp (@$species) {

    my $bind_arg = "^".lc($sp).",|^".lc($sp)."\$|,[ ]{0,1}".lc($sp)."[ ]{0,1},|,[ ]{0,1}".lc($sp)."\$";
    $sth->execute(lc($sp), $bind_arg ); 
    $sth->bind_columns(\$species_id, \$species_name);
    if (my @row = $sth->fetchrow_array()) {
      print "Species $sp is valid (name = " . $species_name . ", ID = " . $species_id . ")\n" if($verbose);
      push @species_ids, $species_id;
    } else {
      print STDERR "Species $sp is not valid; valid species are:\n";
      $self->show_valid_species();
      exit(1);
    }
  }
  return @species_ids;
}

############################################################
# Return 1 if all dependent sources have already been loaded
# else return 0;
############################################################
sub all_dependencies_loaded{
  my ($self, $source_id, $species_id, $s_name, $dep_sth) = @_;
  my $okay = 1;

  $dep_sth->execute($source_id, $species_id, $source_id);
  my ($id, $name);
  $dep_sth->bind_columns(\$id, \$name);
  while($dep_sth->fetch() ){
    print STDERR "dependent source $name ($id) not loaded so cannot process source $s_name\n";
    print "dependent source $name ($id) not loaded so cannot process source $s_name\n";
    $okay = 0;
  }
  return $okay;
}


########################################################################
# Remove potentially problematic characters from string used as file or
# directory names.
########################################################################
sub sanitise {
    my $str = shift;
    $str =~ tr[/:][]d;
    return $str;
}


#######################################################
# Print to stanadrd error the list of species available
#######################################################
sub show_valid_species {
  my ($self) =shift;

  my $dbi = $self->dbi();
  my $sth = $dbi->prepare("SELECT name, aliases FROM species");

  $sth->execute();
  while (my @row = $sth->fetchrow_array()) {
    print STDERR $row[0] . " (aliases: " . $row[1] . ")\n";
  }
  return;
}


#########################################################
# Check if all the source names in a ref to list is valid
#   If they all are; returns 1
#   It not prints the values that are aloud and returns 0.
#########################################################
sub validate_sources {
  my ($self, $speciesref, $sources, $verbose) = @_;

  my $dbi = $self->dbi();
  my $sth = $dbi->prepare("SELECT * FROM source WHERE LOWER(name)=?");

  foreach my $source (@$sources) {

    my $rv = $sth->execute(lc($source));
    if ( $rv > 0 ) {
      print "Source $source is valid\n" if($verbose);
    } else {
      print "\nSource $source is not valid; valid sources are:\n";
      foreach my $sp (@{$speciesref}){
	show_valid_sources($sp);
      }
      return 0;
    }

  }

  return 1;

}


#######################################################
# Print to stanadrd error the list of sources available
#######################################################
sub show_valid_sources {
  my ($self, $species) = @_;

  my $dbi = $self->dbi();
  my $sth = $dbi->prepare("SELECT distinct(name) FROM source s, source_url su WHERE s.download='Y' and s.source_id = su.source_id and su.species_id = $species");

  $sth->execute();
  while (my @row = $sth->fetchrow_array()) {
    print $row[0] . "\n";
  }
  return;
}


####################################################################
# Compute a checksum of a file.  This checksum is not a straight MD5
# hex digest, but instead the file size combined with the first six
# characters of the MD5 hex digest.  This is to save space.
####################################################################
sub md5sum
{
    my $file = shift;

    open my $FH, "<", $file or return;
    binmode($FH);

    my $checksum = sprintf( "%s/%d",
        substr( Digest::MD5->new()->addfile(*$FH)->hexdigest(), 0, 6 ),
        [ stat $FH ]->[7] );

    close($FH);

    return $checksum;
}


####################################
# Set the checksum for a source file
####################################
sub update_source {
  my ($self, $source_url_id, $checksum, $file_name ) = @_;

  my $dbi = $self->dbi();
  my $file = IO::File->new($file_name)
    or croak("Failed to open file '$file_name'");

  my $file_date =
    POSIX::strftime( '%Y%m%d%H%M%S',
		     localtime( [ $file->stat() ]->[9] ) );

  $file->close();

  my $sql =
    "UPDATE source_url SET checksum='$checksum', "
      . "file_modified_date='$file_date', "
	. "upload_date=NOW() "
	  . "WHERE source_url_id=$source_url_id";

  # The release is set by the individual parser by calling the
  # inherited set_release() method.

  $dbi->prepare($sql)->execute() || croak( $dbi->errstr() );
  return;
}


1;
