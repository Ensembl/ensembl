package XrefParser::BaseParser;

use strict;

use Carp;
use DBI;
use Digest::MD5 qw(md5_hex);
use Getopt::Long;
use POSIX qw(strftime);

use File::Basename;
use File::Spec::Functions;
use IO::File;
use Net::FTP;
use URI;
use URI::file;
use Text::Glob qw( match_glob );
use LWP::UserAgent;

use Bio::EnsEMBL::Utils::Exception;

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

my ( $host,             $port,    $dbname,        $user,
     $pass,             $create,  $release,       $cleanup,
     $deletedownloaded, $drop_db, $checkdownload, $dl_path,
     $unzip, $stats, $verbose);


# --------------------------------------------------------------------------------
# Get info about files to be parsed from the database

sub run {
    my $self = shift;

    (  $host,             $port,          $dbname,
       $user,             $pass,          my $speciesr,
       my $sourcesr,      $checkdownload, $create,
       $release,          $cleanup,       $drop_db,
       $deletedownloaded, $dl_path,       my $notsourcesr,
       $unzip, $stats, $verbose
    ) = @_;

    $base_dir = $dl_path if $dl_path;

    my @species    = @$speciesr;
    my @sources    = @$sourcesr;
    my @notsources = @$notsourcesr;

    my $sql_dir = dirname($0);

    if ($create) {
      create( $host, $port, $user, $pass, $dbname, $sql_dir, $drop_db );
    }

    my $dbi = dbi();
    my $sth_c = $dbi->prepare("insert into process_status (status, date) values('xref_created',now())");
    $sth_c->execute;

    # validate species names
    my @species_ids = validate_species(@species);

    # validate source names
    exit(1) if ( !validate_sources(\@species_ids,@sources) );
    exit(1) if ( !validate_sources(\@species_ids,@notsources) );

    # build SQL
    my $species_sql = "";
    if (@species_ids) {
        $species_sql .= " AND su.species_id IN (";
        for ( my $i = 0 ; $i < @species_ids ; $i++ ) {
            $species_sql .= "," if ( $i ne 0 );
            $species_sql .= $species_ids[$i];
        }
        $species_sql .= ") ";
    }

    my $source_sql = "";
    if (@sources) {
        $source_sql .= " AND LOWER(s.name) IN (";
        for ( my $i = 0 ; $i < @sources ; $i++ ) {
            $source_sql .= "," if ( $i ne 0 );
            $source_sql .= "\'" . lc( $sources[$i] ) . "\'";
        }
        $source_sql .= ") ";
    }

    if (@notsources) {
        $source_sql .= " AND LOWER(s.name) NOT IN (";
        for ( my $i = 0 ; $i < @notsources ; $i++ ) {
            $source_sql .= "," if ( $i ne 0 );
            $source_sql .= "\'" . lc( $notsources[$i] ) . "\'";
        }
        $source_sql .= ") ";
    }

    my $sql =
      "SELECT DISTINCT(s.source_id), su.source_url_id, s.name, su.url, "
      . "su.release_url, su.checksum, su.parser, su.species_id "
      . "FROM source s, source_url su, species sp "
      . "WHERE s.download='Y' AND su.source_id=s.source_id "
      . "AND su.species_id=sp.species_id "
      . $source_sql
      . $species_sql
      . "ORDER BY s.ordered";
    #print $sql . "\n";

    my $sth = $dbi->prepare("insert into process_status (status, date) values('parsing_started',now())");
    $sth->execute;

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
    my %sum_list;
    my %sum_syn;


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

        my @files = split( /\s+/, $url );
        my @files_to_parse = ();

        $dir = catdir( $base_dir, sanitise($type) );

        # For summary purposes: If 0 is returned (in
        # $summary{$name}->{$parser}) then it is successful.  If 1 is
        # returned then it failed.  If undef/nothing is returned the we
        # do not know.
        $summary{$name}->{$parser} = 0;

        @files = $self->fetch_files( $dir, @files );
        if ( !@files ) {
            # Fetching failed.
            ++$summary{$name}->{$parser};
            next;
        }
        if ( defined($release_url) ) {
            $release_url =
              $self->fetch_files( $dir, $release_url )->[-1];
        }
	$start_transaction_sth->execute();

        foreach my $file (@files) {

            # Database parsing
            if ( $file =~ /^mysql:/i ) {
                $dsn = $file;
                print "Parsing $dsn with $parser\n" if ($verbose);
                eval "require XrefParser::$parser";
                my $new = "XrefParser::$parser"->new();
                if (
                     $new->run( $dsn,  $source_id, $species_id,
                                $name, undef, $verbose ) )
                {
                    ++$summary{$name}->{$parser};
                }
                next;
            }
	    if ( $file =~ /^script:/i ) {
	      if(!defined($checksum) || $checksum == 0){
   	        print "Parsing $file with $parser\n" if ($verbose);
		eval "require XrefParser::$parser";
		my $new = "XrefParser::$parser"->new();
		if (
		    $new->run_script( $file,  $source_id, $species_id, $verbose ) )
		  {
		    ++$summary{$name}->{$parser};
		  }
                # set the checksum to 1 so that we know the script has been ran
		my $sqlu =
		  "UPDATE source_url SET checksum=1, upload_date=NOW() WHERE source_url_id=$source_url_id";
		
		$dbi->prepare($sqlu)->execute() || croak( $dbi->errstr() );		
	      }
              else{
                print "$file has already been ran with $parser and so will not be ran again\n" if($verbose);
	      }
	      next;
	    }
	    
	    
	    if ( $unzip && ( $file =~ /\.(gz|Z)$/ ) ) {
	      printf( "Uncompressing '%s' using 'gunzip'\n", $file ) if ($verbose);
                system( "gunzip", "-f", $file );
            }
            if ($unzip) { $file =~ s/\.(gz|Z)$// }

            # Compare checksums and parse/upload if necessary need to
            # check file size as some .SPC files can be of zero length

            if ( !defined( $cs = md5sum($file) ) ) {
                printf( "Download '%s'\n", $file ) if($verbose);
                ++$summary{$name}->{$parser};
            } else {
                $file_cs .= ':' . $cs;
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
            $@ && warn( "[ERROR] Cannot require $parser: $@" );
            my $new = "XrefParser::$parser"->new();

            if ( defined $release_url ) {
                # Run with $release_url.
                if (
                     $new->run( $source_id,      $species_id,
                                \@files_to_parse, $release_url, $verbose ) )
                {
                    ++$summary{$name}->{$parser};
                }
            } else {
                # Run without $release_url.
                if (
                     $new->run( $source_id, $species_id,
                                \@files_to_parse, undef, $verbose ) )
                {
                    ++$summary{$name}->{$parser};
                }
            }

            # update AFTER processing in case of crash.
            update_source( $dbi,     $source_url_id,
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
	if($stats){
	# produce summary of what has been added
	  my %sum_line;
	  
	  # first the number of xrefs;
	  my $group_sql = "SELECT count(*), s.name from source s, xref x where s.source_id = x.source_id group by s.name";
	  
	  my $sum_sth = $dbi->prepare($group_sql);
	  $sum_sth->execute();
	  
	  my ($sum_count, $sum_name);
	  $sum_sth->bind_columns(\$sum_count, \$sum_name);
	  
	  while($sum_sth->fetch){
	    if(defined($sum_xrefs{$sum_name})){
	      if($sum_count != $sum_xrefs{$sum_name}){
		my $diff = ($sum_count - $sum_xrefs{$sum_name});
		$sum_line{$sum_name} = [$diff, 0, 0, 0, 0, 0];	    
	      }
#	      else{
#		$sum_line{$sum_name} = [0, 0, 0, 0, 0, 0];
#	      }
	    }
	    else{
	      $sum_line{$sum_name}  = [$sum_count, 0, 0, 0, 0, 0];
	    }
	    $sum_xrefs{$sum_name} = $sum_count;
	  }
	  $sum_sth->finish;
	  

	  # second the number of primary xrefs
	  $group_sql = "SELECT count(*), s.name from source s, primary_xref px, xref x where s.source_id = x.source_id and px.xref_id = x.xref_id group by s.name";
	  
	  my $sum_sth = $dbi->prepare($group_sql);
	  $sum_sth->execute();
	  
	  $sum_sth->bind_columns(\$sum_count, \$sum_name);
	  
	  while($sum_sth->fetch){
	    if($sum_count != $sum_prim{$sum_name}){
	      my $diff = ($sum_count - $sum_prim{$sum_name});
	      $sum_line{$sum_name}[1] = $diff;	    
	    }
	    $sum_prim{$sum_name} = $sum_count;
	  }
	  $sum_sth->finish;
	  
	  
	  
	  # third the number of dependent xrefs
	  $group_sql = "SELECT count(*), s.name from source s, dependent_xref dx, xref x where s.source_id = x.source_id and dx.dependent_xref_id = x.xref_id group by s.name";
	  
	  my $sum_sth = $dbi->prepare($group_sql);
	  $sum_sth->execute();
	  
	  $sum_sth->bind_columns(\$sum_count, \$sum_name);
	  
	  while($sum_sth->fetch){
	    if($sum_count != $sum_dep{$sum_name}){
	      my $diff = ($sum_count - $sum_dep{$sum_name});
	      $sum_line{$sum_name}[2] = $diff;	    
	    }
	    $sum_dep{$sum_name} = $sum_count;
	  }
	  $sum_sth->finish;

	  
	  
	  
	  # fourth,fifth and sixth the number of direct xrefs

	  my $type_count =0;
	  foreach my $type (qw (gene transcript translation)){

	    $group_sql = "SELECT count(*), s.name from source s, ".$type."_direct_xref dx, xref x where s.source_id = x.source_id and dx.general_xref_id = x.xref_id group by s.name";
	  
	    my $sum_sth = $dbi->prepare($group_sql);
	    $sum_sth->execute();
	    
	    $sum_sth->bind_columns(\$sum_count, \$sum_name);
	    
	    while($sum_sth->fetch){
	      $sum_name .= "_$type";
	      if($sum_count != $sum_dir{$sum_name}){
		my $diff = ($sum_count - $sum_dir{$sum_name});
		$sum_line{$sum_name}[3+$type_count] = $diff;	    
	      }
	      $sum_dir{$sum_name} = $sum_count;
	    }
	    $sum_sth->finish;
	    $type_count++;
	  }

	  # seventh the number of coordinate xrefs
	  $group_sql = "SELECT count(*), s.name from source s, coordinate_xref cx  where s.source_id = cx.source_id group by s.name";
	  
	  my $sum_sth = $dbi->prepare($group_sql);
	  $sum_sth->execute();
	  
	  $sum_sth->bind_columns(\$sum_count, \$sum_name);
	  
	  while($sum_sth->fetch){
	    if($sum_count != $sum_coord{$sum_name}){
	      my $diff = ($sum_count - $sum_coord{$sum_name});
	      $sum_line{$sum_name}[6] = $diff;	    
	    }
	    $sum_coord{$sum_name} = $sum_count;
	  }
	  $sum_sth->finish;
	  

	  # eigth the number of synonyms
	  $group_sql = "select count(*), s.name from source s, xref x, synonym o where s.source_id = x.source_id and x.xref_id = o.xref_id group by s.name";
	  
	  my $sum_sth = $dbi->prepare($group_sql);
	  $sum_sth->execute();
	  
	  $sum_sth->bind_columns(\$sum_count, \$sum_name);
	  
	  while($sum_sth->fetch){
	    if($sum_count != $sum_syn{$sum_name}){
	      my $diff = ($sum_count - $sum_syn{$sum_name});
	      $sum_line{$sum_name}[7] = $diff;	    
	    }
	    $sum_syn{$sum_name} = $sum_count;
	  }
	  $sum_sth->finish;


	  print "source                      xrefs\tprim\tdep\tgdir\ttdir\ttdir\tcoord\tsynonyms\n";
	  foreach my $sum_name (keys %sum_line){
	    printf ("%-28s",$sum_name);
	    print join("\t",@{$sum_line{$sum_name}})."\n";
	  }
	  
	} # if ($stats)	
	
	$end_transaction_sth->execute();

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
      my %sum_line;
      
      # first the number of xrefs;
      my $group_sql = "SELECT count(*), s.name from source s, xref x where s.source_id = x.source_id group by s.name";
      
      my $sum_sth = $dbi->prepare($group_sql);
      $sum_sth->execute();
      
      my ($sum_count, $sum_name);
      $sum_sth->bind_columns(\$sum_count, \$sum_name);
      
      while($sum_sth->fetch){
	$sum_line{$sum_name} = [$sum_count, 0, 0, 0, 0, 0];
      }
      $sum_sth->finish;
      
      
      # second the number of primary xrefs
      $group_sql = "SELECT count(*), s.name from source s, primary_xref px, xref x where s.source_id = x.source_id and px.xref_id = x.xref_id group by s.name";
      
      my $sum_sth = $dbi->prepare($group_sql);
      $sum_sth->execute();
      
      $sum_sth->bind_columns(\$sum_count, \$sum_name);
      
      while($sum_sth->fetch){
	$sum_line{$sum_name}[1] = $sum_count;	    
      }
      $sum_sth->finish;
      
      
      
      # third the number of dependent xrefs
      $group_sql = "SELECT count(*), s.name from source s, dependent_xref dx, xref x where s.source_id = x.source_id and dx.dependent_xref_id = x.xref_id group by s.name";
      
      my $sum_sth = $dbi->prepare($group_sql);
      $sum_sth->execute();
      
      $sum_sth->bind_columns(\$sum_count, \$sum_name);
      
      while($sum_sth->fetch){
	$sum_line{$sum_name}[2] = $sum_count;	    
      }
      
      
      
      # fourth,fifth and sixth the number of direct xrefs
      
      my $type_count =0;
      foreach my $type (qw (gene transcript translation)){
	
	$group_sql = "SELECT count(*), s.name from source s, ".$type."_direct_xref dx, xref x where s.source_id = x.source_id and dx.general_xref_id = x.xref_id group by s.name";
	
	my $sum_sth = $dbi->prepare($group_sql);
	$sum_sth->execute();
	
	$sum_sth->bind_columns(\$sum_count, \$sum_name);
	
	while($sum_sth->fetch){
	  $sum_line{$sum_name}[3+$type_count] = $sum_count;	    
	}
	$sum_sth->finish;
	$type_count++;
      }
      
      # seventh the number of coordinate xrefs
      $group_sql = "SELECT count(*), s.name from source s, coordinate_xref cx  where s.source_id = cx.source_id group by s.name";
      
      my $sum_sth = $dbi->prepare($group_sql);
      $sum_sth->execute();
      
      $sum_sth->bind_columns(\$sum_count, \$sum_name);
      
      while($sum_sth->fetch){
	$sum_line{$sum_name}[6] = $sum_count;	    
      }
      $sum_sth->finish;
      
      
      # eigth the number of synonyms
      $group_sql = "select count(*), s.name from source s, xref x, synonym o where s.source_id = x.source_id and x.xref_id = o.xref_id group by s.name";
      
      my $sum_sth = $dbi->prepare($group_sql);
      $sum_sth->execute();
      
      $sum_sth->bind_columns(\$sum_count, \$sum_name);
      
      while($sum_sth->fetch){
	$sum_line{$sum_name}[7] = $sum_count;	    
      }
    
      print "---------------------------------------------------------------------------------------\n";
      print "TOTAL source                xrefs\tprim\tdep\tgdir\ttdir\ttdir\tcoord\tsynonyms\n";
      foreach my $sum_name (keys %sum_line){
	printf ("%-28s",$sum_name);
	print join("\t",@{$sum_line{$sum_name}})."\n";
      }
      
    }

    $sth = $dbi->prepare("insert into process_status (status, date) values('parsing_finished',now())");
    $sth->execute;

    # remove last working directory
    # TODO reinstate after debugging
    #rmtree $dir;

} ## end sub run

# ------------------------------------------------------------------------------

# Given one or several FTP or HTTP URIs, download them.  If an URI is
# for a file or MySQL connection, then these will be ignored.  For
# FTP, standard shell file name globbing is allowed (but not regular
# expressions).  HTTP does not allow file name globbing.  The routine
# returns a list of successfully downloaded local files or an empty list
# if there was an error.

sub fetch_files {
    my $self = shift;

    my ( $dest_dir, @user_uris ) = @_;

    my @processed_files;

    foreach my $user_uri (@user_uris) {
        # Change old-style 'LOCAL:' URIs into 'file:'.
        $user_uri =~ s#^LOCAL:#file:#i;
        my $uri = URI->new($user_uri);

#	print "\n*******$user_uri\n*********\n";
	if ( $uri->scheme() eq 'script'){
	  push( @processed_files, $user_uri );	  
	}elsif ( $uri->scheme() eq 'file' ) {
            # Deal with local files.

            my @local_files;

            $user_uri =~ s/file://;
            if ( -f $user_uri ) {
                push( @processed_files, $user_uri );
            } else {
                printf( "==> Can not find file '%s'\n", $user_uri );
                return ();
            }
        } elsif ( $uri->scheme() eq 'ftp' ) {
            # Deal with FTP files.

            my $file_path =
              catfile( $dest_dir, basename( $uri->path() ) );

            if ( $deletedownloaded && -f $file_path ) {
                printf( "Deleting '%s'\n", $file_path ) if ($verbose);
                unlink($file_path);
            }

            if ( $checkdownload && -f $file_path ) {
                # The file is already there, no need to connect to a FTP
                # server.  This also means no file name globbing was
                # used (for globbing FTP URIs, we always need to connect
                # to a FTP site to see what files are there).

                printf( "File '%s' already exists\n", $file_path ) if ($verbose);
                push( @processed_files, $file_path );
                next;
            }

            printf( "Connecting to FTP host '%s' for file '%s' \n", $uri->host(), $file_path ) if ($verbose);

            my $ftp = Net::FTP->new( $uri->host(), 'Debug' => 0 );
            if ( !defined($ftp) ) {
                printf( "==> Can not open FTP connection: %s\n", $@ );
                return ();
            }

            if ( !$ftp->login( 'anonymous', '-anonymous@' ) ) {
                printf( "==> Can not log in on FTP host: %s\n",
                        $ftp->message() );
                return ();
            }

            if ( !$ftp->cwd( dirname( $uri->path() ) ) ) {
                printf( "== Can not change directory to '%s': %s\n",
                        dirname( $uri->path() ), $ftp->message() );
                return ();
            }

            $ftp->binary();

            foreach my $remote_file ( ( @{ $ftp->ls() } ) ) {
                if (
                     !match_glob( basename( $uri->path() ), $remote_file
                     ) )
                {
                    next;
                }

                $file_path =
                  catfile( $dest_dir, basename($remote_file) );

                if ( $deletedownloaded && -f $file_path ) {
                    printf( "Deleting '%s'\n", $file_path ) if($verbose);
                    unlink($file_path);
                }

                if ( $checkdownload && -f $file_path ) {
                    printf( "File '%s' already exists\n", $file_path ) if ($verbose);
                } else {

                    if ( !-d dirname($file_path) ) {
                        printf( "Creating directory '%s'\n",
                                dirname($file_path) ) if($verbose);
                        if ( !mkdir( dirname($file_path) ) ) {
                            printf(
                                "==> Can not create directory '%s': %s",
                                dirname($file_path), $! );
                            return ();
                        }
                    }

                    printf( "Fetching '%s' (size = %s)\n",
                            $remote_file,
                            $ftp->size($remote_file) || '(unknown)' ) if ($verbose);
                    printf( "Local file is '%s'\n", $file_path ) if($verbose);

                    if ( !$ftp->get( $remote_file, $file_path ) ) {
                        printf( "==> Could not get '%s': %s\n",
                                basename( $uri->path() ),
                                $ftp->message() );
                        return ();
                    }
                } ## end else [ if ( $checkdownload &&...

		if ( $file_path =~ /\.(gz|Z)$/ ) {
		  # Read from zcat pipe
		  #
		  my $cmd = "gzip -t $file_path" ;
		  if(system($cmd) != 0 ){
		    print "system $cmd  failed: $? - Checking of gzip file failed - FILE CORRUPTED ?\n\n";
		    
		    if (-f $file_path ) {
		      print ( "Deleting '%s'\n", $file_path ) if($verbose);
		      unlink($file_path);
		    }
		    return ();
		  }
		  else{
		    print "$file_path passed (gzip -t) corruption test.\n" if($verbose);
		  }
		}
                push( @processed_files, $file_path );
		
	      } ## end foreach my $remote_file ( (...

        } elsif ( $uri->scheme() eq 'http' ) {
            # Deal with HTTP files.

            my $file_path =
              catfile( $dest_dir, basename( $uri->path() ) );

            if ( $deletedownloaded && -f $file_path ) {
                printf( "Deleting '%s'\n", $file_path ) if($verbose);
                unlink($file_path);
            }

            if ( $checkdownload && -f $file_path ) {
                # The file is already there, no need to connect to a
                # HTTP server.

                printf( "File '%s' already exists\n", $file_path ) if ($verbose);
                push( @processed_files, $file_path );
                next;
            }

            if ( !-d dirname($file_path) ) {
                printf( "Creating directory '%s'\n",
                        dirname($file_path) ) if($verbose);
                if ( !mkdir( dirname($file_path) ) ) {
                    printf( "==> Can not create directory '%s': %s",
                            dirname($file_path), $! );
                    return ();
                }
            }

            printf( "Connecting to HTTP host '%s'\n", $uri->host() ) if($verbose);
            printf( "Fetching '%s'\n",                $uri->path() ) if($verbose);

            if ( $checkdownload && -f $file_path ) {
                printf( "File '%s' already exists\n", $file_path ) if($verbose);
            } else {

                printf( "Local file is '%s'\n", $file_path ) if($verbose);

                my $ua = LWP::UserAgent->new();
                $ua->env_proxy();

                my $response = $ua->get( $uri->as_string(),
                                        ':content_file' => $file_path );

                if ( !$response->is_success() ) {
                    printf( "==> Could not get '%s': %s\n",
                            basename( $uri->path() ),
                            $response->content() );
                    return ();
                }
            }

            push( @processed_files, $file_path );

        } elsif ( $uri->scheme() eq 'mysql' ) {
            # Just leave MySQL data untouched for now.
            push( @processed_files, $user_uri );
        } else {
            printf( "==> Unknown URI scheme '%s' in URI '%s'\n",
                    $uri->scheme(), $uri->as_string() );
            return ();
        }
    } ## end foreach my $user_uri (@user_uris)

    return ( wantarray() ? @processed_files : \@processed_files );
} ## end sub fetch_files

# Given a file name, returns a IO::Handle object.  If the file is
# gzipped, the handle will be to an unseekable stream coming out of a
# zcat pipe.  If the given file name doesn't correspond to an existing
# file, the routine will try to add '.gz' to the file name or to remove
# any .'Z' or '.gz' and try again.  Returns undef on failure and will
# write a warning to stderr.

sub get_filehandle
{
    my ($self, $file_name) = @_;

    my $io;

    my $alt_file_name = $file_name;
    $alt_file_name =~ s/\.(gz|Z)$//;

    if ( $alt_file_name eq $file_name ) {
        $alt_file_name .= '.gz';
    }

    if ( !-f $file_name ) {
        carp(   "File '$file_name' does not exist, "
              . "will try '$alt_file_name'" );
        $file_name = $alt_file_name;
    }

    if ( $file_name =~ /\.(gz|Z)$/ ) {
        # Read from zcat pipe
        $io = IO::File->new("zcat $file_name |")
          or carp("Can not open file '$file_name' with 'zcat'");
    } else {
        # Read file normally
        $io = IO::File->new($file_name)
          or carp("Can not open file '$file_name'");
    }

    if ( !defined $io ) { return undef }

    print "Reading from '$file_name'...\n" if($verbose);

    return $io;
}

# ------------------------------------------------------------------------------

sub new
{
    my ($proto) = @_;

    my $class = ref $proto || $proto;
    return bless {}, $class;
}

# --------------------------------------------------------------------------------
# Get source ID for a particular file; matches url field

sub get_source_id_for_filename {

  my ($self, $file) = @_;
  print "FILE $file\n" if($verbose) ; 
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
  
  my ($self, $source_name,$priority_desc) = @_;
  my $sql = "SELECT source_id FROM source WHERE LOWER(name)='" . lc($source_name) . "'";
  if(defined($priority_desc)){
    $sql .= " AND LOWER(priority_description)='".lc($priority_desc)."'";
  }
  my $sth = dbi()->prepare($sql);
  $sth->execute();
  my @row = $sth->fetchrow_array();
  my $source_id;
  if (@row) {
    $source_id = $row[0]; 
  } else {
    print STDERR "WARNING: There is no entity $source_name in the source-table of the xref database.\n" .
      "WARNING:. The external db name ($source_name) is hardcoded in the parser\n";
    warn("WARNING: Couldn't get source ID for source name $source_name\n");

    $source_id = -1;
  }
  return $source_id;
}



# --------------------------------------------------------------------------------
# Get a set of source IDs matching a source name pattern

sub get_source_ids_for_source_name_pattern {

  my ($self, $source_name) = @_;

  my $sql = "SELECT source_id FROM source WHERE upper(name) LIKE '%".uc($source_name)."%'";

  my $sth = dbi()->prepare($sql);
  my @sources;
  $sth->execute();
  while(my @row = $sth->fetchrow_array()){
    push @sources,$row[0];
  }
  $sth->finish;

  return @sources;

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


  my $sql = "select source_id from source where LOWER(name) =?";
  my $sth = dbi()->prepare($sql);
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

  $sql  = "select d.master_xref_id, x2.accession ";
  $sql .= "  from dependent_xref d, xref x1, xref x2 ";
  $sql .= "    where x1.xref_id = d.master_xref_id and";
  $sql .= "          x1.source_id=? and ";
  $sql .= "          x2.xref_id = d.dependent_xref_id and";
  $sql .= "          x2.source_id=? ";
  
  $sth = dbi()->prepare($sql);
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

  $sql  = "select d.general_xref_id, d.ensembl_stable_id, 'Gene', d.linkage_xref, x1.accession ";
  $sql .= "  from gene_direct_xref d, xref x1 ";
  $sql .= "    where x1.xref_id = d.general_xref_id and";
  $sql .= "          x1.source_id=?";
   
  my $sth1 = dbi()->prepare($sql);


  $sql  = "select d.general_xref_id, d.ensembl_stable_id, 'Transcript', d.linkage_xref, x1.accession ";
  $sql .= "  from transcript_direct_xref d, xref x1 ";
  $sql .= "    where x1.xref_id = d.general_xref_id and";
  $sql .= "          x1.source_id=?";
   
  my $sth2 = dbi()->prepare($sql);


  $sql  = "select d.general_xref_id, d.ensembl_stable_id, 'Translation', d.linkage_xref, x1.accession ";
  $sql .= "  from translation_direct_xref d, xref x1 ";
  $sql .= "    where x1.xref_id = d.general_xref_id and";
  $sql .= "          x1.source_id=?";
   
  my $sth3 = dbi()->prepare($sql);

  foreach my $d (@direct_sources){
    $sth1->execute($d);
    while(my @row = $sth1->fetchrow_array()){
      $direct_2_xref{$row[4]} = $row[0]."::".$row[1]."::".$row[2]."::".$row[3];
    }    
    $sth2->execute($d);
    while(my @row = $sth2->fetchrow_array()){
      $direct_2_xref{$row[4]} = $row[0]."::".$row[1]."::".$row[2]."::".$row[3];
    }    
    $sth3->execute($d);
    while(my @row = $sth3->fetchrow_array()){
      $direct_2_xref{$row[4]} = $row[0]."::".$row[1]."::".$row[2]."::".$row[3];
    }
  }

  return \%direct_2_xref;
}

sub label_to_acc{

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
    $sql = "select label, xref_id from xref where species_id = $species_id and source_id = $source";
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



# --------------------------------------------------------------------------------



sub get_existing_mappings {

  my ($self, $from_source_name, $to_source_name, $species_id) =@_;

  my %mappings;

  my $from_source = $self->get_source_id_for_source_name($from_source_name);
  my $to_source = $self->get_source_id_for_source_name($to_source_name);

  my $sql = "SELECT dx.dependent_xref_id, x1.accession as dependent, dx.master_xref_id, x2.accession as master FROM dependent_xref dx, xref x1, xref x2 WHERE x1.xref_id=dx.dependent_xref_id AND x2.xref_id=dx.master_xref_id AND x2.source_id=? AND x1.source_id=? AND x1.species_id=? AND x2.species_id=?";

  my $sth = dbi()->prepare($sql);
  $sth->execute($to_source, $from_source, $species_id, $species_id);
  while(my @row = $sth->fetchrow_array()){
    $mappings{$row[1]} = $row[2];
    #print "mgi_to_uniprot{" . $row[1] . "} = " . $row[2] . "\n";
  }

  print "Got " . scalar(keys(%mappings)) . " $from_source_name -> $to_source_name mappings\n" if($verbose);

  return \%mappings;

}

# --------------------------------------------------------------------------------
# Upload xrefs to the database

sub upload_xref_object_graphs {
  my ($self, $rxrefs) = @_;

  my $dbi = dbi();
  print "count = ".$#$rxrefs."\n" if($verbose);

  if ($#$rxrefs > -1) {

    # remove all existing xrefs with same source ID(s)
#    $self->delete_by_source($rxrefs);

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
			 $xref->{VERSION} || 0,
			 $xref->{LABEL}|| $xref->{ACCESSION},
			 $xref->{DESCRIPTION},
			 $xref->{SOURCE_ID},
			 $xref->{SPECIES_ID},
			 $xref->{INFO_TYPE} || "MISC")){
	 if(!defined($xref->{SOURCE_ID})){
	   print "your xref: $xref->{ACCESSION} does not have a source-id\n";
	   return undef;
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

	if($dbi->err){
	  print STDERR "dbi\t$dbi->err \n$dep{ACCESSION} \n $dep{SOURCE_ID} \n";
	}
	if(!defined($dep_xref_id) || $dep_xref_id ==0 ){
	  print STDERR "acc = $dep{ACCESSION} \nlink = $dep{LINKAGE_SOURCE_ID} \n".$dbi->err."\n";
	  print STDERR "source = $dep{SOURCE_ID}\n";
	}
        $dep_sth->execute( $xref_id, $dep_xref_id,
            $dep{LINKAGE_ANNOTATION},
            $dep{LINKAGE_SOURCE_ID} )
          or croak( $dbi->errstr() );
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
#    print "having now direct-XREF : ".$dr->{ENSEMBL_STABLE_ID}."\t".$dr->{SPECIES_ID}." \n" ;
    my $general_xref_id = get_xref($dr->{ACCESSION},$dr->{SOURCE_ID},$dr->{SPECIES_ID});
    if ($general_xref_id){
      # print "direct_xref:\n$general_xref_id\n$dr->{ENSEMBL_STABLE_ID}\n$dr->{ENSEMBL_TYPE}\t$dr->{LINKAGE_XREF}\n\n";
      $self->add_direct_xref($general_xref_id, $dr->{ENSEMBL_STABLE_ID},$dr->{ENSEMBL_TYPE},$dr->{LINKAGE_XREF});
    }
  }
}

sub add_meta_pair {

  my ($self, $key, $value) = @_;

  my $dbi = dbi();
  my $sth = $dbi->prepare('insert into meta (meta_key, meta_value, date) values("'.$key.'", "'.$value.'", now())');
  $sth->execute;
  $sth->finish;

}



# --------------------------------------------------------------------------------
# Get & cache a hash of all the source names for dependent xrefs (those that are
# in the source table but don't have an associated URL etc)

sub get_dependent_xref_sources {

  my $self = shift;

  if (!%dependent_sources) {

    my $dbi = dbi();
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
  warn( "[DEPRECATED] taxonomy2species_id is a deprecated (unsafe) method. ".
        "Please use species_id2taxonomy instead. Called by ".
        join( ', ', (caller(0))[1..2] ) );

  my $self = shift;

  if (!%taxonomy2species_id) {

    my $dbi = dbi();
    my $sth = $dbi->prepare("SELECT species_id, taxonomy_id FROM species");
    $sth->execute() or croak( $dbi->errstr() );
    while(my @row = $sth->fetchrow_array()) {
      my $species_id = $row[0];
      my $taxonomy_id = $row[1];
      if( my $ori =$taxonomy2species_id{$taxonomy_id} ){
        die( "Taxon $taxonomy_id already used for species $ori. ".
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

    my $dbi = dbi();
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
  warn( "[DEPRECATED] name2species_id is a deprecated (unsafe) method. ".
        "Please use species_id2name instead. Called by ".
        join( ', ', (caller(0))[1..2] ) );

    my $self = shift;

    if ( !%name2species_id ) {

        my $dbi = dbi();
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
            foreach my $name ( split /,\s*/, $row[1] ) {
                if ( my $ori = $name2species_id{$name} ) {
                  die( "Name $name already used for species $ori. ".
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

    my $dbi = dbi();
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
      foreach my $name ( split /,\s*/, $row[1] ) {
        $species_id2name{$species_id} ||= [];
        push @{$species_id2name{$species_id}}, $name;
      }
    }  
  } ## end if ( !%species_id2name)

  return %species_id2name;
} ## end sub species_id2name

# --------------------------------------------------------------------------------
# Update a row in the source table

sub update_source
{
    my ( $dbi, $source_url_id, $checksum, $file_name ) = @_;

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
}


# --------------------------------------------------------------------------------
sub dbi2{

    my $self = shift;
    my ($host, $port, $user, $dbname, $pass) = @_;
    my $dbi2 = undef;

    if ( !defined $dbi2 || !$dbi2->ping() ) {
        my $connect_string =
          sprintf( "dbi:mysql:host=%s;port=%s;database=%s",
            $host, $port, $dbname );

        $dbi2 =
          DBI->connect( $connect_string, $user, $pass,
#            { 'RaiseError' => 1 } )
			)
          or warn( "Can't connect to database: " . $DBI::errstr ) and return undef;
        $dbi2->{'mysql_auto_reconnect'} = 1; # Reconnect on timeout
    }
    
    return $dbi2;
}

# --------------------------------------------------------------------------------

sub dbi
{
    my $self = shift;

    if ( !defined $dbi || !$dbi->ping() ) {
        my $connect_string =
          sprintf( "dbi:mysql:host=%s;port=%s;database=%s",
            $host, $port, $dbname );

        $dbi =
          DBI->connect( $connect_string, $user, $pass,
            { 'RaiseError' => 1 } )
          or croak( "Can't connect to database: " . $DBI::errstr );
        $dbi->{'mysql_auto_reconnect'} = 1; # Reconnect on timeout
    }
    
    return $dbi;
}

# --------------------------------------------------------------------------------

# Compute a checksum of a file.  This checksum is not a straight MD5
# hex digest, but instead the file size combined with the first six
# characters of the MD5 hex digest.  This is to save space.

sub md5sum
{
    my $file = shift;

    if ( !open( FILE, $file ) ) { return undef }
    binmode(FILE);

    my $checksum = sprintf( "%s/%d",
        substr( Digest::MD5->new()->addfile(*FILE)->hexdigest(), 0, 6 ),
        [ stat FILE ]->[7] );

    close(FILE);

    return $checksum;
}

# --------------------------------------------------------------------------------

sub get_xref_id_by_accession_and_source_OLD {

  my ($acc, $source_id, $species_id ) = @_;

  my $dbi = dbi();

  my $sql = '
SELECT xref_id FROM xref WHERE accession=? AND source_id=?';
  if( $species_id ){ $sql .= ' AND species_id=?' }

  my $sth = $dbi->prepare( $sql );

  $sth->execute( $acc, $source_id, ( $species_id ? $species_id : () ) )
    or croak( $dbi->errstr() );

  my @row = $sth->fetchrow_array();
  my $xref_id = $row[0];

  return $xref_id;

}

# --------------------------------------------------------------------------------
# If there was an error, an xref with the same acc & source already exists.
# If so, find its ID, otherwise get ID of xref just inserted

sub insert_or_select {

  my ($self, $sth, $error, $acc, $source, $species) = @_;

  my $id;

  # TODO - check for specific error code rather than for just any error
  if ($error) {

    $id = $self->get_xref($acc, $source, $species);
	
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

}

# --------------------------------------------------------------------------------

sub validate_sources {
  my ($speciesref, @sources) = @_;

  my $dbi = dbi();
  my $sth = $dbi->prepare("SELECT * FROM source WHERE LOWER(name)=?");

  foreach my $source (@sources) {

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

# --------------------------------------------------------------------------------

sub show_valid_sources() {
  my $species = shift;

  my $dbi = dbi();
  my $sth = $dbi->prepare("SELECT distinct(name) FROM source s, source_url su WHERE s.download='Y' and s.source_id = su.source_id and su.species_id = $species");

  $sth->execute();
  while (my @row = $sth->fetchrow_array()) {
    print $row[0] . "\n";
  }

}

# --------------------------------------------------------------------------------

sub validate_species {
  my @species = @_;
  my @species_ids;

  my $dbi = dbi();
  my $sth = $dbi->prepare("SELECT species_id, name FROM species WHERE LOWER(name)=? OR LOWER(aliases) REGEXP ?");
  my ($species_id, $species_name);

  foreach my $sp (@species) {

#    $sth->execute(lc($sp), "%" . lc($sp) . "%");  # no longer allow % as this generates tomany possible errors
    $sth->execute(lc($sp),  "^".lc($sp).",|[ ]".lc($sp)."[,]|^".lc($sp)."\$|[,] ".lc($sp)."\$" );
    $sth->bind_columns(\$species_id, \$species_name);
    if (my @row = $sth->fetchrow_array()) {
      print "Species $sp is valid (name = " . $species_name . ", ID = " . $species_id . ")\n" if($verbose);
      push @species_ids, $species_id;
    } else {
      print STDERR "Species $sp is not valid; valid species are:\n";
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
    print STDERR $row[0] . " (aliases: " . $row[1] . ")\n";
  }

}

sub get_taxonomy_from_species_id{
  my ($self,$species_id) = @_;
  my %hash;

  my $dbi = dbi();
  my $sth = $dbi->prepare("SELECT taxonomy_id FROM species WHERE species_id = $species_id");
  $sth->execute() or croak( $dbi->errstr() );
  while(my @row = $sth->fetchrow_array()) {
    $hash{$row[0]} = 1;
  }   
  $sth->finish;
  return \%hash;
}

sub get_direct_xref{
 my ($self,$stable_id,$type,$link) = @_;

 my $sql = "select general_xref_id from ${type}_direct_xref d where ensembl_stable_id = ?  and linkage_xref= ?";
 my  $direct_sth = $dbi->prepare($sql);

 $direct_sth->execute( $stable_id, $link ) or croak( $dbi->errstr() );
 if(my @row = $direct_sth->fetchrow_array()) {
   return $row[0];
 }
 return undef;
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
  return undef;
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
      dbi->prepare( "INSERT INTO xref "
         . "(accession,version,label,description,source_id,species_id, info_type) "
         . "VALUES(?,?,?,?,?,?,?)" );
  }

  # If the description is more than 255 characters, chop it off and add
  # an indication that it has been truncated to the end of it.

  if ( length($description) > 255 ) {
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

  $direct_xref || die( "Need a direct_xref on which this xref linked too" );
  $acc         || die( "Need an accession of this dependent xref" );
  $version     ||= 0;
  $label       ||= $acc;
  $description ||= undef;
  $linkage     ||= undef;
  $source_id   || die( "Need a source_id for this dependent xref" );
  $species_id  || die( "Need a species_id for this dependent xref" );

  if(!defined($add_xref_sth)){
    $add_xref_sth = dbi->prepare("
INSERT INTO xref 
  (accession,version,label,description,source_id,species_id, info_type)
VALUES
  (?,?,?,?,?,?,?)");
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
}

sub add_to_xrefs{
  my ($self,$master_xref,$acc,$version,$label,$description,$linkage,$source_id,$species_id) = @_;

  $master_xref || die( "Need a master_xref_id on which this xref depends" );
  $acc         || die( "Need an accession of this dependent xref" );
  $version     ||= 0;
  $label       ||= $acc;
  $description ||= undef;
  $linkage     ||= undef;
  $source_id   || die( "Need a source_id for this dependent xref" );
  $species_id  || die( "Need a species_id for this dependent xref" );

  if(!defined($add_xref_sth)){
    $add_xref_sth = dbi->prepare("
INSERT INTO xref 
  (accession,version,label,description,source_id,species_id, info_type)
VALUES
  (?,?,?,?,?,?,?)");
  }
  if(!defined($add_dependent_xref_sth)){
    $add_dependent_xref_sth = dbi->prepare("
INSERT INTO dependent_xref 
  (master_xref_id,dependent_xref_id,linkage_annotation,linkage_source_id)
VALUES
  (?,?,?,?)");
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
      croak(  "Could not find acc $acc in "
            . "xref table source = $source_id of species $species_id\n" );
  }
}


sub add_synonym{
  my ($self, $xref_id, $syn) = @_;

  my $add_synonym_sth;
  
  if(!defined($add_synonym_sth)){
    $add_synonym_sth =  $dbi->prepare("INSERT INTO synonym VALUES(?,?)");
  }

    $add_synonym_sth->execute( $xref_id, $syn )
      or croak( $dbi->errstr() . "\n $xref_id\n $syn\n" );

}



# --------------------------------------------------------------------------------
# Add a single record to the direct_xref table.
# Note that an xref must already have been added to the xref table (xref_id passed as 1st arg)

sub add_direct_xref {

  my ($self, $general_xref_id, $ensembl_stable_id, $ensembl_type, $linkage_type) = @_;

  if (!defined($add_direct_xref_sth{$ensembl_type})){
    my $add_gene_direct_xref_sth = dbi->prepare("INSERT INTO gene_direct_xref VALUES(?,?,?)");
    my $add_tr_direct_xref_sth = dbi->prepare("INSERT INTO transcript_direct_xref VALUES(?,?,?)");
    my $add_tl_direct_xref_sth = dbi->prepare("INSERT INTO translation_direct_xref VALUES(?,?,?)");
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
}

# ------------------------------------------------------------------------------

# Remove potentially problematic characters from string used as file or
# directory names.

sub sanitise {
    my $str = shift;
    $str =~ tr[/:][]d;
    return $str;
}

# ------------------------------------------------------------------------------

# Create database if required. Assumes sql/table.sql and sql/populate_metadata.sql
# are present.

sub create {
  my ( $host, $port, $user, $pass, $dbname, $sql_dir, $drop_db ) = @_;

  my $dbh = DBI->connect( "DBI:mysql:host=$host:port=$port", $user, $pass,
                          {'RaiseError' => 1});

  my $metadata_file =
    catfile( $sql_dir, 'sql', 'populate_metadata.sql' );
  my $ini_file = catfile( $sql_dir, 'xref_config.ini' );

  $| = 1;    # flush stdout

  # Figure out whether to run 'xref_config2sql.pl' or not by comparing
  # the timestamps on 'xref_config.ini' and 'sql/populate_metadata.sql'.
  my $ini_tm  = ( stat $ini_file )[9];
  my $meta_tm = ( stat $metadata_file )[9];

  if ( !defined($meta_tm) || $ini_tm > $meta_tm ) {
    printf( "==> Your copy of 'xref_config.ini' is newer than '%s'\n",
            catfile( 'sql', 'populate_metadata.sql' ) );
    print("==> Should I re-run 'xref_config2sql.pl' for you? [y/N]: ");

    my $reply = <STDIN>;
    chomp $reply;

    if ( lc( substr( $reply, 0, 1 ) ) eq 'y' ) {
      my $cmd = sprintf( "perl %s %s >%s",
                         catfile( $sql_dir, 'xref_config2sql.pl' ),
                         $ini_file, $metadata_file );

      if ( system($cmd) == 0 ) {
        print("==> Done.\n") if($verbose);
      } else {
        if ( $? == -1 ) {
          croak("Failed to execute: $!\n");
        } elsif ( $? & 127 ) {
          croak(
                 sprintf( "Command died with signal %d, %s coredump\n",
                          ( $? & 127 ),
                          ( $? & 128 ) ? 'with' : 'without'
                 ) );
        } else {
          croak( sprintf( "Command exited with value %d\n", $? >> 8 ) );
        }
      }

    }
  } ## end if ( !defined($meta_tm...

  # check to see if the database already exists
  my %dbs = map {$_->[0] => 1} @{$dbh->selectall_arrayref('SHOW DATABASES')};

  if ($dbs{$dbname}) {

    if ( $drop_db ) {     
	$dbh->do( "DROP DATABASE $dbname" );
	print "Database $dbname dropped\n" if($verbose) ; 
    }
  
    if ( $create && !$drop_db ) {
      print "WARNING: about to drop database $dbname on $host:$port; yes to confirm, otherwise exit: ";
      my $p = <STDIN>;
      chomp $p;
      if ($p eq "yes") {
	$dbh->do( "DROP DATABASE $dbname" );
	print "Removed existing database $dbname\n" if($verbose);
      } else {
	print "$dbname NOT removed\n";
	exit(1);
      }
    } elsif ( !$create) {
      croak(  "Database $dbname already exists. "
            . "Use -create option to overwrite it." );
    }
  }

  $dbh->do( "CREATE DATABASE " . $dbname );

  my $table_file = catfile( $sql_dir, 'sql', 'table.sql' );

  printf( "Creating %s from %s\n", $dbname, $table_file ) if($verbose);
  if ( !-e $table_file ) {
    croak( "Cannot open  " . $table_file );
  }

  my $cmd =
    "mysql -u $user -p'$pass' -P $port -h $host $dbname < $table_file";
  system($cmd) == 0
    or croak("Cannot execute the following command (exit $?):\n$cmd\n");

  printf( "Populating metadata in %s from %s\n",
          $dbname, $metadata_file ) if($verbose);
  if ( !-e $metadata_file ) {
    croak( "Cannot open " . $metadata_file );
  }

  $cmd = "mysql -u $user -p'$pass' -P $port -h $host "
    . "$dbname < $metadata_file";
  system($cmd) == 0
    or croak("Cannot execute the following command (exit $?):\n$cmd\n");
}

sub get_label_to_acc{
  my ($self, $name, $species_id, $prio_desc) = @_;
  my %hash1=();

  my $dbi = dbi();
  my $sql = "select xref.accession, xref.label from xref, source where source.name like '$name%' and xref.source_id = source.source_id";
  if(defined($prio_desc)){
    $sql .= " and source.priority_description like '".$prio_desc."'";
  }	
  if(defined($species_id)){
    $sql .= " and xref.species_id  = $species_id";
  }
  my $sub_sth = dbi->prepare($sql);    

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
  my $sub_sth = dbi->prepare($sql);    

  $sub_sth->execute();
  while(my @row = $sub_sth->fetchrow_array()) {
    $hash1{$row[1]} = $row[0];
  }   
 

  return \%hash1;
}


sub get_label_to_desc{
  my ($self, $name, $species_id, $prio_desc) = @_;
  my %hash1=();

  my $dbi = dbi();
  my $sql = "select xref.description, xref.label from xref, source where source.name like '$name%' and xref.source_id = source.source_id";
  if(defined($prio_desc)){
    $sql .= " and source.priority_description like '".$prio_desc."'";
  }	
  if(defined($species_id)){
    $sql .= " and xref.species_id  = $species_id";
  }
  my $sub_sth = dbi->prepare($sql);    

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

sub set_release
{
    my $self = shift;
    my ( $source_id, $release ) = @_;

    my $dbi = dbi();

    my $sth =
      $dbi->prepare(
        "UPDATE source SET source_release=? WHERE source_id=?");

    print "Setting release to '$release' for source ID '$source_id'\n" if($verbose);

    $sth->execute( $release, $source_id );
}

sub get_dependent_mappings {
  my $self = shift;
  my $source_id = shift;
  my $dbi = dbi();
  
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

}

sub get_hgnc_synonyms{
  my $self = shift;
  my %hgnc_syns;
  my %seen;          # can be in more than once fro each type of hgnc.

  my $sql = (<<SYN);
SELECT  x.accession, x.label, sy.synonym 
  FROM xref x, source so, synonym sy
    WHERE x.xref_id = sy.xref_id
      AND so.source_id = x.source_id
      AND so.name like "HGNC"
SYN

  my $dbi = $self->dbi();
  my $sth = $dbi->prepare($sql);    

  $sth->execute;
  my ($acc, $label, $syn);
  $sth->bind_columns(\$acc, \$label, \$syn);

  my $count = 0;
  while($sth->fetch){
    if(!defined($seen{$acc.":".$syn})){
      push @{$hgnc_syns{$acc}}, $syn;
      push @{$hgnc_syns{$label}}, $syn;
      $count++;
    }
    $seen{$acc.":".$syn} = 1;
  }
  $sth->finish;

  return \%hgnc_syns;

}

#
# Store data needed to beable to revert to same stage as after parsing
#

sub parsing_finished_store_data{
  my $self = shift;

# Store max id for 

# gene/transcript/translation_direct_xref     general_xref_id  #Does this change??

# xref                                        xref_id

# dependent_xref                              object_xref_id is all null

# go_xref                                     object_xref_id
# object_xref                                 object_xref_id
# identity_xref                               object_xref_id

  my %table_and_key = ('xref' => "xref_id", 'object_xref' => "object_xref_id");

  my $dbi = $self->dbi();
  foreach my  $table (keys %table_and_key){
#    print "select MAX(".$table_and_key{$table}.") from $table\n";
    my $sth = $dbi->prepare("select MAX(".$table_and_key{$table}.") from $table");
    $sth->execute;
    my $max_val;
    $sth->bind_columns(\$max_val);
    $sth->fetch;
    $self->add_meta_pair("PARSED_".$table_and_key{$table}, $max_val);
    $sth->finish;
  }
  
}



sub reset_to_just_parsed{
  my $self= shift;

 # for dependent_xref set object_xref_id to NULL


 # delete all from gene_transcript_translation, gene/transcript/translation_stable_id, object_xref, identity_xref, go_xref

 # delete from xref where xref_id > MAX


 # set process_status to "parsing_finished"


}


sub reset_to_mapping_finished{
  my $self= shift;
 $self->reset_to_parsed();
#set process_status to "mapping_finished"
}

# --------------------------------------------------------------------------------
1;

