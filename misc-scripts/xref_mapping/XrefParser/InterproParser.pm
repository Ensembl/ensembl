package XrefParser::InterproParser;
  
use strict;
use POSIX qw(strftime);
use File::Basename;
  
use base qw( XrefParser::BaseParser );
  
my $xref_sth ;
my $dep_sth;
  
 
  
# --------------------------------------------------------------------------------
# Parse command line and run if being run directly
  
if (!defined(caller())) {
  
  if (scalar(@ARGV) != 1) {
    print "\nUsage: InterproParser.pm file <source_id> <species_id>\n\n";
    exit(1);
  }
  
  run(@ARGV);
}
  
 
sub run {
  my $self = shift if (defined(caller(1)));

  my $source_id = shift;
  my $species_id = shift;
  my $files_ref = shift;
  my $release_file = shift;
  my $verbose = shift;


  my $file = @{$files_ref}[0];

  if(!defined($source_id)){
    $source_id = $self->get_source_id_for_filename($file);
  }
  if(!defined($species_id)){
    $species_id = $self->get_species_id_for_filename($file);
  }

  my $add_interpro_sth =
    $self->dbi()
    ->prepare("INSERT INTO interpro (interpro, pfam, dbtype ) VALUES(?,?,?)");

  my $get_interpro_sth =
    $self->dbi()
    ->prepare( "SELECT interpro FROM interpro "
        . "WHERE interpro = ? AND pfam = ?" );

  my $add_xref_sth =
    $self->dbi()
    ->prepare( "INSERT INTO xref "
        . "(accession,version,label,description,source_id,species_id, info_type) "
        . "VALUES(?,?,?,?,?,?,?)" );

  my $dir = dirname($file);

  my %short_name;
  my %description;
  my %pfam;

  my $xml_io = $self->get_filehandle($file);

  if ( !defined $xml_io ) {
    print "ERROR: Can't open hugo interpro file $file\n";
    return 1;    # 1= error
  }

  #<interpro id="IPR001023" type="Family" short_name="Hsp70" protein_count="1556">
  #    <name>Heat shock protein Hsp70</name>
  #     <db_xref protein_count="18" db="PFAM" dbkey="PF01278" name="Omptin" />
  #      <db_xref protein_count="344" db="TIGRFAMs" dbkey="TIGR00099" name="Cof-subfamily" />
  
  my %count;
  local $/ = "</interpro>";

  my $last = "";
  my $i =0;

  while ( $_ = $xml_io->getline() ) {
    my ($interpro)   = $_ =~ /interpro id="(\S+)"/;
    my ($short_name) = $_ =~ /short_name="(\S+)"/;
    my ($name)       = $_ =~ /<name>(.*)<\/name>/;

    if ($interpro) {
        #      print $interpro."\n";
        if ( !$self->get_xref( $interpro, $source_id, $species_id ) ) {
            $count{INTERPRO}++;
            if (
                !$add_xref_sth->execute(
                    $interpro, '',         $short_name,
                    $name,     $source_id, $species_id, 'MISC'
                )
              )
            {
                print STDERR "Problem adding '$interpro'\n";
                return 1;    # 1 is an error
            }
        }

        my ($members) = $_ =~ /<member_list>(.+)<\/member_list>/s;

        while ( $members =~
/db="(PROSITE|PFAM|PREFILE|PROFILE|TIGRFAMs|PRINTS|PIRSF|SMART|SSF)"\s+dbkey="(\S+)"/cgm
          )
        {
            my ( $db_type, $id ) = ( $1, $2 );
            if( $db_type eq 'SSF' ){ $id =~ s/^SSF// } # Strip SSF prefix

#            if ( !$self->get_xref( $interpro, $id, $species_id ) ) {#no idea why if was checking for something with source of the name it would never work
	    $add_interpro_sth->execute( $interpro, $id, $db_type );
	    $count{$db_type}++;
#            }
        }
    }
  }

  $xml_io->close();

    for my $db ( keys %count ) {
        print "\t" . $count{$db} . " $db loaded.\n" if($verbose);
    }

    if ( defined $release_file ) {
        # Parse the second file that we got.  This is assumed to be the
        # HTML file that will contain the release information.

        my $release;
        my $release_io = $self->get_filehandle($release_file);
        while ( defined( my $line = $release_io->getline() ) ) {
            chomp $line;
            if ( $line =~ m#(Release [0-9.]+, .*)# ) {
                $release = $1;
                last;
            }
        }
        $release_io->close();

        if ( defined $release ) {
            print "Interpro release is '$release'\n" if($verbose);
            $self->set_release( $source_id, $release );
        } else {
            print "Did not find release info in '$release_file'\n" if($verbose);
        }
    }

    return 0;
}

1;
