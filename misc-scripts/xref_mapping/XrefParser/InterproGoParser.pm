package XrefParser::InterproGoParser;

use strict;
use XrefParser::BaseParser;
use vars qw(@ISA);
@ISA = qw(XrefParser::BaseParser);

sub run {
  my $self = shift if (defined(caller(1)));

  my $source_id = shift;
  my $species_id = shift;
  my $files       = shift;
  my $release_file   = shift;
  my $verbose       = shift;

  my $file = @{$files}[0];

  my $file_io = $self->get_filehandle($file)
      || ( print( "ERROR: Cannot open $file\n" ) && return 1 );


  my %interpros = %{$self->get_valid_codes("interpro",$species_id)};

  scalar( keys %interpros )
      || ( print STDERR  "ERROR: No InterPro xrefs found in DB"  && return 1 );


  #get the "main" GO source id.
  $source_id = $self->get_source_id_for_source_name("GO","main");

  # get the mapping that are already there so that we don't get lots of duplicates.
  # stored in the global hash xref_dependent_mapped.
  $self->get_dependent_mappings($source_id);


  # Process the file
  my( $skip_count, $dependent_xref_count ) = (0,0);
  while( my $line = $file_io->getline() ){
    next if $line =~ /^!/; # Skip comments

    # Example line
    # InterPro:IPR000003 Retinoid X receptor > GO:DNA binding ; GO:0003677
    if( $line =~ m/^InterPro:(\S+)\s+(.+)\s+>\s+GO:(.+)\s+;\s+(GO:\d+)/ ){
      my $ipro_id = $1;
      my $go_desc = $2;
      my $go_term = $3;
      my $go_id   = $4;
      
      if(defined($interpros{$ipro_id})){
	$self->add_to_xrefs($interpros{$ipro_id},$go_id,1,$go_id,$go_term,'IEA',
			    $source_id,$species_id);
	$dependent_xref_count++;
      }
      else{
	$skip_count++;
      }
    }

  }    
  print "Parsed identifiers from $file\n".
        "\tadded $dependent_xref_count GO xrefs dependent on InterPro\n".
        "\tskipped $skip_count GO terms due to missing InterPros\n" if($verbose);

  return 0;
}

sub get_valid_codes{

  my ($self, $source_name,$species_id) =@_;

  # First cache synonyms so we can quickly add them later
  my %synonyms;
  my $syn_sth = $self->dbi()->prepare("SELECT xref_id, synonym FROM synonym");
  $syn_sth->execute();

  my ($xref_id, $synonym);
  $syn_sth->bind_columns(\$xref_id, \$synonym);
  while ($syn_sth->fetch()) {

    push @{$synonyms{$xref_id}}, $synonym;

  }

  my %valid_codes;
  my @sources;

  my $sql = "select source_id from source where upper(name) like '%".uc($source_name)."%'";
  my $sth = $self->dbi()->prepare($sql);
  $sth->execute();
  while(my @row = $sth->fetchrow_array()){
    push @sources,$row[0];
  }
  $sth->finish;

  foreach my $source (@sources){
    $sql = 'select x.accession, x.xref_id from xref x, interpro i where i.interpro = x.accession and x.species_id = '.$species_id.
           ' and x.source_id = '.$source.' and i.dbtype ne "PRINTS"';
    my $sth = $self->dbi()->prepare($sql);
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

sub new{
  my $self = {};
  bless $self, __PACKAGE__;
  return $self;
}
