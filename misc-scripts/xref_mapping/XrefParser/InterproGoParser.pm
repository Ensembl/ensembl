=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

package XrefParser::InterproGoParser;

use strict;
use Carp;
use XrefParser::BaseParser;
use vars qw(@ISA);
@ISA = qw(XrefParser::BaseParser);

sub run {
  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $files        = $ref_arg->{files};
  my $verbose      = $ref_arg->{verbose};

  if((!defined $source_id) or (!defined $species_id) or (!defined $files) ){
    croak "Need to pass source_id, species_id and files as pairs";
  }
  $verbose |=0;

  my $file = @{$files}[0];

  my $file_io = $self->get_filehandle($file);
  unless ($file_io) {
      print( "ERROR: Cannot open $file\n" ); 
      return 1 ;
  }

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
	$self->add_dependent_xref({ master_xref_id => $interpros{$ipro_id},
				    acc            => $go_id,
				    version        => 1,
				    label          => $go_id,
				    desc           => $go_term,
				    linkage        => 'IEA',
				    source_id      => $source_id,
				    species_id     => $species_id} );
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
           ' and x.source_id = '.$source.' and i.dbtype != "PRINTS"';
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

