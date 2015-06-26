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

package XrefParser::InterproParser;
  
use strict;
use warnings;
use POSIX qw(strftime);
use File::Basename;
use Carp;
use base qw( XrefParser::BaseParser );

sub run {
  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $files        = $ref_arg->{files};
  my $release_file = $ref_arg->{rel_file};
  my $verbose      = $ref_arg->{verbose};

  if((!defined $source_id) or (!defined $species_id) or (!defined $files) or (!defined $release_file)){
    croak "Need to pass source_id, species_id, files and rel_file as pairs";
  }
  $verbose |=0;

  my $file = @{$files}[0];

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
        . "(accession,label,description,source_id,species_id, info_type) "
        . "VALUES(?,?,?,?,?,?)" );

  my $dir = dirname($file);

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
                    $interpro, $short_name,
                    $name, $source_id, $species_id, 'MISC'
                )
              )
            {
                print STDERR "Problem adding '$interpro'\n";
                return 1;    # 1 is an error
            }
        }

        my ($members) = $_ =~ /<member_list>(.+)<\/member_list>/s;

        while ( $members =~
/db="(PROSITE|PFAM|PREFILE|PROFILE|TIGRFAMs|PRINTS|PIRSF|SMART|SSF|Gene3D|Panther)"\s+dbkey="(\S+)"/cgm
          )
        {
            my ( $db_type, $id ) = ( $1, $2 );

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
