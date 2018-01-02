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

# $Id$

package XrefParser::FlybaseParser;

use strict;
use warnings;

use Carp;

use base qw( XrefParser::BaseParser );
my $verbose;

# The object types we'd like to parse.
our %object_types = ( gene       => 1,
                      mRNA       => 1,
                      protein    => 1,
                      pseudogene => 1,
                      miRNA      => 1,
                      ncRNA      => 1,
                      pre_miRNA  => 1,
                      rRNA       => 1,
                      snoRNA     => 1,
                      snRNA      => 1,
                      tRNA       => 1 );

# The Dbxref name 'flybase_annotation_id' will be associated with the
# source_name FlyBaseCGID_{gene,transcript,translation} depending on
# the type of 'ID' of the line.
#
# Likewise, the source_names FlyBaseName_{gene,transcript,translation}
# will be associated with the 'Name' of each entry depending on the type
# of 'ID'.
#
# ... and the source_names flybase_{gene,transcript,translation}_id will
# be associated with the 'ID' of each entry depending on the type of
# 'ID'.

# This hash will translate the Dbxref names in the data file into source
# names known by the Xref system.
# The protein-level FlyBase annotations (UniProt, SwissProt, Interpro) are
# not imported; they're attached to genes, which means that when
# we run our protein pipeline (at the translation level), those results are
# shifted to the gene level, which messes up the web display.
our %source_name_map = ( 'FlyBase'            => 'flybase_annotation_id',
                         'BIOGRID'            => 'BioGRID',
                         'EPD'                => 'EPD',
                         'flyexpress'         => 'FlyExpress',
                         'FlyReactome'        => 'FlyReactome',
                         'GenomeRNAi'         => 'GenomeRNAi',
                         'INTERACTIVEFLY'     => 'InteractiveFly',
                         'MIR'                => 'miRBase',
                         'MITODROME'          => 'MitoDrome',
                         'Rfam'               => 'Rfam',
                         'TF'                 => 'TransFac',
);

# This is for source_ids that depend on the type of 'ID' of the line.
our %special_source_name_map = (
                                'gene' => {
                                   'Dbxref' => 'FlyBaseCGID_gene',
                                   'Name'   => 'FlyBaseName_gene',
                                   'ID'     => 'flybase_gene_id'
                                },
                                'transcript' => {
                                   'Dbxref' => 'FlyBaseCGID_transcript',
                                   'Name'   => 'FlyBaseName_transcript',
                                   'ID'     => 'flybase_transcript_id'
                                },
                                'translation' => {
                                  'Dbxref' => 'FlyBaseCGID_translation',
                                  'Name'   => 'FlyBaseName_translation',
                                  'ID'     => 'flybase_translation_id'
                                } );

# This hash will eventually be populated with the source_id for the
# sources above.
our %source_id;

sub get_source_id_for_source_name {
  my ($self, $source_name, $priority_desc) = @_;

  if ( !defined( $source_id{$source_name} ) ) {
    $source_id{$source_name} =
      $self->SUPER::get_source_id_for_source_name($source_name, $priority_desc);

    printf( "source_id for source '%s' is %d\n",
      $source_name, $source_id{$source_name} ) if ($verbose);
  }

  if ( !defined( $source_id{$source_name} ) || $source_id{$source_name} < 0 )
  {
    carp( sprintf( "Can not find source_id for source '%s'", $source_name ) );
  }

  return $source_id{$source_name};
}

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

  # Note: The import of the GO terms from the FlyBase GFF has been removed.
  # Only Dmel is annotated with evidence codes by FlyBase, the other flies
  # are inferred from Interpro analysis - so can be handled equally well by
  # the GOParser (which maps them to translations rather than genes too).
  # In addition, the evidence codes for Dmel are not even in the GFF
  # file, and have to be patched across further down the line. A new Dmel-
  # specific section has been added to GOParser to automate this, in the same
  # way that C. elegans is done, for example.

  print "-------------------------\n";
  print "FlybaseParser::run species_id $species_id\n";
  print "-------------------------\n\n";

  my $data_file = @{$files}[0];

  my %xref_ids;

  my $data_io = $self->get_filehandle($data_file);

  my ( $count_read, $count_skipped, $last_count_read ) = ( 0, 0, 0 );

  my $status_interval = 30;
  local $SIG{ALRM} = sub {
    printf( "%d lines read, %d skipped, %d parsed; %d lines/s\n",
            $count_read, $count_skipped,
            $count_read - $count_skipped,
            ( $count_read - $last_count_read )/$status_interval ) if($verbose);
    $last_count_read = $count_read;
    alarm($status_interval);
  };
  alarm($status_interval);

  while ( defined( my $line = $data_io->getline() ) ) {
    ++$count_read;

    # Skip comment lines at the start of the file.
    if ( substr( $line, 0, 1 ) eq '#' ) { ++$count_skipped; next }

    chomp($line);

    # Split each line into fields.
    my @fields = split( /\t/, $line );

    # Only pick out the interesting lines.
    if (
         !(    defined( $fields[1] )
            && $fields[1] eq 'FlyBase'
            && defined( $fields[2] )
            && exists( $object_types{ $fields[2] } ) ) )
    {
      ++$count_skipped;
      next;
    }

    # Go though each attribute (from the 9th field), split them up into
    # key-value pairs and store them.
    my %attributes;
    foreach my $attribute ( split( /;/, $fields[8] ) ) {
      my ( $key, $value ) = split( /=/, $attribute );
      if ( $key ne '' && $value ne '' ) {
        $attributes{$key} = $value;
      }
    }

    my $id = $attributes{'ID'};

    my $type;
    if    ( substr( $id, 0, 4 ) eq 'FBgn' ) { $type = 'gene' }
    elsif ( substr( $id, 0, 4 ) eq 'FBtr' ) { $type = 'transcript' }
    elsif ( substr( $id, 0, 4 ) eq 'FBpp' ) { $type = 'translation' }
    else                                     { $type = 'unknown' }

    if ( exists( $attributes{'Dbxref'} ) ) {
      my %tmphash;
      foreach my $subattribute ( split( /,/, $attributes{'Dbxref'} ) ) {
        my ( $key, $value ) = split( /:/, $subattribute, 2 );
        push( @{ $tmphash{$key} }, $value );
      }

      # Replace the attribute entry with the hash.
      $attributes{'Dbxref'} = \%tmphash;
    }

    # For the 'Alias' attributes, we split them up by commas
    # but we can't divide them in to key-value. So, we'll create
    # a fake key Alias.
    # Aliases will be stored as synonyms and will comprise secondary
    # IDs from FlyBase to keep tracks of split/merged annotations.
    my $alias_key = 'Alias';
    if ( exists( $attributes{$alias_key} ) ) {
      my @tmp_array = split( /,/, $attributes{$alias_key} );
      $attributes{$alias_key} =\@tmp_array;
    }

    #----------------------------------------------------------------------
    # Store Xrefs and Direct Xrefs for all the interesting Dbxref entries.
    #----------------------------------------------------------------------
    my $dbxref = $attributes{'Dbxref'};
    foreach my $dbxref_name ( keys( %{$dbxref} ) ) {
      if ( exists( $source_name_map{$dbxref_name} ) ) {
        my $source_name = $source_name_map{$dbxref_name};
        my $source_id = $self->get_source_id_for_source_name($source_name);

        foreach my $accession ( @{ $dbxref->{$dbxref_name} } ) {
          my $xref_id;
          if ( exists( $xref_ids{$source_name}{$accession} ) ) {
            $xref_id = $xref_ids{$source_name}{$accession};
          } else {
            $xref_id =
              $self->add_xref({ acc        => $accession,
                                label      => $accession,
                                source_id  => $source_id,
                                species_id => $species_id,
                                info_type  => 'DIRECT'}
            );
            $xref_ids{$source_name}{$accession} = $xref_id;
          }
          $self->add_direct_xref( $xref_id, $id, $type, '' );
        }
      }
    }

    #-------------------------------------------------------------------
    # Store Xrefs and Direct Xrefs for the 'FlyBase_Annotation_IDs'
    # Dbxref entry (depends on type of 'ID').
    #-------------------------------------------------------------------
    if ( exists( $dbxref->{'FlyBase_Annotation_IDs'} ) ) {
      my $source_name = $special_source_name_map{$type}{'Dbxref'};
      my $source_id = $self->get_source_id_for_source_name($source_name);

      foreach my $accession ( @{ $dbxref->{'FlyBase_Annotation_IDs'} } ) {
        my $xref_id;
        if ( exists( $xref_ids{$source_name}{$accession} ) ) {
          $xref_id = $xref_ids{$source_name}{$accession};
        } else {
          $xref_id =
            $self->add_xref({ acc        => $accession,
                              label      => $accession,
                              source_id  => $source_id,
                              species_id => $species_id,
                              info_type  => 'DIRECT'}
          );
          $xref_ids{$source_name}{$accession} = $xref_id;
        }
        $self->add_direct_xref( $xref_id, $id, $type, '' );
      }
    }

    #----------------------------------------------------------------------
    # Store Xref and Direct Xref for the 'Name' (depends on type of 'ID').
    #----------------------------------------------------------------------
    {
      my $source_name = $special_source_name_map{$type}{'Name'};
      my $source_id = $self->get_source_id_for_source_name($source_name);

      my $accession = $attributes{'Name'};

			# Names other than D. melanogaster start with D...\ (like Dper\β3galt6)
			$accession =~ s/^D...\\//;

			my $description = (defined($attributes{'fullname'})) ? $attributes{'fullname'} : '';

			# FlyBase use %2C to distinguish from the , separator in the GFF dump;
			# we have to put it back
			$description =~ s/%2C/,/g;

			# Embedded newlines wreak havoc further down the line
			$description =~ s/[\n\r]//gm;
			# And slashes to ensure that slashes aren't mistakenly interpreted as control characters
			$description =~ s/\\/\\\\/gm;

      my $xref_id;

      if ( exists( $xref_ids{$source_name}{$accession} ) ) {
        $xref_id = $xref_ids{$source_name}{$accession};
      } else {
        $xref_id =
          $self->add_xref({ acc =>  $id,
                            label => $accession,
                            desc => $description,
                            source_id => $source_id,
                            species_id => $species_id,
                            info_type => 'DIRECT'}
        );
        $xref_ids{$source_name}{$accession} = $xref_id;
      }
      $self->add_direct_xref( $xref_id, $id, $type, '' );
    }

    #-------------------------------------------------------------------
    # Store Xref and Direct Xref for the 'ID' (depends on type of 'ID').
    #-------------------------------------------------------------------
    {
      my $source_name = $special_source_name_map{$type}{'ID'};
      my $source_id = $self->get_source_id_for_source_name($source_name);

      my $accession = $id;
      my $xref_id;

      if ( exists( $xref_ids{$source_name}{$accession} ) ) {
        $xref_id = $xref_ids{$source_name}{$accession};
      } else {
        $xref_id =
          $self->add_xref({ acc        => $accession,
                            label      => $accession,
                            source_id  => $source_id,
                            species_id => $species_id,
                            info_type  => 'DIRECT'}
        );
        $xref_ids{$source_name}{$accession} = $xref_id;
      }
      $self->add_direct_xref( $xref_id, $id, $type, '' );

			#-------------------------------------------------------------------
			# Now, if we have aliases for this gene/transcript/translation
			# Store them in the external_synonym table.
			#-------------------------------------------------------------------

			if (defined ($attributes{$alias_key})) {
        foreach my $alias (@{$attributes{$alias_key}}) {
          # Skip synonyms with non-ASCII characters
          next unless $alias =~ /^[\x00-\x7F]+$/;
          # Embedded newlines wreak havoc further down the line
          $alias =~ s/[\n\r]//gm;
          $self->add_synonym($xref_id, $alias);
        }
			}
	  }

  }
  $data_io->close();

  alarm(0);

  if ($verbose) {
    print("FlybaseParser Summary:\n");
    print("--------------------------------------------------------------\n");
    foreach my $label ( sort( keys(%xref_ids) ) ) {
      my $accessions = $xref_ids{$label};
      printf( "\t%-32s %6d\n", $label, scalar( keys( %{$accessions} ) ) );
    }
    print("--------------------------------------------------------------\n");
  }

  return 0;
}

1;
