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
											pre_miRNA  => 1,
                      miRNA      => 1,
                      ncRNA      => 1,
                      protein    => 1,
                      pseudogene => 1,
                      rRNA       => 1,
                      snRNA      => 1,
                      snoRNA     => 1,
                      tRNA       => 1 );

# This is some statistics from the 5.4 file 'dmel-all-r5.4.gff.gz',
# looking at the FlyBase (2nd column) object types (3rd column) gene,
# mRNA, miRNA, ncRNA, protein, pseudogene, rRNA, snRNA, snoRNA, and tRNA
# only:
#
# cnt   Dbxref name             source_name in Xref database
# 57884 FlyBase_Annotation_IDs  (special, see below)
# 33441 GB_protein              protein_id
# 14324 GB                      EMBL
# 13951 flygrid                 FlyGrid
# 13265 FlyBase                 flybase_annotation_id
# 12768 dedb                    dedb
# 11745 UniProt/TrEMBL          Uniprot/SPTREMBL
# 10076 INTERPRO                Interpro
# 8077  orthologs               SKIPPED
# 2089  UniProt/Swiss-Prot      Uniprot/SWISSPROT
# 1596  bdgpinsituexpr          bdgpinsituexpr
# 1207  hybrigenics             SKIPPED
# 787   if                      SKIPPED
# 290   MITODROME               SKIPPED
# 153   TF                      SKIPPED
# 82    EPD                     SKIPPED
# 80    MIR                     SKIPPED
# 61    PDB                     SKIPPED
# 56    MEROPS                  SKIPPED
# 17    GCR                     SKIPPED
# 15    Rfam                    SKIPPED
# 12    NRL_3D                  SKIPPED
# 11    GO                      SKIPPED
#
# The Dbxref name 'FlyBase_Annotation_IDs' will be associated with the
# source_names FlyBaseCGID_{gene,transcript,translation} depending on
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
our %source_name_map = ( 'FlyBase'    => 'flybase_annotation_id',
                         'GB'         => 'EMBL',
                         'GB_protein' => 'protein_id',
                         'INTERPRO'   => 'Interpro',
#                         'UniProt/Swiss-Prot' => 'Uniprot/SWISSPROT',
#                         'UniProt/TrEMBL'     => 'Uniprot/SPTREMBL',
                         'bdgpinsituexpr'     => 'BDGP_insitu_expr',
                         'dedb'               => 'DEDb',
                         'flygrid'            => 'FlyGrid',
												 'TF'                 => 'TransFac',
                         'EPD'                => 'EPD',
												 'MIR'                => 'miRBase',
												 'MEROPS'             => 'MEROPS',
												 'BIOGRID'            => 'BioGRID',
												 'FlyReactome'        => 'FlyReactome',
												 'GenomeRNAi_gene'    => 'GenomeRNAi',
												 'INTERACTIVEFLY'     => 'InteractiveFly',
												 'MITODROME'          => 'MitoDrome',
												 'flyexpress'         => 'FlyExpress',
												 'Rfam'               => 'RFAM',
			 #'FlyAtlas'           => 'FlyAtlas',
			 #'GCR'                => 'GPCR',
			 #'GLEANR'             => 'GLEAN-R',
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

  if ( !defined( $source_id{$source_name} )
       || $source_id{$source_name} < 0 )
  {
    carp(
       sprintf( "Can not find source_id for source '%s'", $source_name )
    );
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


	# Create a go source id for GO terms extracted from the GFF file
	# The reason being to separate from any other sources for GO terms
	# like dependent xrefs (GOA UniProt).

	my $go_source_id = $self->SUPER::get_source_id_for_source_name('GO','flybasego');

  print "-------------------------\n";
  print "FlybaseParser::run species_id $species_id\n";
  print "-------------------------\n\n";

  my $data_file = @{$files}[0];

  # Fetch hashes of already stored Uniprot and Interpro accessions.
  my %pre_xref_ids = (
        'Uniprot'  => $self->get_valid_codes( 'uniprot',  $species_id ),
        'Interpro' => $self->get_valid_codes( 'interpro', $species_id )
  );

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
    else                                    { $type = 'unknown' }

    # For the 'Dbxref' and 'Ontology_term' attributes, split them up on
    # commas, divide into key-value pairs, and store them.
    foreach my $attribute_key ( 'Dbxref', 'Ontology_term' ) {
      if ( exists( $attributes{$attribute_key} ) ) {
        my %tmphash;
        foreach
          my $subattribute ( split( /,/, $attributes{$attribute_key} ) )
        {
						# For GO term, we keep the form GO:0004080
						if ($subattribute =~ /^GO/) {
								#print "$attribute_key Storing GO term: $subattribute for $id\n";
								push( @{ $tmphash{'GO'} }, $subattribute );
						} else {
								my ( $key, $value ) = split( /:/, $subattribute, 2 );
								push( @{ $tmphash{$key} }, $value );
						}
        }

        # Replace the attribute entry with the hash.
        $attributes{$attribute_key} = \%tmphash;
      }
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
		

    my $dbxref = $attributes{'Dbxref'};

    #-------------------------------------------------------------------
    # Store Xrefs and Direct Xrefs for all the interesting Dbxref
    # entries.
    #-------------------------------------------------------------------
    foreach my $dbxref_name ( keys( %{$dbxref} ) ) {
      if ( exists( $source_name_map{$dbxref_name} ) ) {
        my $source_name = $source_name_map{$dbxref_name};
        my $source_id =
          $self->get_source_id_for_source_name($source_name);

        # Treat Uniprot and Interpro differently.
        my ($pre_source) = ( $source_name =~ /^(Uniprot|Interpro)/ );

        if ( defined($pre_source) ) {
          foreach my $accession ( @{ $dbxref->{$dbxref_name} } ) {
            if ( exists( $pre_xref_ids{$pre_source}{$accession} ) ) {
	      foreach my $xref_id (@{ $pre_xref_ids{$pre_source}{$accession} }){
		$self->add_direct_xref($xref_id, $id, $type, '' );
		$xref_ids{$pre_source}{$accession} = $xref_id;
	      }	
            } else {
              $xref_ids{ $pre_source . ' (missed)' }{$accession} = -1;
            }
          }
        } else {
          foreach my $accession ( @{ $dbxref->{$dbxref_name} } ) {
            my $xref_id;
            if ( exists( $xref_ids{$source_name}{$accession} ) ) {
              $xref_id = $xref_ids{$source_name}{$accession};
            } else {
              # The Dbxref 'bdgpinsituexpr' needs case sensitivity, just
              # like the FlyBase Names, so use the ID as the accession
              # for this source.
              if ( $dbxref_name eq 'bdgpinsituexpr' ) {
                $xref_id =
                  $self->add_xref({ acc        => $id,
				    label      => $accession,
				    source_id  => $source_id,
				    species_id => $species_id,
				    info_type  =>'DIRECT'} );
              } else {
                $xref_id =
                  $self->add_xref({ acc        => $accession,
				    label      => $accession,
				    source_id  => $source_id,
				    species_id => $species_id,
				    info_type  => 'DIRECT'} );
              }
              $xref_ids{$source_name}{$accession} = $xref_id;
					}

            $self->add_direct_xref( $xref_id, $id, $type, '' );
          }
        }
      } ## end if ( exists( $source_name_map...
    } ## end foreach my $dbxref_name ( keys...

    #-------------------------------------------------------------------
    # Store Xrefs and Direct Xrefs for the GO 'Ontology_term' entries.
    #-------------------------------------------------------------------

    if ( exists( $attributes{'Ontology_term'}{'GO'} ) ) {
      my $source_name = 'GO';
      #my $source_id =
      #  $self->get_source_id_for_source_name($source_name);


      foreach my $accession ( @{ $attributes{'Ontology_term'}{'GO'} } )
      {
					my $xref_id;
					if ( exists( $xref_ids{$source_name}{$accession} ) ) {
							$xref_id = $xref_ids{$source_name}{$accession};
					} else {
							#print "FlyBaseParser\t$source_name\tadd_xref: $accession, $go_source_id\n";  
							$xref_id =
									$self->add_xref({ acc        => $accession,
																		label      => $accession,
																		source_id  => $go_source_id,
																		species_id => $species_id,
																		info_type  => 'DIRECT'} );
							$xref_ids{$source_name}{$accession} = $xref_id;
					}
					
					#print "Add GO direct xref for $id with $type \n";
					$self->add_direct_xref( $xref_id, $id, $type, 'IEA' );
      }
    }

    #-------------------------------------------------------------------
    # Store Xrefs and Direct Xrefs for the 'FlyBase_Annotation_IDs'
    # Dbxref entry (depends on type of 'ID').
    #-------------------------------------------------------------------
    if ( exists( $dbxref->{'FlyBase_Annotation_IDs'} ) ) {
      my $source_name = $special_source_name_map{$type}{'Dbxref'};
      my $source_id =
        $self->get_source_id_for_source_name($source_name);

      foreach my $accession ( @{ $dbxref->{'FlyBase_Annotation_IDs'} } )
      {
        my $xref_id;
        if ( exists( $xref_ids{$source_name}{$accession} ) ) {
          $xref_id = $xref_ids{$source_name}{$accession};
        } else {
          $xref_id =
            $self->add_xref({ acc        => $accession,
			      label      => $accession,
			      source_id  => $source_id,
			      species_id => $species_id,
			      info_type  => 'DIRECT'} );
          $xref_ids{$source_name}{$accession} = $xref_id;
        }

        $self->add_direct_xref( $xref_id, $id, $type, '' );
      }

    }

    #-------------------------------------------------------------------
    # Store Xref and Direct Xref for the 'Name' (depends on type of
    # 'ID').
    #-------------------------------------------------------------------
    {
      my $source_name = $special_source_name_map{$type}{'Name'};
      my $source_id =
        $self->get_source_id_for_source_name($source_name);

      my $accession = $attributes{'Name'};

			# every fly Names but d. melanogaster ones start with D...\ (like Dper\β3galt6)
			# we remove the prefix. 
			if ($accession =~ m/D...\\(.+)/) {
					$accession = $1;
			}

			my $description = (defined($attributes{'fullname'})) ? $attributes{'fullname'} : '';

			# because FlyBase use %2C to distinguish from the , separator in the GFF dump
			# we have to put it back
			$description =~ s/%2C/,/g; 

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
			    info_type => 'DIRECT'} );

        $xref_ids{$source_name}{$accession} = $xref_id;
      }

      $self->add_direct_xref( $xref_id, $id, $type, '' );
    }

    #-------------------------------------------------------------------
    # Store Xref and Direct Xref for the 'ID' (depends on type of 'ID').
    #-------------------------------------------------------------------
    {
      my $source_name = $special_source_name_map{$type}{'ID'};
      my $source_id =
        $self->get_source_id_for_source_name($source_name);

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
			    onfo_type  => 'DIRECT'} );
        $xref_ids{$source_name}{$accession} = $xref_id;
      }

      $self->add_direct_xref( $xref_id, $id, $type, '' );
 


			#-------------------------------------------------------------------
			# Now, if we have aliases for this gene/transcript/translation
			# Store them in the external_synonym table.
			#-------------------------------------------------------------------
			
			if (defined ($attributes{$alias_key})) {
					foreach my $alias (@{$attributes{$alias_key}}) {
							$self->add_synonym($xref_id, $alias);
					}
			}
	  }
	
  } ## end while ( defined( my $line...
  $data_io->close();

  alarm(0);

  print("FlybaseParser Summary:\n") if($verbose);
	print("--------------------------------------------------------------\n") if($verbose);
  foreach my $label ( sort( keys(%xref_ids) ) ) {
    my $accessions = $xref_ids{$label};
    printf( "\t%-32s %6d\n", $label, scalar( keys( %{$accessions} ) ) ) if($verbose);
  }
	print("--------------------------------------------------------------\n") if($verbose);



  return 0;
} ## end sub run

1;
