package Bio::EnsEMBL::Utils::SchemaConversion;

=head1 NAME

Bio::EnsEMBL::Utils::SchemaConversion - Utility module for Vega schema conversion script

=head1 SYNOPSIS

    my $serverroot = '/path/to/ensembl';
    my $conversion = Bio::EnsEMBL::Utils::ConversionSupport->new($serverroot);

    # parse common options
    $conversion->conv_usage->parse_common_options;

    # convert from schema 19 to 20+
    $conversion->do_conversion()


=head1 DESCRIPTION

This module is a helper module for database conversion, for both vega-vega and 
ensembl-vega schemas. It provides a wrapper around SeqStoreConverter::BasicConverter
and the species specific methods therein. Also provides access to helper functions in
Bio::EnsEMBL::Utils::ConversionSupport

=head1 LICENCE

This code is distributed under an Apache style licence:
Please see http://www.ensembl.org/code_licence.html for details

=head1 AUTHOR

Steve Trevanion <st3@sanger.ac.uk>

=head1 CONTACT

Post questions to the EnsEMBL development list ensembl-dev@ebi.ac.uk

=cut

use  Bio::EnsEMBL::Utils::ConversionSupport;
use strict;
use warnings;

use Data::Dumper;

=head2 new

  Example     : $conversion->Bio::EnsEMBL::Utils::SchemaConversion->new($serverroot);
  Description : Constructor, including an instance of a Bio::EnsEMBL::Utils::ConversionSupport object
  Return type : Bio::EnsEMBL::Utils::SchemaConversion object 
  Exceptions  : thrown if $Siteroot not passed over
  Caller      : $Siteroot/utils/vega_schema_conversion

=cut

sub new {
	my $class = shift;
	my $support = shift;
	my $self = {};
	bless ($self,$class);
	$self->{config} =  Bio::EnsEMBL::Utils::ConversionSupport->new($support);
	return $self;
}

=head2 conv_support

  Example     : $conversion->conv_support; 
  Description : This method provides access to Bio::EnsEMBL::Utils::ConversionSupport methods
  Return type : Bio::EnsEMBL::Utils::ConversionSuppor object
  Exceptions  : none
  Caller      : general

=cut

sub conv_support {
	my $self = shift;
	return $self->{config};
}

=head2 species_alias

  Example     : $self->species_alias
  Description : lookup table for common species names and schema conversion modules
  Return type : string
  Exceptions  : die if wrong species name used
  Caller      : $self

=cut

sub species_alias {
	my $self=shift;
	my $name = shift;
	my %aliases = (
				   dog       => 'CanisFamiliaris',
				   human     => 'HomoSapiens',
				   mouse     => 'MusMusculus',
				   zebrafish => 'DanioRerio'
			  );
	return $aliases{$name} || die "invalid species name, please check configuration file";
}

=head2 choose_conversion_type

  Example     : $conversion->choose_conversion_type
  Description : compares conversion type (ensembl or vega) and species type with 
                available modules and chooses that to use for the conversion. Stores
                a converter object within the caller
  Return type : none
  Exceptions  : none
  Caller      : $Siteroot/utils/vega_schema_conversion

=cut

sub choose_conversion_type {
	my $self = shift;
	my $siteroot = $self->conv_support->serverroot;
	my $converter;
	my $species;
	$species = $self->species_alias($self->conv_support->param('species'));
	if ($self->conv_support->param('do_vega_schema_conversion')) {
		$species = "vega::".$species;
		eval "require SeqStoreConverter::$species";
		if($@) {
			warn("Could not require conversion module SeqStoreConverter::$species\ for vega conversion\n" .
				 "Using SeqStoreConverter::BasicConverter instead:\n$@");
			require SeqStoreConverter::BasicConverter;
			$species = "BasicConverter";
		}
		else {
			warn "Using conversion module SeqStoreConverter::$species for vega conversion\n";
		}
	}
	else {
		eval "require SeqStoreConverter::$species";
		if($@) {
			warn("Could not require conversion module SeqStoreConverter::$species for Ensembl conversion\n" .
				 "Using SeqStoreConverter::BasicConverter instead:\n$@");
			require SeqStoreConverter::BasicConverter;
			$species = "BasicConverter";
		}
		else {
			warn "Using conversion module SeqStoreConverter::$species for Ensembl conversion\n";
		}
	}
	$converter = "SeqStoreConverter::$species"->new
		( $self->conv_support->param('user'), 
		  $self->conv_support->param('pass'), 
		  $self->conv_support->param('host'), 
		  $self->conv_support->param('source_db'), 
		  $self->conv_support->param('target_db'), 
		  $siteroot.$self->conv_support->param('core_sql'), 
		  $siteroot.$self->conv_support->param('vega_sql'), 
		  $self->conv_support->param('force'), 
		  $self->conv_support->param('verbose'),
		  '',
		);
	
	$self->{'converter_object'} = $converter;
	warn "2. ",Dumper($self);

}

=head2 do_conversion

  Example     : $conversion->do_conversion
  Description : does the database conversion 
  Return type : none
  Exceptions  : none
  Caller      : $Siteroot/utils/vega_schema_conversion

=cut


sub do_conversion {
	my $self= shift;
	$self->{converter_object}->debug( "\n\n*** converting " . $self->{converter_object}->source . " to " . 
					   $self->{converter_object}->target() . " ***");
	$self->{converter_object}->transfer_meta();
	$self->{converter_object}->create_coord_systems();
	$self->{converter_object}->create_seq_regions();
	$self->{converter_object}->create_assembly();
	$self->{converter_object}->create_attribs();
	$self->{converter_object}->set_top_level();
	
	$self->{converter_object}->transfer_dna();
	$self->{converter_object}->transfer_genes();
	$self->{converter_object}->transfer_prediction_transcripts();
	
	if ($self->conv_support->param('do_features')) {
		$self->{converter_object}->transfer_features();
	}
	if ($self->conv_support->param('do_vega_schema_conversion')) {
		$self->{converter_object}->transfer_vega_stable_ids();
	} else {
		$self->{converter_object}->transfer_stable_ids();
	}
	
	$self->{converter_object}->copy_other_tables();
	$self->{converter_object}->copy_repeat_consensus();
	$self->{converter_object}->create_meta_coord();
	if ($self->conv_support->param('do_vega_schema_conversion')) {
		
		$self->{converter_object}->update_clone_info();
		$self->{converter_object}->remove_supercontigs();
		$self->{converter_object}->copy_internal_clone_names();
	}
	print STDERR "*** All finished ***\n";
}

=head2 conv_usage

  Example     : $conversion->conv_usage("message")
  Description : prints usage information and exits
  Return type : none
  Exceptions  : none
  Caller      : $Siteroot/utils/vega_schema_conversion

=cut

sub conv_usage {
	my $self = shift;
	my $msg = shift;

	print STDERR "\nMSG: $msg\n" if($msg);
	
	print STDERR <<EOF;

usage: ./conversion_densities.pl  <options>

options: -pass <password>       the mysql user's password (required)

         -conf <conf_file>      configuration file (required):

                                   fields:
                                      do_vega_conversion (0 or 1)
                                      do_ensembl_conversion (0 or 1)
                                      user (a mysql db user with read/write priveleges)
                                      host (plus port eg ecs3d:3307)
                                      species (common name - choose from dog, mouse, human and zebrafish)
                                      source_db (schema 19 source database)
                                      target_db (schema 20+ trget database)
                                      core_sql (location of ensembl schema creation script eg ensembl/sql/table.sql)
                                      vega_sql (location of creation script for additinoal vega tables eg ensembl/sql/vega_specific_tables.sql)
                                      force (overwrite existing target database? 0 or 1)
                                      verbose (print out debug statements 0 or 1)
                                      do_features (0 or 1 - transfer dna- and protein-align features, for debugging)

         -help                  display this message

EOF
  exit;

}

1;



