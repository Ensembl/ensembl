package XrefParser::CoreXrefParser;

use strict;
use warnings;
use Carp;
use DBI;

use base qw( XrefParser::BaseParser );
use Bio::EnsEMBL::Registry;

sub run_script {

  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $file         = $ref_arg->{file};
  my $verbose      = $ref_arg->{verbose};
  my $mapper       = $ref_arg->{mapper};

  if((!defined $source_id) or (!defined $species_id) or (!defined $file) ){
    croak "Need to pass source_id, species_id and file as pairs";
  }
  if(!defined $mapper) {
    croak "Need to connect to a core database; please use a mapper based configuration";
  }
  $verbose |=0;

  my $opts = $self->parse_opts_from_file($file);
  my $biotype = $opts->{biotype};
  my $object_type = lc($opts->{object_type}) || 'gene';
  my $copy_description_from_object = $opts->{copy_description_from_object};

  my $external_db_name = $self->get_source_name_for_source_id($source_id);

  #copy object xrefs from core
  my %valid_object_types = (
    gene => 'Gene', transcript => 'Transcript', translation => 'Translation',
  );

  if (!exists($valid_object_types{$object_type}) ) {
      die("Unsupported object type value. Supported values: ", join(',', keys %valid_object_types) );
  }

  if ($biotype && $object_type ne 'gene' && $object_type ne 'transcript') {
      die("Incorrect parser argument values: expecting gene or transcript object type when biotype provided.\n");
  }

  my $dba = $mapper->core()->dba();
  my $object_adaptor = Bio::EnsEMBL::Registry->get_adaptor($dba->species(), 'core', $object_type);

  my @objects;

  if ($biotype) {
      @objects = @{$object_adaptor->fetch_all_by_biotype($biotype)};
  } else {
      @objects = @{$object_adaptor->fetch_all()};
  }

  my %added_xref;
  my $direct_count;

  foreach my $object (@objects) {

      my @xrefs = @{$object->get_all_DBEntries($external_db_name)};

      foreach my $xref (@xrefs) {

	  my $xref_id;

	  if (!exists($added_xref{$xref->primary_id()})) {

	      my $description = $xref->description();

	      if ($copy_description_from_object && !$description) {

		  if ($object->description()) {
                      #populate xref description with object description stripping the [Source: .. part
		      ($description) = $object->description() =~ /([^\[]+)/;
		      #trim trailing spaces
		      $description =~ s/\s+$//; 
		  }
	      }
	      
	      $xref_id = $self->add_xref({ acc        => $xref->primary_id(),
				      version    => $xref->version(),
				      label      => $xref->display_id(),
				      desc       => $description,
				      source_id  => $source_id,
				      species_id => $species_id,
				      info_type  => "DIRECT"} );


	      $added_xref{$xref->primary_id()} = $xref_id;
	  } 

	  if (!$xref_id) {
	      $xref_id = $added_xref{$xref->primary_id()};
	  }
	  
	  $self->add_direct_xref($xref_id, $object->stable_id(), $valid_object_types{$object_type}, "");
	  $direct_count++;
      }
  }

  my $xref_count = scalar(keys %added_xref);

  print "Added $xref_count $external_db_name xrefs and $direct_count $object_type direct xrefs\n" if($verbose);
  if ( !$xref_count ) {
      return 1;    # 1 error
  }

  return 0; # successfull
 

}

1;
