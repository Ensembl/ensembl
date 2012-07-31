package XrefParser::ArrayExpressParser;

## Parsing format looks like:
# ggallus_gene_ensembl
# hsapiens_gene_ensembl
# mmulatta_gene_ensembl
# mmusculus_gene_ensembl
# osativa_eg_gene
# ptrichocarpa_eg_gene
#

use strict;
use warnings;
use Carp;
use base qw( XrefParser::BaseParser );

sub meta_key {
  my ($self) = @_;
  return "array_express.exported";
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
  $verbose ||=0;

  my %species_id_to_names = $self->species_id2name();
  my $species_id_to_names = \%species_id_to_names;
  my $names = $species_id_to_names->{$species_id};
  my $contents_lookup = $self->_get_contents($files, $verbose);
  my $active = $self->_is_active($contents_lookup, $names, $verbose);

  if ($active && $verbose) {
      print "ArrayExpress xrefs will be created when running xref_mapper.pl/DirectXrefs.pm as gene stable ids are required to create the xrefs\n";
  }
  $self->add_meta_pair($self->meta_key(),$active);
  $self->add_meta_pair('species_id',$species_id);

  return;
}

sub _get_contents {
  my ($self, $files, $verbose) = @_;
  my %lookup;
  my $fh = $self->get_filehandle($files->[0]);
  while(my $line = <$fh>) {
    chomp $line;
    my ($species, $remainder) = $line =~ /^([a-z|A-Z]+)_(.+)$/;
    croak "The line '$line' is not linked to a gene set. This is unexpected." if $remainder !~ /gene/;
    $lookup{$species} = 1;
  }
  close ($fh);
  if($verbose) {
    printf("ArrayExpress is using the species [%s]\n", join(q{, }, keys %lookup));
  }
  return \%lookup;
}

sub _is_active {
  my ($self, $contents_lookup, $names, $verbose) = @_;
  #Loop through the names and aliases first. If we get a hit then great
  my $active = 0;
  foreach my $name (@{$names}) {
    if($contents_lookup->{$name}) {
      printf('Found ArrayExpress has declared the name "%s". This was an alias'."\n", $name) if $verbose;
      $active = 1;
      last;
    }
  }
  
  #Last ditch using the default name and messing around with the name
  if(!$active) {
    my $default_name = $names->[0];
    my ($letter, $remainder) = $default_name =~ /^(\w).+_(.+)$/;
    my $new_name = join(q{}, $letter, $remainder);
    if($contents_lookup->{$new_name}) {
      printf('Found ArrayExpress has declared the name "%s". We have constructed this from the default name'."\n", $new_name) if $verbose;
      $active = 1;
    }
  }
  
  return $active;
}

#this method is called from XrefMapper/DirectXrefs.pm

sub create_xrefs {
    my $self = shift;
    my $verbose = shift;

    my $array_xrefs_meta_key = $self->meta_key();
  
    if ($array_xrefs_meta_key) {
      my $active = $self->get_meta_value($array_xrefs_meta_key);
      if ($active) {
	   #create ArrayExpress direct xrefs
	   my $source_name = 'ArrayExpress';
	   my $source_id = $self->get_source_id_for_source_name($source_name);
	   
	   my $species_id = $self->get_meta_value('species_id');
	   #get gene stable_ids
	   my $gene_id_sth = $self->dbi()->prepare("select stable_id from gene_stable_id order by stable_id");
	   $gene_id_sth->execute();
	   my $gene_stable_id;
	   $gene_id_sth->bind_columns(\$gene_stable_id);
	   my $xref_count = 0;
	   while ($gene_id_sth->fetch()) {
	       
	       my $xref_id = $self->add_xref({ acc        => $gene_stable_id,
					 label      => $gene_stable_id,
					 source_id  => $source_id,
					 species_id => $species_id,
					 info_type => "DIRECT"} );
	
	       $self->add_direct_xref( $xref_id, $gene_stable_id, 'gene', '');
	       if ($xref_id) {
		   $xref_count++;
	       }
	   }
	   $gene_id_sth->finish();
	   
	   if ($xref_count > 0) {
	       print "Loaded $xref_count $source_name DIRECT xrefs\n" if $verbose;	       
	   } else {

	       print "Warning: 0 $source_name DIRECT xrefs loaded even though $source_name is active for the species.\n";
	   }
      }

    }
}


1;
