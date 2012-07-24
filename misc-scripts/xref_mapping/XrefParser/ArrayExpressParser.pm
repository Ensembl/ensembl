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

  my $species_id_to_names = $self->species_id2name();
  my $names = $species_id_to_names->{$species_id};
  my $contents_lookup = $self->_get_contents($files, $verbose);
  my $active = $self->_is_active($contents_lookup, $names, $verbose);
  $self->_insert_meta($active);
  
  return;
}

sub _get_contents {
  my ($self, $files, $verbose) = @_;
  my %lookup;
  my $fh = $self->get_filehandle($files->[0]);
  while(my $line = <$fh>) {
    chomp $line;
    my ($species, $remainder) = $line =~ /^(\w+)_(.+)$/;
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

sub _insert_meta {
  my ($self, $active) = @_;
  my $sth = $self->dbi->prepare('INSERT INTO meta (meta_key, meta_value) values (?,?)');
  $sth->execute($self->meta_key(), $active);
  $sth->finish();
  return;
}

1;