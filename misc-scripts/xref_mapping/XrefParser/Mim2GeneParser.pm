package XrefParser::Mim2GeneParser;

use strict;
use warnings;
use Carp;
use POSIX qw(strftime);
use File::Basename;

use base qw( XrefParser::BaseParser );

sub run {

 my ($self, $ref_arg) = @_;
  my $general_source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $files        = $ref_arg->{files};
  my $verbose      = $ref_arg->{verbose};

  if((!defined $general_source_id) or (!defined $species_id) or (!defined $files) ){
    croak "Need to pass source_id, species_id and files as pairs";
  }
  $verbose |=0;


  my $file = @{$files}[0];

  my $eg_io = $self->get_filehandle($file);
  if ( !defined $eg_io ) {
    print STDERR "ERROR: Could not open $file\n";
    return 1;    # 1 is an error
  }

  my $entrez_source_id =
      $self->get_source_id_for_source_name(
        'EntrezGene');

  my (%mim_gene)   = %{XrefParser::BaseParser->get_valid_codes("MIM_GENE",$species_id)};
  my (%mim_morbid) = %{XrefParser::BaseParser->get_valid_codes("MIM_MORBID",$species_id)};
  my (%entrez)     = %{XrefParser::BaseParser->get_valid_codes("EntrezGene",$species_id)};
 
  my $dbi = $self->dbi();
  my $add_dependent_xref_sth = $self->dbi->prepare("INSERT INTO dependent_xref  (master_xref_id,dependent_xref_id, linkage_source_id) VALUES (?,?, $entrez_source_id)");

  my $missed_entrez = 0;
  my $missed_omim   = 0;
  my $diff_type     = 0;
  my $count;

  $eg_io->getline(); # do not need header
  while ( $_ = $eg_io->getline() ) {
    $count++;
    chomp;
    my ($omim_id, $entrez_id, $type) = split;

    if(!defined($entrez{$entrez_id})){
      $missed_entrez++;
      next;
    }
    
    if((!defined $mim_gene{$omim_id} ) and (!defined $mim_morbid{$omim_id} ) ){
      $missed_omim++;
      next;
    }

    if($type eq "gene"){
      if(defined($mim_gene{$omim_id})){
	foreach my $ent_id (@{$entrez{$entrez_id}}){
	  foreach my $mim_id (@{$mim_gene{$omim_id}}){
	    $add_dependent_xref_sth->execute($ent_id, $mim_id);
	  }
	}
	# $add_dependent_xref_sth->execute($entrez{$entrez_id}, $mim_gene{$omim_id});
      }
      else{
	$diff_type++;
	foreach my $ent_id (@{$entrez{$entrez_id}}){
	  foreach my $mim_id (@{$mim_morbid{$omim_id}}){
	    $add_dependent_xref_sth->execute($ent_id, $mim_id);
	  }
	}
	# $add_dependent_xref_sth->execute($entrez{$entrez_id}, $mim_morbid{$omim_id});	
      }
    }
    elsif($type eq "phenotype"){
      if(defined($mim_morbid{$omim_id})){
	foreach my $ent_id (@{$entrez{$entrez_id}}){
	  foreach my $mim_id (@{$mim_morbid{$omim_id}}){
	    $add_dependent_xref_sth->execute($ent_id, $mim_id);
	  }
	}
	#	$add_dependent_xref_sth->execute($entrez{$entrez_id}, $mim_morbid{$omim_id});
      }
      else{
	$diff_type++;
	foreach my $ent_id (@{$entrez{$entrez_id}}){
	  foreach my $mim_id (@{$mim_gene{$omim_id}}){
	    $add_dependent_xref_sth->execute($ent_id, $mim_id);
	  }
	}
	#	$add_dependent_xref_sth->execute($entrez{$entrez_id}, $mim_gene{$omim_id});	
      }
    }
    else{
      print "WARNING unknown type $type\n";
      next;
    }

  }
  $add_dependent_xref_sth->finish;

  print $missed_entrez." EntrezGene entries could not be found.\n";
  print $missed_omim." Omim entries could not be found.\n";
  print $diff_type." had different types out of $count Entries.\n";

  return 0;
}

1;
