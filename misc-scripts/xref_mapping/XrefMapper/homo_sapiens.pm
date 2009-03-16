package XrefMapper::homo_sapiens;

use  XrefMapper::BasicMapper;
use  XrefMapper::SubmitMapper;
use strict;
use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };

sub get_set_lists {

  return [["ExonerateGappedBest1", ["homo_sapiens","*"]]];

}

sub gene_description_filter_regexps {

  return ('^BA\S+\s+\(NOVEL PROTEIN\)\.?',
	  '^DJ\S+\s+\(NOVEL PROTEIN\)\.?',
	  '^LOC\d+\s*(PROTEIN)?\.?',
	  '^ORF.*',
	  '^PROTEIN C\d+ORF\d+\.*',
	  '\(CLONE \S+\)\s+',
	  '^BC\d+\_\d+\.?',
	  '^CGI\-\d+ PROTEIN\.?\;?',
	  '[0-9A-Z]{10}RIK PROTEIN[ \.]',
	  'R\d{5}_\d[ \.,].*',
	  'PROTEIN KIAA\d+[ \.].*',
	  'RIKEN CDNA [0-9A-Z]{10}[ \.]',
	  '^\(*HYPOTHETICAL\s+.*',
	  '^UNKNOWN\s+.*',
	  '^DKFZP[A-Z0-9]+\s+PROTEIN[\.;]?.*',
	  '^CHROMOSOME\s+\d+\s+OPEN\s+READING\s+FRAME\s+\d+\.?.*',
	  '^FKSG\d+\.?.*',
	  '^HSPC\d+\s+PROTEIN\.?.*',
	  '^KIAA\d+\s+PROTEIN\.?.*',
	  '^KIAA\d+\s+GENE\s+PRODUCT\.?.*',
	  '^HSPC\d+.*',
	  '^PRO\d+\s+PROTEIN\.?.*',
	  '^PRO\d+\.?.*',
	  '^FLJ\d+\s+PROTEIN.*',
	  '^PRED\d+\s+PROTEIN.*',
	  '^WUGSC:.*\s+PROTEIN\.?.*',
	  '^SIMILAR TO GENE.*',
	  '^SIMILAR TO PUTATIVE[ \.]',
	  '^SIMILAR TO HYPOTHETICAL.*',
	  '^SIMILAR TO (KIAA|LOC).*',
	  '^SIMILAR TO\s+$',
          '^WUGSC:H_.*',
          '^\s*\(?PROTEIN\)?\.?\s*$',
	  '^\s*\(?FRAGMENT\)?\.?\s*$',
          '^\s*\(?GENE\)?\.?\s*$',
	  '^\s*\(\s*\)\s*$',
          '^\s*\(\d*\)\s*[ \.]$');

}


sub get_official_name{
   return "HGNC";
}

sub get_canonical_name{
   return "HGNC";
}

sub species_specific_cleanup{
  my $self = shift;
  my $dbname = $self->get_canonical_name;

  print "Removing all $dbname from object_xref not on a Gene\n";
  my $remove_old_ones = (<<JSQL);
delete ox 
  from object_xref ox, xref x, external_db e
    where e.db_name like "$dbname" and 
          ox.ensembl_object_type != "Gene" and
          ox.xref_id = x.xref_id and
	  x.external_db_id = e.external_db_id;
JSQL

  #
  # First Delete all the hgnc object_xrefs not on a gene. (i.e these are copys).
  #

  my $sth = $self->core->dbc->prepare($remove_old_ones);

  $sth->execute() || die "Could not execute: \n$remove_old_ones \n";

  $sth->finish;

}


# For human we want to make a copy of the HGNC references on the genes and put them on 
# the "canonical" transcripts

#sub species_specific_pre_attributes_set{
#  my $self  = shift;
#  $self->official_naming();
#}



1;
