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

package XrefParser::RGDParser;

use strict;
use warnings;
use Carp;
use POSIX qw(strftime);
use File::Basename;
use Text::CSV;

use parent qw( XrefParser::BaseParser );

sub run {

  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $files        = $ref_arg->{files};
  my $verbose      = $ref_arg->{verbose};
  my $dbi          = $ref_arg->{dbi};
  $dbi = $self->dbi unless defined $dbi;

  if((!defined $source_id) or (!defined $species_id) or (!defined $files) ){
    croak "Need to pass source_id, species_id and files as pairs";
  }
  $verbose |=0;

  my $source_sql = "select source_id from source where name = 'RGD' and priority_description = 'direct_xref'";
  my $sth = $dbi->prepare($source_sql);
  $sth->execute();
  my ($direct_source_id);
  $sth->bind_columns(\$direct_source_id);
  $sth->fetch();
  $sth->finish();

  my $file = @{$files}[0];

  # Used to assign dbIDs for when RGD Xrefs are dependent on RefSeq xrefs
  my (%preloaded_refseq) = %{$self->get_valid_codes("refseq",$species_id, $dbi)};
  
  my $rgd_io = $self->get_filehandle($file);

  if ( !defined $rgd_io ) {
    carp "Could not open $file when trying to parse RGD";
  }
  my $csv = Text::CSV->new({ sep => "\t", blank_is_undef => 1, strict => 0, auto_diag => 1, binary => 1, allow_loose_quotes => 1})
    or carp "Cannot use CSV: ".Text::CSV->error_diag ();
  # WARNING - Text::CSV does not like the GENES-RAT.txt file. It is improperly formatted and contains a non-ASCII character
  # Make sure binary is turned on or it silently fails and you get 1/3rd of the records.
  # strict is turned off to prevent failure on a blank line at the end

  my $line = '#';
  until (substr($line,0,1) ne '#') {
    $line = $rgd_io->getline;
  }
  $csv->parse($line);
  $csv->column_names($csv->fields);
  # Columns we want
  #  GENE_RGD_ID => 0,
  #  SYMBOL => 1,
  #  NAME => 2,
  #  GENBANK_NUCLEOTIDE => 23,
  #  OLD_SYMBOL => 29,
  #  ENSEMBL_ID => 37

  my $count= 0;
  my $ensembl_count = 0;
  my $mismatch = 0;
  my $syn_count = 0;
  my $cols; # Digested columns from CSV
  while ( $cols = $csv->getline_hr($rgd_io) ) {
    next if $csv->is_missing (0) || (exists $cols->{GENE_RGD_ID} && $cols->{GENE_RGD_ID} eq '' || !defined $cols->{GENE_RGD_ID});

    my @nucs;
    @nucs = split(/\;/,$cols->{GENBANK_NUCLEOTIDE}) if defined $cols->{GENBANK_NUCLEOTIDE};
    my $done = 0;
    # @nucs are sorted in the file in alphabetical order. Filter them to down
    # to a higher quality subset, then add dependent Xrefs where possible
    foreach my $nuc ($self->sort_refseq_accessions(@nucs)){
      if(!$done && exists $preloaded_refseq{$nuc}){
        foreach my $xref (@{$preloaded_refseq{$nuc}}){
          my $xref_id = $self->add_dependent_xref({ 
                      master_xref_id => $xref,
                      acc            => $cols->{GENE_RGD_ID},
                      label          => $cols->{SYMBOL},
                      desc           => $cols->{NAME},
                      source_id      => $source_id,
                      dbi            => $dbi,
                      species_id     => $species_id} );
          $count++;
          $syn_count += $self->add_synonyms($cols->{OLD_SYMBOL},$xref_id,$dbi);
          $done = 1;
        }
      }
    }

    if (defined $cols->{ENSEMBL_ID}) {
      my @ensembl_ids = split(/\;/, $cols->{ENSEMBL_ID});
      foreach my $id (@ensembl_ids) {
        $ensembl_count++;
        $self->add_to_direct_xrefs({ stable_id => $id,
                                     type => 'gene',
                                     acc => $cols->{GENE_RGD_ID},
                                     label => $cols->{SYMBOL},
                                     desc => $cols->{NAME},
                                     dbi  => $dbi,
                                     source_id => $direct_source_id,
                                     species_id => $species_id} );
        my $xref_id = $self->get_xref($cols->{GENE_RGD_ID}, $direct_source_id, $species_id, $dbi);
        $syn_count += $self->add_synonyms($cols->{OLD_SYMBOL},$xref_id,$dbi);
        $done = 1;
      }
    }
    if(!$done){
      $self->add_xref({ 
        acc        => $cols->{GENE_RGD_ID},
        label      => $cols->{SYMBOL},
        desc       => $cols->{NAME},
        source_id  => $source_id,
        species_id => $species_id,
        dbi        => $dbi,
        info_type  => "MISC"} );
      $mismatch++;
    }

  }
  if (! $csv->eof) {
    confess 'Failed to finish parsing RGD file: '.$csv->error_diag();
  }
  $rgd_io->close();

  if($verbose){
    print "\t$count xrefs succesfully loaded and dependent on refseq\n";
    print "\t$mismatch xrefs added but with NO dependencies\n";
    print "\t$ensembl_count direct xrefs successfully loaded\n";
    print "\tTried to add $syn_count synonyms, including duplicates\n";
  }
  return 0;
}

# Predefined importance levels for the most valued RefSeq accession types
my %refseq_priorities = (
  NM => 1,
  NP => 1,
  NR => 1,
  XM => 2,
  XP => 2,
  XR => 2,
);

# Filter out any accessions which are not in the "normal" set of genomic features
# The column in question contains EMBL accessions as well as other things, and we don't 
# have the ability to make Xrefs to those. It wouldn't be useful to do so in any case
sub sort_refseq_accessions {
  my ($self,@accessions) = @_;
  @accessions = sort { $refseq_priorities{substr $a, 0,2} <=> $refseq_priorities{substr $b, 0,2} } 
                grep { exists $refseq_priorities{substr $_, 0,2} } @accessions;
  return @accessions;
}


sub add_synonyms {
  my ($self,$synonym_string,$xref_id,$dbi) = @_;
  my $syn_count = 0;
  return $syn_count if (!defined $synonym_string || !defined $xref_id || !defined $dbi);
  my $add_syn_sth;

  if (! $self->{_syn_dbh_cache}) {
    my $sql = "insert ignore into synonym (xref_id, synonym) values (?, ?)";
    $add_syn_sth = $dbi->prepare($sql);
    $self->{_syn_dbh_cache} = $add_syn_sth;
  }
  $add_syn_sth = $self->{_syn_dbh_cache};

  my @syns;
  @syns = split(/\;/,$synonym_string);
  foreach my $syn(@syns){
    $add_syn_sth->execute($xref_id, $syn);
    $syn_count++;
  }
  return $syn_count;
}

1;
