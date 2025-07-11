=head1 LICENSE

See the NOTICE file distributed with this work for additional information
regarding copyright ownership.

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

package XrefParser::ZFINParser;

use strict;
use warnings;
use Carp;
use File::Basename; # provides dirname
use File::Spec::Functions;
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

  my $file = @{$files}[0];
  my $dir = dirname($file);

  # Get the ZFIN source ids
  my $direct_src_id = $self->get_source_id_for_source_name('ZFIN_ID', 'direct', $dbi);
  my $dependent_src_id = $self->get_source_id_for_source_name('ZFIN_ID', 'uniprot/refseq', $dbi);
  my $description_src_id = $self->get_source_id_for_source_name('ZFIN_ID', 'description_only', $dbi);

  # Get the ZFIN descriptions
  my %description;

  my $sth = $dbi->prepare("select accession, description from xref where source_id=?");
  $sth->execute($description_src_id);
  my ($acc, $desc);
  my $zfin_loaded_count = 0;
  $sth->bind_columns(\$acc, \$desc);
  while (my @row = $sth->fetchrow_array()) {
    $description{$acc} = $desc if(defined($desc));
    $zfin_loaded_count++;
  }
  $sth->finish;

  # Get the Uniprot and RefSeq accessions
  my (%swiss) = %{$self->get_valid_codes("uniprot/swissprot",$species_id, $dbi)};
  my (%refseq) = %{$self->get_valid_codes("refseq",$species_id, $dbi)};

  # Process ZFIN to ensEMBL mappings
  my %zfin;
  my $zfin_io = $self->get_filehandle(catfile($dir, 'ensembl_1_to_1.txt'));
  if (!defined($zfin_io)) {
    croak "ERROR: Could not open " . catfile($dir, 'ensembl_1_to_1.txt') . "\n";
  }

  my $zfin_csv = Text::CSV->new({
    sep_char       => "\t",
    empty_is_undef => 1,
    strict         => 1,
  }) or croak "Could not use zfin file: " . Text::CSV->error_diag();

  $zfin_csv->column_names(['zfin', 'so', 'label', 'ensembl_id']);

  while (my $zfin_line = $zfin_csv->getline_hr($zfin_io)) {
    my ($zfin_acc, $so, $label, $ensembl_id) = @{$zfin_line}{qw(zfin so label ensembl_id)};

    $self->add_to_direct_xrefs({
      stable_id  => $ensembl_id,
      type       => 'gene',
      acc        => $zfin_acc,
      label      => $label,
      desc       => $description{$zfin_acc},
      dbi        => $dbi,
      source_id  => $direct_src_id,
      species_id => $species_id
    });

    $zfin{$zfin_acc} = 1;
  }

  $zfin_io->close();

  my $spcount =0;
  my $rscount =0;
  my $mismatch=0;

  # Process ZFIN to Uniprot mappings
  my $swissprot_io = $self->get_filehandle( catfile( $dir, 'uniprot.txt' ) );
  if ( !defined $swissprot_io ) {
    croak "ERROR: Could not open " . catfile( $dir, 'uniprot.txt' ). "\n" ;
  }

  my $swissprot_csv = Text::CSV->new({
    sep_char       => "\t",
    empty_is_undef => 1,
    strict         => 1,
  }) or croak "Could not use swissprot file $file: " . Text::CSV->error_diag();

  $swissprot_csv->column_names([ 'zfin', 'so', 'label', 'acc' ]);

  #swissprot file format (in uniprot.txt)
  #ZDB-GENE-000112-47      SO:0000704      ppardb  Q90Z66
  #ZDB-GENE-000125-12      SO:0000704      igfbp2a Q9PTH3
  #ZDB-GENE-000125-4       SO:0000704      dlc     B3DFM3

  while ( my $swissprot_line = $swissprot_csv->getline_hr( $swissprot_io ) ) {
    my ($zfin_acc, $so, $label, $acc) = @{$swissprot_line}{qw(zfin so label acc)};

    if(defined($swiss{$acc}) && !defined($zfin{$zfin_acc})){
      foreach my $xref_id (@{$swiss{$acc}}){
        $self->add_dependent_xref({
          master_xref_id => $xref_id,
          acc            => $zfin_acc,
          label          => $label,
          desc           => $description{$zfin_acc},
          source_id      => $dependent_src_id,
          dbi            => $dbi,
          species_id     => $species_id
        });
        $spcount++;
      }
    } else {
      $mismatch++;
    }
  }

  $swissprot_io->close();

  # Process ZFIN to RefSeq mappings
  my $refseq_io = $self->get_filehandle( catfile( $dir, 'refseq.txt' ) );
  if ( !defined $refseq_io ) {
    croak "ERROR: Could not open " . catfile( $dir, 'refseq.txt' ),"\n" ;
  }

  my $refseq_csv = Text::CSV->new({
    sep_char       => "\t",
    empty_is_undef => 1,
    strict         => 1,
  }) or croak "could not use refseq file $file: " . Text::CSV->error_diag();

  $refseq_csv->column_names([ 'zfin', 'so', 'label', 'acc' ]);

  #ZDB-GENE-000125-12      SO:0000704      igfbp2a NP_571533
  #ZDB-GENE-000125-4       SO:0000704      dlc     NM_130944
  #ZDB-GENE-000125-4       SO:0000704      dlc     NP_571019
  #ZDB-GENE-000128-11      SO:0000704      dbx1b   NM_131178

  while ( my $refseq_line = $refseq_csv->getline_hr( $refseq_io ) ) {
    my ($zfin_acc, $so, $label, $acc) = @{$refseq_line}{qw(zfin so label acc)};
    # Ignore mappings to predicted RefSeq
    if ($acc =~ /^XP_/ || $acc =~ /^XM_/ || $acc =~ /^XR_/) { next; }

    if(defined($refseq{$acc}) && !defined($zfin{$zfin_acc})){
      foreach my $xref_id (@{$refseq{$acc}}){
        $self->add_dependent_xref({
          master_xref_id => $xref_id,
          acc            => $zfin_acc,
          label          => $label,
          desc           => $description{$zfin_acc},
          source_id      => $dependent_src_id,
          dbi            => $dbi,
          species_id     => $species_id
        });
        $rscount++;
      }
    } else {
      $mismatch++;
    }
  }

  $refseq_io->close();

  # Get the added ZFINs again (with deps)
  (%zfin) = %{$self->get_valid_codes("zfin", $species_id, $dbi)};

  # Process the synonyms
  my $aliases_io = $self->get_filehandle( catfile( $dir, 'aliases.txt' ) );
  if ( !defined $aliases_io ) {
    croak "ERROR: Could not open " . catfile( $dir, 'aliases.txt' ), "\n" ;
  }

  my $aliases_csv = Text::CSV->new({
    sep_char       => '\t',
    empty_is_undef => 1,
    strict         => 1,
  }) or croak "could not use zfin file $file: " . Text::CSV->error_diag();

  $aliases_csv->column_names([ 'acc', 'cur_name', 'cur_symbol', 'syn', 'so' ]);

  #DB-ALT-000717-2        zc1Tg   zc1Tg   zc1     SO:0001218
  #ZDB-ALT-000717-4        zc3Tg   zc3Tg   Tg(NBT:MAPT-GFP)        SO:0001218

  my $syncount = 0;

  $sth = $dbi->prepare('SELECT source_id from source where name like "ZFIN_ID"');
  $sth->execute;
  my $s1;
  $sth->bind_columns(\$s1);
  my $sources;
  while($sth->fetch()){
    push @$sources, $s1;
  }
  $sth->finish;

  while ( my $aliases_line = $aliases_csv->getline_hr( $aliases_io ) ) {
    my ($acc, $syn) = @{$aliases_line}{qw(acc syn)};
    if(defined($zfin{$acc})){
      $self->add_to_syn_for_mult_sources($acc, $sources, $syn, $species_id, $dbi);
      $syncount++;
    }
  }

  $aliases_io->close();

  if($verbose){
    print "\t$spcount xrefs from UniProt and\n";
    print "\t$rscount xrefs from RefSeq succesfully loaded\n";
    print "\t$syncount synonyms loaded\n";
    print "\t$mismatch xrefs ignored\n";
  }
  return 0;
}

1;
