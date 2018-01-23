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


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::IdMapping::Archiver - create gene_archive and peptide_archive

=head1 SYNOPSIS

  my $archiver = Bio::EnsEMBL::IdMapping::Archiver->new(
    -LOGGER => $logger,
    -CONF   => $conf,
    -CACHE  => $cache
  );

  # create gene and peptide archive
  $archiver->create_archive($mapping_session_id);

  # dump existing archive tables to file
  my $num_entries =
    $archiver->dump_table_to_file( 'source', 'gene_archive',
    'gene_archive_existing.txt', 1 );

=head1 DESCRIPTION

This module creates the gene_archive and peptide_archive
tables. Data is written to a file as tab-delimited text for
loading into a MySQL database (this can be done manually, or using
StableIdmapper->upload_file_into_table()).

An archive entry for a given source gene is created if no target
gene exists, or if any of its transcripts or their translations
changed. Non-coding transcripts only have an entry in gene_archive (i.e.
without a corresponding peptide_archive entry).

=head1 METHODS

  create_archive
  dump_gene
  dump_tuple
  dump_nc_row
  mapping_session_id

=cut


package Bio::EnsEMBL::IdMapping::Archiver;

use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::IdMapping::BaseObject;
our @ISA = qw(Bio::EnsEMBL::IdMapping::BaseObject);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::ScriptUtils qw(path_append);
use Digest::MD5 qw(md5_hex);


# instance variables
my $pa_id;


=head2 create_archive

  Arg[1]      : Int $mapping_session_id - the mapping_session_id for this run
  Example     : $archiver->create_archive($stable_id_mapper->mapping_session_id);
  Description : Creates the gene_archive and peptide_archive tables and writes
                the data to a tab-delimited file. The decision as to what to
                archive is deferred to dump_gene(), see documentation there for
                details.
  Return type : none
  Exceptions  : Thrown on missing argument.
  Caller      : id_mapping.pl
  Status      : At Risk
              : under development

=cut

sub create_archive {
  my $self = shift;
  my $mapping_session_id = shift;

  # argument check
  unless ($mapping_session_id) {
    $self->logger->warning("No mapping_session_id set.");
  }

  $self->mapping_session_id($mapping_session_id);

  # get filehandles to write gene and peptide archive
  my $ga_fh = $self->get_filehandle('gene_archive_new.txt', 'tables');
  my $pa_fh = $self->get_filehandle('peptide_archive_new.txt', 'tables');

  # get the currently highest peptide_archive_id from the source db
  my $s_dba = $self->cache->get_DBAdaptor('source');
  my $s_dbh = $s_dba->dbc->db_handle;
  my $sql = qq(SELECT MAX(peptide_archive_id) FROM peptide_archive);
  $pa_id = $self->fetch_value_from_db($s_dbh, $sql);

  unless ($pa_id) {
    $self->logger->warning("No max(peptide_archive_id) found in db.\n", 1);
    $self->logger->info("That's ok if this is the first stable ID mapping for this species.\n", 1);
  }
  
  $pa_id++;
  $self->logger->debug("Starting with peptide_archive_id $pa_id.\n");

  # lookup hash of target gene stable IDs
  my %target_genes = map { $_->stable_id => $_ }
    values %{ $self->cache->get_by_name("genes_by_id", 'target') };

  # loop over source genes and dump to archive (dump_gene() will decide whether
  # to do full or partial dump)
  foreach my $source_gene (values %{ $self->cache->get_by_name("genes_by_id",
    'source') }) {

    $self->dump_gene($source_gene, $target_genes{$source_gene->stable_id}, 
      $ga_fh, $pa_fh);
  }
  
  close($ga_fh);
  close($pa_fh);
}


=head2 dump_gene

  Arg[1]      : Bio::EnsEMBL::IdMapping::TinyGene $s_gene - source gene
  Arg[2]      : Bio::EnsEMBL::IdMapping::TinyGene $t_gene - target gene
  Arg[3]      : Filehandle $ga_fh - filehandle for writing gene_archive data
  Arg[4]      : Filehandle $pa_fh - filehandle for writing peptide_archive data
  Example     : my $target_gene = $gene_mappings{$source_gene->stable_id};
                $archiver->dump_gene($source_gene, $target_gene, $ga_fh, $pa_fh);
  Description : Given a source gene, it will write a gene_achive and
                peptide_achive entry for it if no target gene exists, or if any
                of its transcripts or their translation changed. 
  Return type : none
  Exceptions  : none
  Caller      : create_archive()
  Status      : At Risk
              : under development

=cut

sub dump_gene {
  my ($self, $s_gene, $t_gene, $ga_fh, $pa_fh) = @_;

  # private method, so no argument check done for performance reasons
  
  # loop over source transcripts
  foreach my $s_tr (@{ $s_gene->get_all_Transcripts }) {

    my $changed_flag = 1;
    my $t_tl;
    my $s_tl = $s_tr->translation;

    if ( $t_gene ) {

      foreach my $t_tr (@{ $t_gene->get_all_Transcripts }) {

        $t_tl = $t_tr->translation;

        # If there is a translation, there should also be a target translation
        # and both transcript and translation should have same stable id and version
        if ( $s_tl ) {

          if ( $t_tl ) {

            if ( $s_tl->stable_id eq $t_tl->stable_id
                 and $s_tl->version == $t_tl->version
                 and $s_tr->stable_id eq $t_tr->stable_id
                 and $s_tr->version == $t_tr->version ) {

              $changed_flag = 0;

            }
          }
        } else {

        # nothing changed if source transcript has a target transcript with same stable id and version
          if ( $s_tr->stable_id eq $t_tr->stable_id
               and $s_tr->version == $t_tr->version ) {

            $changed_flag = 0;

          }
        }
      }
    }

    if ($changed_flag) {
      $self->dump_tuple($s_gene, $s_tr, $s_tl, $ga_fh, $pa_fh);
    }
  }
}


=head2 dump_tuple

  Arg[1]      : Bio::EnsEMBL::IdMapping::TinyGene $gene - gene to archive
  Arg[2]      : Bio::EnsEMBL::IdMapping::TinyTrancript $tr - its transcript
  Arg[3]      : Bio::EnsEMBL::IdMapping::TinyTranslation $tl - its translation
  Arg[4]      : Filehandle $ga_fh - filehandle for writing gene_archive data
  Arg[5]      : Filehandle $pa_fh - filehandle for writing peptide_archive data
  Example     : $archive->dump_tuple($s_gene, $s_tr, $s_tl, $ga_fh, $pa_fh);
  Description : Writes entry lines for gene_archive and peptide_archive.
  Return type : none
  Exceptions  : none
  Caller      : dump_gene()
  Status      : At Risk
              : under development

=cut

sub dump_tuple {
  my ($self, $gene, $tr, $tl, $ga_fh, $pa_fh) = @_;
  
  # private method, so no argument check done for performance reasons

  # gene archive
  print $ga_fh join("\t",
                    $gene->stable_id,
                    $gene->version,
                    $tr->stable_id,
                    $tr->version
  );

  print $ga_fh "\t";

  if ( $tl ) {
    print $ga_fh join("\t",
                    $tl->stable_id,
                    $tl->version,
                    $pa_id,
                    $self->mapping_session_id
    );

    # peptide archive
    my $pep_seq = $tl->seq;
    print $pa_fh join("\t", $pa_id, md5_hex($pep_seq), $pep_seq);
    print $pa_fh "\n";
    # increment peptide_archive_id
    $pa_id++;
  } else {
    print $ga_fh join ("\t",
                       '\N',
                       '\N',
                       '\N',
                       $self->mapping_session_id
    );
  }
  print $ga_fh "\n";
}


=head2 mapping_session_id

  Arg[1]      : (optional) Int - mapping_session_id to set
  Example     : my $msi = $archiver->mapping_session_id;
  Description : Getter/setter for mapping_session_id.
  Return type : Int
  Exceptions  : none
  Caller      : create_archive()
  Status      : At Risk
              : under development

=cut

sub mapping_session_id {
  my $self = shift;
  $self->{'_mapping_session_id'} = shift if (@_);
  return $self->{'_mapping_session_id'};
}


1;

