package Bio::EnsEMBL::IdMapping::Archiver;

=head1 NAME


=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS


=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http:#www.ensembl.org/info/about/code_licence.html for details.

=head1 AUTHOR

Patrick Meidl <meidl@ebi.ac.uk>, Ensembl core API team

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<ensembl-dev@ebi.ac.uk>

=cut


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


sub dump_gene {
  my ($self, $s_gene, $t_gene, $ga_fh, $pa_fh) = @_;

  # private method, so no argument check done for performance reasons
  
  # deal with ncRNA differently
  # hope this simple biotype regex is accurate enough...
  my $is_ncRNA = 0;
  $is_ncRNA = 1 if ($s_gene->biotype =~ /RNA/);

  # loop over source transcripts
  foreach my $s_tr (@{ $s_gene->get_all_Transcripts }) {
    my $s_tl = $s_tr->translation;

    # we got a coding transcript
    if ($s_tl) {
      
      # do a full dump of this gene if no target gene exists
      if (! $t_gene) {
        $self->dump_tuple($s_gene, $s_tr, $s_tl, $ga_fh, $pa_fh);

      # otherwise, only dump if translation of this transcript changed
      } else {

        my $changed_flag = 1;
      
        foreach my $t_tr (@{ $t_gene->get_all_Transcripts }) {
          my $t_tl = $t_tr->translation;
          next unless ($t_tl);

          if (($s_tr->stable_id eq $t_tr->stable_id) and
              ($s_tl->stable_id eq $t_tl->stable_id) and
              ($s_tl->seq eq $t_tl->seq)) {
            $changed_flag = 0;
          }
        }

        if ($changed_flag) {
          $self->dump_tuple($s_gene, $s_tr, $s_tl, $ga_fh, $pa_fh);
        }
      }

    # now deal with ncRNAs (they don't translate but we still want to archive
    # them)
    } elsif ($is_ncRNA) {
    
      if (! $t_gene) {
    
        $self->dump_nc_row($s_gene, $s_tr, $ga_fh);
    
      } else {
        
        my $changed_flag = 1;
      
        foreach my $t_tr (@{ $t_gene->get_all_Transcripts }) {
          $changed_flag = 0 if ($s_tr->stable_id eq $t_tr->stable_id);
        }

        if ($changed_flag) {
          $self->dump_nc_row($s_gene, $s_tr, $ga_fh);
        }
      
      }
    }
  }
}


sub dump_tuple {
  my ($self, $gene, $tr, $tl, $ga_fh, $pa_fh) = @_;
  
  # private method, so no argument check done for performance reasons

  # gene archive
  print $ga_fh join("\t",
                    $gene->stable_id,
                    $gene->version,
                    $tr->stable_id,
                    $tr->version,
                    $tl->stable_id,
                    $tl->version,
                    $pa_id,
                    $self->mapping_session_id
  );
  print $ga_fh "\n";

  # peptide archive
  my $pep_seq = $tl->seq;
  print $pa_fh join("\t", $pa_id, md5_hex($pep_seq), $pep_seq);
  print $pa_fh "\n";

  # increment peptide_archive_id
  $pa_id++;
}


sub dump_nc_row {
  my ($self, $gene, $tr, $ga_fh) = @_;
  
  # private method, so no argument check done for performance reasons

  # gene archive
  print $ga_fh join("\t",
                    $gene->stable_id,
                    $gene->version,
                    $tr->stable_id,
                    $tr->version,
                    '\N',
                    '\N',
                    '\N',
                    $self->mapping_session_id
  );
  print $ga_fh "\n";
}


sub mapping_session_id {
  my $self = shift;
  $self->{'_mapping_session_id'} = shift if (@_);
  return $self->{'_mapping_session_id'};
}


1;

