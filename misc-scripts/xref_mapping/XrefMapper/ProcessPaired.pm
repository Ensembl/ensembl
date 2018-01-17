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

package XrefMapper::ProcessPaired;
use strict;

use vars '@ISA';
@ISA = qw{ XrefMapper::BasicMapper };

use warnings;
use XrefMapper::BasicMapper;

use Cwd;
use DBI;
use File::Basename;
use IPC::Open3;

sub new {
  my($class, $mapper) = @_;

  my $self ={};
  bless $self,$class;
  $self->xref($mapper->xref);
  $self->verbose($mapper->verbose);
  return $self;
}


sub process{
  my ($self) = @_;

  #get all 'DUMP_OUT' transcript RefSeq object xrefs (ccds priority or refseq sequence matched)
  #foreach transcript refseq find its protein pair, check if it's matched to the corresponding translation 
  #if it's not add 'INFERRED_PAIR' object xref

  #set ox_status to 'MULTI_DELETE' for all translation RefSeq object xrefs if a better object_xref exists (better object_xref is one whose corresponding transcript is linked to the paired RefSeq_mRNA)
  

  print "Process Pairs\n" if($self->verbose);
  my $dbi = $self->xref->dbc;
  my $object_xref_id;

  #this query gives us transcript RefSeq_mRNA% object xrefs, and the paired RefSeq_peptide% accession as well as the translation id for the transcript 
  my $transcr_obj_xrefs_sth = $dbi->prepare("select gtt.translation_id, p.source_id, p.accession1, ix.query_identity, ix.target_identity from object_xref ox join xref x on (ox.xref_id = x.xref_id and ox.ox_status = 'DUMP_OUT') join source s on (x.source_id = s.source_id and s.name like 'RefSeq\_mRNA%') join pairs p on (x.accession = p.accession2) join gene_transcript_translation gtt on (gtt.transcript_id = ox.ensembl_id) join identity_xref ix using(object_xref_id)");

  #this query is used to check if and object_xref exists for the related translation and paired RefSeq_peptide% with a status of 'DUMP_OUT'
  my $ox_translation_sth =  $dbi->prepare("select ox.object_xref_id, ox.xref_id from object_xref ox join xref x using(xref_id) where ox.ox_status in ('DUMP_OUT', 'FAILED_PRIORITY') and ox.ensembl_object_type = 'Translation' and ox.ensembl_id = ? and x.source_id = ? and x.accession = ?");
 
  my $ox_insert_sth = $dbi->prepare("insert into object_xref (xref_id, ensembl_id, ensembl_object_type, linkage_type, ox_status) values(?, ?, ?, 'INFERRED_PAIR', 'DUMP_OUT')");
  my $get_object_xref_id_sth = $dbi->prepare("select object_xref_id from object_xref where xref_id = ? and ensembl_id = ? and ensembl_object_type = ? and linkage_type = 'INFERRED_PAIR' and ox_status = 'DUMP_OUT'");

  my $xref_sth =  $dbi->prepare("select xref_id from xref where accession = ? and source_id = ?");

  my $xref_update_sth =  $dbi->prepare("update xref set info_type = 'INFERRED_PAIR' where xref_id = ?");
  my $identity_update_sth = $dbi->prepare("insert into identity_xref (object_xref_id, query_identity, target_identity) values(?, ?, ?)");

  my $transl_object_xrefs_sth = $dbi->prepare("select ox.object_xref_id, ox.ensembl_id, x.accession, gtt.transcript_id from gene_transcript_translation gtt join object_xref ox on (gtt.translation_id = ox.ensembl_id and ox.ensembl_object_type = 'Translation') join xref x on (ox.xref_id = x.xref_id and ox.ox_status = 'DUMP_OUT' and ox.ensembl_object_type = 'Translation') join source s on (x.source_id = s.source_id and s.name like 'RefSeq\_peptide%')");

  my $ox_mark_delete_sth = $dbi->prepare("update object_xref set ox_status  = 'MULTI_DELETE' where object_xref_id = ?");

  $transcr_obj_xrefs_sth->execute();
  
  my %change;

  #this hash stores all the translations linked to RefSeq_peptide% xrefs whose transcript is also linked to the paired RefSeq_mRNA - this will be needed to get rid of RefSeq_peptide object xrefs which don't have the additional support of transcripts linked to paired RefSeq_mRNAs; keyed on RefSeq_peptide% accession
  my %RefSeq_pep_translation;

  while(my ($translation_id, $pep_source_id, $pep_accession, $query_identity, $target_identity) = $transcr_obj_xrefs_sth->fetchrow_array() ){
  
      #check if translation is linked to the paired RefSeq peptide

      if ($translation_id) {

	  $ox_translation_sth->execute($translation_id, $pep_source_id, $pep_accession);
	  my ($transl_object_xref_id, $xref_id) = $ox_translation_sth->fetchrow_array();

	  #if it's already linked we don't have to do anything

	  if (!$transl_object_xref_id) {

		  #add a new object xref 
		  $xref_sth->execute($pep_accession, $pep_source_id);
		  ($xref_id) = $xref_sth->fetchrow_array();
		  if (!$xref_id) {
		      die("Xref not found for accession $pep_accession source_id $pep_source_id");
		  }
		  $ox_insert_sth->execute($xref_id, $translation_id, "Translation") || die "Could not insert object xref $object_xref_id: xref_id $xref_id, translation_id $translation_id" ;
                  $get_object_xref_id_sth->execute($xref_id, $translation_id, 'Translation');
                  $object_xref_id = ($get_object_xref_id_sth->fetchrow_array())[0];
		  $xref_update_sth->execute($xref_id)|| die "Could not update xref_id $xref_id";

		  if ($query_identity && $target_identity) {
		      $identity_update_sth->execute($object_xref_id, $query_identity, $target_identity);
		  }
		  	  
		  $change{'translation object xrefs added'}++;
		  $transl_object_xref_id = $object_xref_id;
	      
	  }

	  if ($transl_object_xref_id) {
	      push @{$RefSeq_pep_translation{$pep_accession}}, $translation_id;
	  }

      }

  }

  $transcr_obj_xrefs_sth->finish();
  $ox_translation_sth->finish();
  $ox_insert_sth->finish();
  $xref_update_sth->finish();
  $identity_update_sth->finish();
  $xref_sth->finish();

  #go through RefSeq_peptide% object_xrefs 
  $transl_object_xrefs_sth->execute();
  while (my ($translation_object_xref_id, $translation_id, $pep_accession, $transcript_id) =  $transl_object_xrefs_sth->fetchrow_array() ) {

      if (exists($RefSeq_pep_translation{$pep_accession}) ) {

	  my $found = 0;
	  foreach my $tr_id (@{$RefSeq_pep_translation{$pep_accession}}) {
	      if ($tr_id == $translation_id) {
		  $found = 1;
	      }
	  }
	  if (!$found) {
	      #this translations's transcript is not matched with the paired RefSeq_mRNA%,
	      #change the status to 'MULTI_DELETE'
	      $ox_mark_delete_sth->execute($translation_object_xref_id) || die("Failed to update status to 'MULTI_DELETE for object_xref_id $translation_object_xref_id");
              # Process all dependent xrefs as well
              $self->process_dependents($translation_object_xref_id, $translation_id, $transcript_id, $dbi);

	      $change{'translation object xrefs removed'}++;
	  }
	  
      }
  }
 
  $transl_object_xrefs_sth->finish();
  $ox_mark_delete_sth->finish();
  
  foreach my $key (keys %change){
      print "$key:\t".$change{$key}."\n" if($self->verbose);
  }

  #update process status
  my $sth_stat = $dbi->prepare("insert into process_status (status, date) values('processed_pairs',now())");
  $sth_stat->execute();
  $sth_stat->finish;
}

sub process_dependents {
  my ($self, $translation_object_xref_id, $translation_id, $transcript_id, $dbi) = @_;

  my $dep_tl_sth        = $dbi->prepare("select distinct dependent_ox.object_xref_id from object_xref master_ox, object_xref dependent_ox, xref dependent, xref master, dependent_xref dx where dependent.xref_id = dx.dependent_xref_id and master.xref_id = dx.master_xref_id and dependent.xref_id = dependent_ox.xref_id and master.xref_id = master_ox.xref_id and master_ox.object_xref_id = ? and dependent_ox.master_xref_id = master.xref_id and dependent_ox.ensembl_id = ? and dependent_ox.ensembl_object_type = 'Translation' and dependent_ox.ox_status = 'DUMP_OUT' ");
  my $dep_tr_sth       =  $dbi->prepare("select distinct dependent_ox.object_xref_id from object_xref master_ox, object_xref dependent_ox, xref dependent, xref master, dependent_xref dx where dependent.xref_id = dx.dependent_xref_id and master.xref_id = dx.master_xref_id and dependent.xref_id = dependent_ox.xref_id and master.xref_id = master_ox.xref_id and master_ox.object_xref_id = ? and dependent_ox.master_xref_id = master.xref_id and dependent_ox.ensembl_id = ? and dependent_ox.ensembl_object_type = 'Transcript' and dependent_ox.ox_status = 'DUMP_OUT' ");
  my $ox_dx_delete_sth = $dbi->prepare("update object_xref set ox_status = 'MULTI_DELETE' where object_xref_id = ?");

  my @master_object_xrefs;
  my $new_master_object_xref_id;
  push @master_object_xrefs, $translation_object_xref_id;
  my %master_object_xref_id;
  $master_object_xref_id{$translation_object_xref_id} = 1;

  while (my $master_object_xref_id = pop(@master_object_xrefs)) {
    my $dependent_object_xref_id;
    $dep_tl_sth->execute($master_object_xref_id, $translation_id);
    $dep_tl_sth->bind_columns(\$dependent_object_xref_id);
    while ($dep_tl_sth->fetch()) {
      $ox_dx_delete_sth->execute($dependent_object_xref_id);
      if (!defined $master_object_xref_id{$dependent_object_xref_id}) {
        $master_object_xref_id{$dependent_object_xref_id} = 1;
        push @master_object_xrefs, $dependent_object_xref_id;
      }
    }
    $dep_tr_sth->execute($master_object_xref_id, $transcript_id);
    $dep_tr_sth->bind_columns(\$dependent_object_xref_id);
    while ($dep_tr_sth->fetch()) {
      $ox_dx_delete_sth->execute($dependent_object_xref_id);
      if (!defined $master_object_xref_id{$dependent_object_xref_id}) {
        $master_object_xref_id{$dependent_object_xref_id} = 1;
        push @master_object_xrefs, $dependent_object_xref_id;
      }
    }
  }
}

1;
