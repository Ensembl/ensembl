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

package XrefMapper::CoreInfo;
use strict;
use warnings;

use vars '@ISA';
@ISA = qw{ XrefMapper::BasicMapper };

use XrefMapper::BasicMapper;

use Cwd;
use DBI;
use File::Basename;
use IPC::Open3;

# Get info from the core database.

# Need to load tables:-
#
# gene_transcript_translation 
# gene_stable_id
# transcript_stable_id
# translation_stable_id


sub new {
  my($class, $mapper) = @_;

  my $self ={};
  bless $self,$class;
  $self->core($mapper->core);
  $self->xref($mapper->xref);
  $self->verbose($mapper->verbose);
  return $self;
}



sub get_core_data {
  my $self = shift;

  # gene_transcript_translation 
  # gene_stable_id
  # transcript_stable_id
  # translation_stable_id


  # load table gene_transcript_translation 

  $self->load_gene_transcript_translation();

  # load table xxx_stable_id

  $self->load_stable_ids();


  my $sth = $self->xref->dbc->prepare("insert into process_status (status, date) values('core_data_loaded',now())");
  $sth->execute();
  $sth->finish;


  return;
}


sub load_gene_transcript_translation{
  my ($self) = shift;
  
  my $ins_sth =  $self->xref->dbc->prepare("insert ignore into gene_transcript_translation (gene_id, transcript_id, translation_id) values (?, ?, ?)"); 

  my $sql = "select tn.gene_id, tn.transcript_id, tl.translation_id from transcript tn left join translation tl on tl.transcript_id = tn.transcript_id";
  my $sth = $self->core->dbc->prepare($sql);
  $sth->execute();
  my  ($gene_id, $transcript_id, $translation_id);
  $sth->bind_columns(\$gene_id, \$transcript_id, \$translation_id); 
  while($sth->fetch()){
    $ins_sth->execute($gene_id, $transcript_id, $translation_id);
  }
  $ins_sth->finish;
  $sth->finish;
  return;
}

sub load_stable_ids{
  my ($self) = shift;

  my ($id, $stable_id, $biotype);
  foreach my $table (qw(gene translation)){
 
    my $sth = $self->core->dbc->prepare("select ".$table."_id, stable_id from ".$table);
    my $ins_sth = $self->xref->dbc->prepare("insert ignore into ".$table."_stable_id (internal_id, stable_id) values(?, ?)");
    $sth->execute();
    $sth->bind_columns(\$id, \$stable_id);
    while($sth->fetch){
      $ins_sth->execute($id, $stable_id);
    }
    $ins_sth->finish;
    $sth->finish;
  }

  #populate transcript_stable_id table incuding the biotype column
  my $table = "transcript";
  my $sth = $self->core->dbc->prepare("select ".$table."_id, stable_id, biotype from ".$table);
  my $ins_sth = $self->xref->dbc->prepare("insert ignore into ".$table."_stable_id (internal_id, stable_id, biotype) values(?, ?, ?)");
  $sth->execute();
  $sth->bind_columns(\$id, \$stable_id, \$biotype);
  while($sth->fetch){
    $ins_sth->execute($id, $stable_id, $biotype);
  }
  $ins_sth->finish;
  $sth->finish;

  return;
}
1;
