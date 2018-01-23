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


use strict;
use warnings;

package SeqStoreConverter::HomoSapiens;

use SeqStoreConverter::BasicConverter;

use vars qw(@ISA);

@ISA = qw(SeqStoreConverter::BasicConverter);


sub create_attribs {
  my $self = shift;

  #
  # Human clones need their htg phase information copied
  #

  my $source = $self->source();
  my $target = $self->target();
  my $dbh    = $self->dbh();

  $self->SUPER::create_attribs();
  
  $self->debug("HomoSapiens specific: Creating HTG Phase seq_region attribs");

  $dbh->do
    ("INSERT INTO $target.attrib_type( code, name, description ) " .
     "VALUES ('htg_phase', 'HTG Phase', 'High Throughput Genome Phase')");


  $dbh->do
    ("INSERT INTO $target.seq_region_attrib( seq_region_id, attrib_type_id, " .
                                            "value) " .
     "SELECT tmp_cln.new_id, attrib_type.attrib_type_id, cln.htg_phase " .
     "FROM   $target.tmp_cln_map tmp_cln, $target.attrib_type attrib_type, " .
     "       $source.clone cln " .
     "WHERE  cln.clone_id = tmp_cln.old_id " .
     "AND    attrib_type.code = 'htg_phase'");

  return;
}


1;
