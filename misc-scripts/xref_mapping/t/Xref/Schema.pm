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


=head1 DESCRIPTION

A DBIC schema base class used by Xref::DB

=cut

package Xref::Schema;

use strict;
use warnings;
use utf8;

use base 'DBIx::Class::Schema';

__PACKAGE__->load_namespaces;

=head2 sqlt_deploy_hook

Set all tables engine to MyISAM

=cut

sub sqlt_deploy_hook {
  my ($self, $sqlt_schema) = @_;

  for my $table ($sqlt_schema->get_tables) {
    $table->extra(mysql_table_type => 'MyISAM');
  }

  return;
}

1;
