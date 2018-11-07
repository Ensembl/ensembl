use utf8;
package Xref::Schema::Result::HavanaStatus;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Xref::Schema::Result::HavanaStatus

=cut

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 TABLE: C<havana_status>

=cut

__PACKAGE__->table("havana_status");

=head1 ACCESSORS

=head2 stable_id

  data_type: 'varchar'
  is_nullable: 1
  size: 128

=head2 status

  data_type: 'enum'
  extra: {list => ["KNOWN","NOVEL","PUTATIVE","PREDICTED","KNOWN_BY_PROJECTION","UNKNOWN"]}
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "stable_id",
  { data_type => "varchar", is_nullable => 1, size => 128 },
  "status",
  {
    data_type => "enum",
    extra => {
      list => [
        "KNOWN",
        "NOVEL",
        "PUTATIVE",
        "PREDICTED",
        "KNOWN_BY_PROJECTION",
        "UNKNOWN",
      ],
    },
    is_nullable => 1,
  },
);

=head1 UNIQUE CONSTRAINTS

=head2 C<status_idx>

=over 4

=item * L</stable_id>

=back

=cut

__PACKAGE__->add_unique_constraint("status_idx", ["stable_id"]);


# Created by DBIx::Class::Schema::Loader v0.07042 @ 2018-10-23 11:58:10
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:2yXKql+v1IkSzSU+H3mrkw


# You can replace this text with custom code or comments, and it will be preserved on regeneration
1;
