use utf8;
package Xref::Schema::Result::Interpro;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Xref::Schema::Result::Interpro

=cut

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 TABLE: C<interpro>

=cut

__PACKAGE__->table("interpro");

=head1 ACCESSORS

=head2 interpro

  data_type: 'varchar'
  is_nullable: 0
  size: 255

=head2 pfam

  data_type: 'varchar'
  is_nullable: 0
  size: 255

=head2 dbtype

  data_type: 'enum'
  extra: {list => ["PROSITE","PFAM","PREFILE","PROFILE","TIGRFAMs","PRINTS","PIRSF","SMART","SSF"]}
  is_nullable: 0

=cut

__PACKAGE__->add_columns(
  "interpro",
  { data_type => "varchar", is_nullable => 0, size => 255 },
  "pfam",
  { data_type => "varchar", is_nullable => 0, size => 255 },
  "dbtype",
  {
    data_type => "enum",
    extra => {
      list => [
        "PROSITE",
        "PFAM",
        "PREFILE",
        "PROFILE",
        "TIGRFAMs",
        "PRINTS",
        "PIRSF",
        "SMART",
        "SSF",
      ],
    },
    is_nullable => 0,
  },
);


# Created by DBIx::Class::Schema::Loader v0.07042 @ 2018-10-23 11:58:10
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:OgvhrYx2Gx3fC3dC4qaSjQ


# You can replace this text with custom code or comments, and it will be preserved on regeneration
1;
