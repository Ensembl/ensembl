# EnsEMBL module for ArchiveStableIdAdaptor
# Copyright EMBL-EBI/Sanger center 2003
#
#
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::ArchiveStableIdAdaptor

=head1 SYNOPSIS

ArchiveStableIdAdaptor does all SQL to create ArchiveStableIds and works of 
 stable_id_event
 mapping_session
 archive_peptide
 archive_gene

tables inside EnsEMBL core.

=head1 DESCRIPTION
  

=cut



package Bio::EnsEMBL::ArchiveStableIdAdaptor;


use warnings;
use strict;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use vars qw(@ISA);


@ISA = qw( Bio::EnsEMBL::DBSQL::BaseAdaptor );


sub fetch_by_stable_id_version {
}

sub fetch_by_stable_id_dbname {
}

sub fetch_all_pre_by_stable_id_dbname {
}

sub fetch_all_succ_by_stable_id_dbname {
}



# list oldest to newest dbname in the system
sub list_dbnames {
  my $self = shift;
  
  if( ! defined $self->{'dbnames'} ) {
    my $sql = qq {
      SELECT old_db_name, new_db_name
        FROM mapping_session
       ORDER BY created
     };
    my $sth = $self->prepare( $sql );
    $sth->execute();
    my ( $old_db_name, $new_db_name );
    
    my $first = 1;

    $sth->bind_columns( \$old_db_name, \$new_db_name );
    while( $sth->fetch() ) {
      if( $old_db_name eq "ALL" ) {
	next;
      }
      if( $first ) {
	$self->{'dbnames'} = [];
	push( @{$self->{'dbnames'}}, $old_db_name );
	$first = 0;
      }
      push( @{$self->{'dbnames'}}, $new_db_name );
    }
  }

  return $self->{'dbnames'};
}

# find the petide sequence if you can ...
sub _peptide {
  my $self = shift;
  my $arch_id = shift;

  my $sql = qq {
    SELECT peptide_sequence
      FROM archive_peptide
     WHERE translation_stable_id = ?
       AND version = ?
     };
}


# given an ArchiveStableId thats missing the version but has the db_name
# this should fill in the version (not trivial...)
sub _lookup_version {
  my $self = shift;
  my $arch_id = shift;


  # lookup in archive_gene join with mapping session
  # if type is not set, have to 
}


  
