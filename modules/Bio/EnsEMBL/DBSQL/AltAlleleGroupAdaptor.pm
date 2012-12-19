
=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::DBSQL::AltAlleleGroupAdaptor - Adaptor for the manipulation of
Alternative allele groupings

=head1 SYNOPSIS


=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::DBSQL::AltAlleleGroupAdaptor;

use strict;
use warnings;

use base qw/Bio::EnsEMBL::DBSQL::BaseAdaptor/;

use Bio::EnsEMBL::AltAlleleGroup;
use DBI qw( :sql_types );

=head2 fetch_all_Groups

  Arg[1]      : (optional) String - type of group
  Description : Fetches all the alt-allele groups, creates objects to represent
                them and returns them in a list          
  Returntype  : Listref of Bio::EnsEMBL::AltAlleleGroup
=cut

sub fetch_all_Groups {
    my $self = shift;
    my $type = shift;

    $type = uc($type);
    
    my @group_list;
    my @members;
    
    my $get_all_sql = q(
        SELECT DISTINCT alt_allele_id FROM alt_alelle
    );
    
    if ($type) {$get_all_sql .= qq(WHERE type = ?);}
    
    my $sth = $self->prepare($get_all_sql);
    $sth->bind_param(1,$type, SQL_INTEGER) if ($type);
    $sth->execute();
    
    my ($allele_id, @alt_alleles);    
    $sth->bind_columns( \$allele_id );
    
    while ( $sth->fetch() ) {
        my $aag = $self->fetch_Group_by_id($allele_id);
        push @group_list, $aag;
    }
    return \@group_list;
}

=head2 fetch_all_Groups_by_type

  Arg[1]      : String - type of group
  Description : Convenience method for restricting group fetches to just one
                type. Technically it selects out mixed-annotation groups where 
                a single member contains that type.       
  Returntype  : Listref of Bio::EnsEMBL::AltAlleleGroup
=cut

sub fetch_all_Groups_by_type {
    my $self = shift;
    my $type = shift; # refers to alt_allele.type
    
    my $group_list = $self->fetch_all_Groups($type);
    return $group_list;
}

=head2 fetch_Group_by_id

  Arg[1]      : Alt-allele db ID.
  Description : Creates and returns an AltAlleleGroup for the given alt-allele
                
  Returntype  : Bio::EnsEMBL::AltAlleleGroup

=cut

sub fetch_Group_by_id {
    my $self = shift;
    my $allele_id = shift;
    
    my @members;
    
    my $get_alt_allele_sql = q(
        SELECT gene_id, is_ref, type FROM alt_allele
        WHERE alt_allele_id = ?
    );
    my $sth = $self->prepare($get_alt_allele_sql);
    
    $sth->bind_param(1,$allele_id, SQL_INTEGER);
    
    $sth->execute();
    my ($gene_id, $is_ref, $type);
    $sth->bind_columns( \($gene_id,$is_ref,$type) );
    while ($sth->fetch()) {
        push @members,\[$gene_id,$is_ref,$type];
    }
    my $aag = Bio::EnsEMBL::AltAlleleGroup->new(
        -MEMBERS => @members,
    );
    return $aag;
}

1;