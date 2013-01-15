
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
use Bio::EnsEMBL::Utils::Exception;
use DBI qw( :sql_types );

=head2 fetch_all_Groups

  Arg[1]      : (optional) String - type of group
  Description : Fetches all the alt-allele groups, creates objects to represent
                them and returns them in a list
                Multispecies support is triggered by the is_multispecies flag
                and species_id of the DBAdaptor.         
  Returntype  : Listref of Bio::EnsEMBL::AltAlleleGroup
=cut

sub fetch_all_Groups {
    my $self = shift;
    my $type = shift;

    $type = uc($type) if ($type);
    
    my @group_list;
    my @members;
    
    my $species_id;
    my $get_all_sql;
    if ($self->db->is_multispecies()) {
        # multispecies databases must be restricted in their treatment
        $species_id = $self->db->species_id;
        $get_all_sql = q(
            SELECT DISTINCT alt_allele_id FROM alt_allele a
            JOIN (gene g, seq_region s, coord_system c)
            ON (
                c.coord_system_id = s.coord_system_id 
            AND s.seq_region_id = g.seq_region_id
            AND g.gene_id = a.gene_id
            )
            WHERE c.species_id = ? 
        );
        if ($type) {$get_all_sql .= q( AND a.type = ? )}
    } else {
        $get_all_sql = q(SELECT DISTINCT alt_allele_id FROM alt_allele);
        if ($type) {$get_all_sql .= q( WHERE type = ?);}
    }
    
    my $sth = $self->prepare($get_all_sql);
    
    my $x = 1;
    if ($self->db->is_multispecies()) {
        $sth->bind_param($x,$species_id, SQL_INTEGER);
        $x++;
    }
    $sth->bind_param($x,$type, SQL_VARCHAR) if ($type);
    
    eval { $sth->execute() };
    if ($@) {
        throw("Query error in AltAlleleGroupAdaptor: $@");
    }
    
    
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
        push @members,[$gene_id,$is_ref,$type];
    }
    if ($allele_id && scalar(@members) > 0) {
        my $aag = Bio::EnsEMBL::AltAlleleGroup->new(
            -dbID => $allele_id,
            -MEMBERS => \@members,
            -ADAPTOR => $self,
        );
        return $aag;
    }
    return;
}

sub fetch_Group_by_Gene_dbID {
    my $self = shift;
    my $gene_id = shift;
    
    my $gene_id_sql = q(
        SELECT alt_allele_id FROM alt_allele
        WHERE gene_id = ?
    );
    my $sth = $self->prepare($gene_id_sql);
    $sth->bind_param(1,$gene_id, SQL_INTEGER);
    
    my $group_id;
    $sth->execute();
    $sth->bind_columns(\$group_id);
    $sth->fetch;
    
    if (!$@ && $group_id) {
        return $self->fetch_Group_by_id($group_id);
    }
}

=head2 store

  Arg[0]     : Bio::EnsEMBL::AltAlleleGroup
  Description: Used for persisting new groups to the database.
               It updates the dbID of the object handed to it to match the
               database.
  Returntype : Alt Allele Group id

=cut

sub store {
    my $self = shift;
    my $allele_group = shift;
    
    if (ref($allele_group) ne "Bio::EnsEMBL::AltAlleleGroup") {
        throw ("Can only store Bio::EnsEMBL::AltAlleleGroup objects.");   
    } else {
        if ($allele_group->size < 2) {
            warning('At least 2 genes must be provided to construct alternative alleles. Ignoring.');
            return;
        }
        
        # first dbID goes in as undef for new records. Mysql swallows this and autoincrements
        my $dbID = $allele_group->dbID;
        
        my $sth = $self->prepare("INSERT INTO alt_allele (alt_allele_id, gene_id, is_ref, type) VALUES (?,?,?,?)");
        
        foreach my $allele (@{ $allele_group->get_all_members() }) {
            $sth->bind_param(1, $dbID, SQL_INTEGER);
            $sth->bind_param(2, $allele->[0], SQL_INTEGER);
            $sth->bind_param(3, $allele->[1], SQL_BOOLEAN);
            $sth->bind_param(4, $allele->[2], SQL_VARCHAR);
            $sth->execute();
            $dbID = $sth->{'mysql_insertid'}; # all alleles get added to the same alt_allele_id group
        }
        if ($@) {throw ("Problem inserting new AltAlleleGroup into database: $@");}
        $sth->finish;
        
        $allele_group->dbID($dbID);
        return $dbID;
    }
}

=head2 update

  Arg [1]    : AltAlleleGroup 
  Description: Removes the existing DB record of an AltAlleleGroup and stores 
               the altered version.
  Returntype : Integer - the return value of the store method, viz. whether the
               insert was successful.
=cut

sub update {
    my $self = shift;
    my $allele_group = shift;
    
    if (ref($allele_group) ne "Bio::EnsEMBL::AltAlleleGroup") {
        throw ("Can only update Bio::EnsEMBL::AltAlleleGroup objects.");   
    } elsif (!$allele_group->dbID) {
        throw ("Cannot update AltAlleleGroup with missing dbID. AltAlleleGroups should be fetched from the DB prior to updating them");
    }    
    $self->remove($allele_group);
    return $self->store($allele_group);
}

=head2 remove

  Arg [1]    : The AltAlleleGroup to delete or dbID thereof 
  Example    : $aaga->remove($alt_allele_group);
  Description: This removes an AltAlleleGroup from the database. 
  Exceptions : None
  
=cut

sub remove {
    my $self = shift;
    my $allele_group = shift;
    
    my $sth = $self->prepare(
        "DELETE FROM alt_allele WHERE alt_allele_id = ?"
    );
    if (ref($allele_group) eq "Bio::EnsEMBL::AltAlleleGroup") {
        $sth->bind_param(1,$allele_group->dbID,SQL_INTEGER);
    } else {
        $sth->bind_param(1,$allele_group,SQL_INTEGER);
    }
    $sth->execute;
    $sth->finish;
}


1;