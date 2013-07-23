
=head1 LICENSE

  Copyright (c) 1999-2013 The European Bioinformatics Institute and
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

  use Bio::EnsEMBL::AltAlleleGroup;
  use Bio::EnsEMBL::DBSQL::AltAlleleGroupAdaptor;
  
  my $aag_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor("Human","core","AltAlleleGroup");
  
  # For a known Gene, find the reference alternative allele
  my $aag = $aag_adaptor->fetch_Group_by_dbID($gene->dbID);
  my $reference_gene = $aag->get_ref_Gene;
  
  # Get a list of AltAlleleGroups
  my $list = $aag_adaptor->fetch_all_Groups_by_type('IS_REPRESENTATIVE');
  $list = $aag_adaptor->fetch_all_Groups();
  
  my $dbID = $aag_adaptor->store($aag);
  
  $aag = $aag_adaptor->fetch_Group_by_id($dbID);
  $aag_adaptor->remove($aag);

=head1 DESCRIPTION

  The AltAlleleGroupAdaptor provides CRUD for AltAlleleGroup objects. It allows
  groups of alleles to be retrieved by group and gene ids.

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
                Specifying a group type identifies all groups containing a
                member of this type. It does not filter out the other members
  Returntype  : Listref of Bio::EnsEMBL::AltAlleleGroup
=cut

sub fetch_all_Groups {
    my $self = shift;
    my $type = shift;

    $type = uc($type) if ($type);
    
    my @group_list = ();
    my @members;
    
    my $species_id;
    my $get_all_sql;
    if ($self->db->is_multispecies()) {
        # multispecies databases must be restricted in their treatment
        $species_id = $self->db->species_id;
        
        if ($type) {
            $get_all_sql = q(
                SELECT DISTINCT alt_allele_group_id FROM alt_allele a
                JOIN (gene g, seq_region s, coord_system c, alt_allele_attrib b)
                ON (
                    c.coord_system_id = s.coord_system_id 
                    AND s.seq_region_id = g.seq_region_id
                    AND g.gene_id = a.gene_id
                    AND a.alt_allele_id = b.alt_allele_id
                )
                WHERE c.species_id = ? AND b.attrib = ?
            );
        }
        $get_all_sql = q(
            SELECT DISTINCT alt_allele_group_id FROM alt_allele a
            JOIN (gene g, seq_region s, coord_system c)
            ON (
                c.coord_system_id = s.coord_system_id 
                AND s.seq_region_id = g.seq_region_id
                AND g.gene_id = a.gene_id
            )
            WHERE c.species_id = ? 
        );
    } else {
        if ($type) {
            $get_all_sql = q(SELECT DISTINCT alt_allele_group_id 
                FROM alt_allele a, alt_allele_attrib b
                WHERE a.alt_allele_id = b.alt_allele_id
                AND b.attrib = ?);
        } else {
            $get_all_sql = q(SELECT DISTINCT alt_allele_group_id FROM alt_allele);
        }
       
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
    
    
    my $group_id;    
    $sth->bind_col(1, \$group_id );
    
    while ( $sth->fetch() ) {
        my $aag = $self->fetch_Group_by_id($group_id);
        push @group_list, $aag;
    }
    $sth->finish;
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
    my $type = shift; # refers to alt_allele_attrib type
    
    my $group_list = $self->fetch_all_Groups($type);
    return $group_list;
}

=head2 fetch_Group_by_id

  Arg[1]      : AltAlleleGroup dbID.
  Description : Creates and returns an AltAlleleGroup for the given group id
                
  Returntype  : Bio::EnsEMBL::AltAlleleGroup

=cut

sub fetch_Group_by_id {
    my $self = shift;
    my $group_id = shift;
    
    my @members;
    
    my $get_alt_allele_sql = q(
        SELECT alt_allele_id, gene_id FROM alt_allele
        WHERE alt_allele_group_id = ? ORDER BY alt_allele_id
    );
    my $sth = $self->prepare($get_alt_allele_sql);
    
    $sth->bind_param(1,$group_id, SQL_INTEGER);
    
    $sth->execute();
    my ($alt_allele_id, $gene_id);
    $sth->bind_columns( \($alt_allele_id,$gene_id) );
    
    my $attrib_fetch = q(
        SELECT attrib FROM alt_allele_attrib WHERE alt_allele_id = ?
    );
    my $attrib_sth = $self->prepare($attrib_fetch);
    my $attrib;
    
    while ($sth->fetch()) {
        # fetch alt_allele attributes
        $attrib_sth->execute($alt_allele_id);
        $attrib_sth->bind_col(1,\$attrib);
        my %attrib_list;
        while ($attrib_sth->fetch) {
            $attrib_list{$attrib} = 1;
        }
        push @members,[$gene_id, \%attrib_list];
    }
    $attrib_sth->finish;
    $sth->finish;
    
    if ($group_id && scalar(@members) > 0) {
        my $aag = Bio::EnsEMBL::AltAlleleGroup->new(
            -dbID => $group_id,
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
        SELECT alt_allele_group_id FROM alt_allele
        WHERE gene_id = ?
    );
    my $sth = $self->prepare($gene_id_sql);
    $sth->bind_param(1,$gene_id, SQL_INTEGER);
    
    my $group_id;
    $sth->execute();
    $sth->bind_col(1,\$group_id);
    $sth->fetch;
    $sth->finish;
    if (!$@ && $group_id) {
        return $self->fetch_Group_by_id($group_id);
    }
    return;
}

=head2 store

  Arg[0]     : Bio::EnsEMBL::AltAlleleGroup
  Description: Used for persisting new groups to the database.
               It updates the dbID of the object handed to it to match the
               database.
  Returntype : Integer Alt Allele Group id

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
        
        my $dbID = $allele_group->dbID;
        
        my $new_group_sth = $self->prepare("INSERT INTO alt_allele_group (alt_allele_group_id) VALUES (?)");
        my $group_sth = $self->prepare("SELECT alt_allele_group_id FROM alt_allele_group WHERE alt_allele_group_id = ?");
        my $altered_rows;
        
        # Do not create a new group ID if one already exists, such as when updating a group.
        my $existing_rows = $group_sth->execute($dbID);
        if ($existing_rows == 0) {
            $altered_rows = $new_group_sth->execute($dbID);
            
            if ($altered_rows > 0) {
                $dbID = $self->last_insert_id(undef,undef,undef,'alt_allele_group');
                $allele_group->dbID($dbID);
            } 
        }
        
        my $sth = $self->prepare("INSERT INTO alt_allele (alt_allele_id, alt_allele_group_id, gene_id) VALUES (?,?,?)");
        my $attrib_sth = $self->prepare("INSERT INTO alt_allele_attrib (alt_allele_id,attrib) VALUES (?,?)");
        
        foreach my $allele (@{ $allele_group->get_all_members() }) {
            my $gene_id = $allele->[0];
            my %flags = %{$allele->[1]};
            
            $sth->bind_param(1, undef, SQL_INTEGER);
            $sth->bind_param(2, $dbID, SQL_INTEGER);
            $sth->bind_param(3, $gene_id, SQL_INTEGER);
            my $altered_rows = $sth->execute();
            my $allele_id;
            if ($altered_rows > 0) {
                $allele_id = $self->last_insert_id(); # all alleles get added to the same alt_allele_id group
            } else {
                throw("Creation of new alt_allele failed: $@");
            }
            
                
            foreach my $flag (keys %flags) {
                $attrib_sth->bind_param(1, $allele_id);
                $attrib_sth->bind_param(2, $flag);
                $attrib_sth->execute();
            }
            
            if (! $dbID) {
                $group_sth->bind_param(1, $dbID);
                $group_sth->execute($allele_id);
            }
        }
        if ($@) {throw ("Problem inserting new AltAlleleGroup into database: $@");}
        $sth->finish;
        $attrib_sth->finish;
        $group_sth->finish;
        
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

  Arg [1]    : The AltAlleleGroup to remove.
  Example    : $aaga->remove($alt_allele_group);
  Description: This removes an AltAlleleGroup from all tables of the database. 
  Exceptions : None
  
=cut

sub remove {
    my $self = shift;
    my $allele_group = shift;
    
    my $group_id;
    if (ref($allele_group) eq "Bio::EnsEMBL::AltAlleleGroup") {
        $group_id = $allele_group->dbID;
    } else {
        throw("Cannot remove a non-AltAlleleGroup.");
    }
    
    my $allele_sth = $self->prepare("SELECT alt_allele_id FROM alt_allele WHERE alt_allele_group_id = ?");
    my $attrib_sth = $self->prepare("DELETE FROM alt_allele_attrib WHERE alt_allele_id IN (". join(',', ('?')x$allele_group->size) .")" );

    $allele_sth->execute($group_id);
    
    my $allele_id;
    $allele_sth->bind_columns(\$allele_id);
    my @ids;
    while ($allele_sth->fetch) {
        push @ids,$allele_id;
        $attrib_sth->bind_param($#ids+1,$ids[$#ids],SQL_INTEGER);
    }
    $attrib_sth->execute();

    my $sth = $self->prepare("DELETE FROM alt_allele WHERE alt_allele_group_id = ?");
    
    $sth->bind_param(1,$group_id,SQL_INTEGER);
    
    $sth->execute;
    $sth->finish;
    $allele_sth->finish;
    $attrib_sth->finish;
}

sub _tables {
    return (['alt_allele', 'a'], ['alt_allele_group', 'g'], ['alt_allele_attrib', 'b']);
}

1;