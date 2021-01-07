=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

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



=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::DBSQL::AltAlleleGroupAdaptor - Adaptor for the manipulation of
Alternative allele groupings

=head1 SYNOPSIS

  use Bio::EnsEMBL::AltAlleleGroup;
  use Bio::EnsEMBL::DBSQL::AltAlleleGroupAdaptor;
  
  my $aag_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor("Human","core","AltAlleleGroup");
  
  # For a known Gene, find the reference alternative allele
  my $aag = $aag_adaptor->fetch_by_gene_id($gene->dbID);
  my $reference_gene = $aag->get_ref_Gene;
  
  # Get a list of AltAlleleGroups
  my $list = $aag_adaptor->fetch_all_('IS_REPRESENTATIVE');
  $list = $aag_adaptor->fetch_all();
  
  my $dbID = $aag_adaptor->store($aag);
  
  $aag = $aag_adaptor->fetch_by_dbID($dbID);
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
use Bio::EnsEMBL::Utils::Exception qw/throw warning/;
use Bio::EnsEMBL::Utils::Scalar qw/assert_ref/;
use DBI qw( :sql_types );

=head2 fetch_all

  Arg[1]      : (optional) String - type of group
                Restrict group fetches to just one type. Technically it selects 
                out mixed-annotation groups where a single member contains that type.
  Description : Fetches all the alt-allele groups, creates objects to represent
                them and returns them in a list. Specifying a group type 
                identifies all groups containing a member of this type. It 
                does not filter out the other members
                
                Multispecies support is triggered by the is_multispecies flag
                and species_id of the DBAdaptor.
  Returntype  : ArrayRef of Bio::EnsEMBL::AltAlleleGroup

=cut

sub fetch_all {
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
                SELECT DISTINCT alt_allele_group_id
                FROM alt_allele        a
                JOIN gene              g  ON g.gene_id         = a.gene_id
                JOIN seq_region        s  ON s.seq_region_id   = g.seq_region_id
                JOIN coord_system      c  ON c.coord_system_id = s.coord_system_id
                JOIN alt_allele_attrib b  ON a.alt_allele_id   = b.alt_allele_id
                WHERE c.species_id = ? AND b.attrib = ?
            );
        } else {
            $get_all_sql = q(
                SELECT DISTINCT alt_allele_group_id
                FROM alt_allele a
                JOIN gene              g  ON g.gene_id         = a.gene_id
                JOIN seq_region        s  ON s.seq_region_id   = g.seq_region_id
                JOIN coord_system      c  ON c.coord_system_id = s.coord_system_id
                WHERE c.species_id = ?
            );
        }
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
        my $aag = $self->fetch_by_dbID($group_id);
        push @group_list, $aag;
    }
    $sth->finish;
    return \@group_list;
}

=head2 fetch_by_dbID

  Arg[1]      : AltAlleleGroup dbID.
  Description : Creates and returns an AltAlleleGroup for the given group id
                
  Returntype  : Bio::EnsEMBL::AltAlleleGroup

=cut

sub fetch_by_dbID {
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

=head2 fetch_by_gene_id

  Arg[1]      : Integer Gene ID of the member to query by
  Description : Creates and returns an AltAlleleGroup which contains
                the specified gene member                
  Returntype  : Bio::EnsEMBL::AltAlleleGroup

=cut

sub fetch_by_gene_id {
    my ($self, $gene_id) = @_;

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
        return $self->fetch_by_dbID($group_id);
    }
    return;
}

=head2 store

  Arg[1]     : Bio::EnsEMBL::AltAlleleGroup
  Description: Used for persisting new groups to the database.
               It updates the dbID of the object handed to it to match the
               database.
  Returntype : Integer Alt Allele Group id

=cut

sub store {
    my $self = shift;
    my $allele_group = shift;

    assert_ref($allele_group, 'Bio::EnsEMBL::AltAlleleGroup', 'allele_group');
    if ($allele_group->size < 2) {
        warning('At least 2 genes must be provided to construct alternative alleles. Ignoring.');
        return;
    }
    
    my $helper = $self->dbc()->sql_helper();
    my $dbID = $allele_group->dbID;
    
    my $new_group_sql = 'INSERT INTO alt_allele_group (alt_allele_group_id) VALUES (?)';
    my $existing_group_sql = 'SELECT count(*) FROM alt_allele_group WHERE alt_allele_group_id = ?';

    my $already_exists = $helper->execute_single_result(-SQL => $existing_group_sql, -PARAMS => [[$dbID, SQL_INTEGER]]);
    
    # If the ID is not already there then we need to add one
    if($already_exists == 0) {
        $helper->execute_update(-SQL => $new_group_sql, -CALLBACK => sub {
            my ($sth, $dbh, $rv) = @_;
            if($rv) {
                my $id = $dbh->last_insert_id(undef, undef, 'alt_allele_group', 'alt_allele_group_id');
                $dbID = $id;
            }
            return;
        });
    }
    
    my $sth = $self->prepare("INSERT INTO alt_allele (alt_allele_id, alt_allele_group_id, gene_id) VALUES (?,?,?)");
    my $attrib_sth = $self->prepare("INSERT INTO alt_allele_attrib (alt_allele_id,attrib) VALUES (?,?)");
    my $check_exists_sth = $self->prepare("SELECT alt_allele_id FROM alt_allele WHERE gene_id = ?");

    foreach my $allele (@{ $allele_group->get_all_members() }) {
        my $gene_id = $allele->[0];
        my %flags = %{$allele->[1]};
        my $allele_id;

# Check if gene is not already stored
# Return allele_id if it is
        $check_exists_sth->bind_param(1, $gene_id, SQL_INTEGER);
        $check_exists_sth->execute();
        $check_exists_sth->bind_col(1, \$allele_id);
        if ($check_exists_sth->fetch() ) {
          return $allele_id;
        }
        
        $sth->bind_param(1, undef, SQL_INTEGER);
        $sth->bind_param(2, $dbID, SQL_INTEGER);
        $sth->bind_param(3, $gene_id, SQL_INTEGER);
        my $altered_rows = $sth->execute();
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
    }
    if ($@) {throw ("Problem inserting new AltAlleleGroup into database: $@");}
    $sth->finish;
    $attrib_sth->finish;
    $check_exists_sth->finish;
    
    $allele_group->dbID($dbID);
    
    return $dbID;
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
    assert_ref($allele_group, 'Bio::EnsEMBL::AltAlleleGroup', 'allele_group');
    throw "Cannot update an AltAlleleGroup without a dbID. AltAlleleGroups should be fetched from the DB prior to updating them" if ! $allele_group->dbID();
    my $keep_group = 1;
    $self->remove($allele_group, $keep_group);
    return $self->store($allele_group);
}

=head2 remove

  Arg [1]    : The AltAlleleGroup to remove.
  Arg [2]    : Boolean indicates if the entry in alt_allele_group should be retained or remove. Defaults to removing the entry
  Example    : $aaga->remove($alt_allele_group);
  Description: This removes an AltAlleleGroup from all tables of the database. 
  Exceptions : None
  
=cut

sub remove {
    my ($self, $allele_group, $keep_group) = @_;
    assert_ref($allele_group, 'Bio::EnsEMBL::AltAlleleGroup', 'allele_group');

    my $helper = $self->dbc()->sql_helper();
    my $delete_attribs_sql;
    if ($self->dbc->driver() eq 'mysql') {
      $delete_attribs_sql = q{
        DELETE aaa 
        FROM alt_allele_attrib aaa 
        JOIN alt_allele aa using (alt_allele_id) 
        where alt_allele_group_id =?
      };
    }
    else {
      $delete_attribs_sql = q{
        DELETE FROM alt_allele_attrib WHERE alt_allele_id IN (
            SELECT alt_allele_id FROM alt_allele WHERE alt_allele_group_id = ?
        )
      };
    }
    my $delete_alt_alleles_sql = 'DELETE FROM alt_allele where alt_allele_group_id =?';
    my $delete_group_sql = 'DELETE from alt_allele_group where alt_allele_group_id =?';
    my $params = [[$allele_group->dbID, SQL_INTEGER]];

    $helper->execute_update(-SQL => $delete_attribs_sql, -PARAMS => $params);
    $helper->execute_update(-SQL => $delete_alt_alleles_sql, -PARAMS => $params);
    if(! $keep_group) {
        $helper->execute_update(-SQL => $delete_group_sql, -PARAMS => $params);
    }

    return;
}

sub _tables {
    return (['alt_allele', 'a'], ['alt_allele_group', 'g'], ['alt_allele_attrib', 'b']);
}

1;
