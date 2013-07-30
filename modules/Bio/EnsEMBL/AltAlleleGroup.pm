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

  Bio::EnsEMBL::AltAlleleGroup

=head1 SYNOPSIS

  use Bio::EnsEMBL::AltAlleleGroup;
  use Bio::EnsEMBL::DBSQL::AltAlleleGroupAdaptor;
  
  my $aag_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor("Human","core","AltAlleleGroup");
  
  # For a known Gene, find the reference alternative allele
  my $aag = $aag_adaptor->fetch_Group_by_dbID($gene->dbID);
  my $reference_gene = $aag->get_representative_Gene;
  
  # Get a list of AltAlleleGroups
  my $list = $aag_adaptor->fetch_all_Groups_by_type('HAS_CODING_POTENTIAL');
  $list = $aag_adaptor->fetch_all_Groups();
  
  while ($aag = shift @$list) {
      $aag->get_all_Genes;
      # Do your important things ...
  }
  
  # Creating and editing an AltAlleleGroup
  
  my $type_flags = [qw(IS_MOST_COMMON_ALLELE AUTOMATICALLY_ASSIGNED)];
  
  $aag = Bio::EnsEMBL::AltAlleleGroup->new(
     -MEMBERS => [ [$gene_id,$type_flags ] ],
  );
  $aag->remove_all_members;
  $aag->add_member([$gene_id,$type_flags]);
  
  my $dbID = $aag_adaptor->store($aag);
  

=head1 DESCRIPTION

    Alt allele groups keep track of which alleles are tied to a particular Gene
    They allow related genes to be located. This class allows fetching of both
    IDs and fully fledged Gene objects.
    
    AltAlleleGroup members are assigned types to differentiate them by their
    origin. These types are set as flags, allowing you to select the union of
    types as well as by individual ones.
        
    No flags set denotes a situation of no information.
    Valid flags are as follows:
    'IS_REPRESENTATIVE',
    'IS_MOST_COMMON_ALLELE',
    'IN_CORRECTED_ASSEMBLY',
    'HAS_CODING_POTENTIAL',
    'IN_ARTIFICIALLY_DUPLICATED_ASSEMBLY',
    'IN_SYNTENIC_REGION',
    'HAS_SAME_UNDERLYING_DNA_SEQUENCE',
    'IN_BROKEN_ASSEMBLY_REGION',
    'IS_VALID_ALTERNATE',
    'SAME_AS_REPRESENTATIVE',
    'SAME_AS_ANOTHER_ALLELE',
    'MANUALLY_ASSIGNED',
    'AUTOMATICALLY_ASSIGNED'
=cut

package Bio::EnsEMBL::AltAlleleGroup;

use strict;
use warnings;
#use constant {
#    IS_REPRESENTATIVE => 1,
#    IS_MOST_COMMON_ALLELE => 2,
#    IN_CORRECTED_ASSEMBLY => 3,
#    HAS_CODING_POTENTIAL => 4,
#    IN_ARTIFICIALLY_DUPLICATED_ASSEMBLY => 5,
#    IN_SYNTENIC_REGION => 6,
#    HAS_SAME_UNDERLYING_DNA_SEQUENCE => 7,
#    IN_BROKEN_ASSEMBLY_REGION => 8,
#    IS_VALID_ALTERNATE => 9,
#    SAME_AS_REPRESENTATIVE => 10,
#    SAME_AS_ANOTHER_ALLELE => 11,
#    MANUALLY_ASSIGNED => 12,
#    AUTOMATICALLY_ASSIGNED => 13,
#};


use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(warning throw);

use base qw/Bio::EnsEMBL::Storable/;

=head2 new

  Arg [-MEMBERS]: A list reference of [gene_id,type_flags]
                : gene_id is a dbID for Gene (consistent only within one release)
                : type_flags is a hash ref of attributes for this member
  Example    : $aag = Bio::EnsEMBL::AltAlleleGroup->new(
                   -MEMBERS => [ [1,{$type} ], [2,{$other_type}],[3,{$type}],
               );
  Description: Creates a new alt-allele group object
  Returntype : Bio::EnsEMBL::AltAlleleGroup
  Exceptions : none
  Caller     : general

=cut

sub new {
    my $caller = shift;

    my $class = ref($caller) || $caller;
    my $self = $class->SUPER::new(@_);
    my ( $list ) = rearrange( [ 'MEMBERS'], @_ );
    
    $self->{'MEMBERS'} = $list;
    
    return $self;
}

=head2 add_member

  Arg [1]     : Gene dbID
  Arg [2]     : Type List, used for assigning type flags of this member, see Description above
  Description : Adds a record of one new member to the AltAlleleGroup. Once a
                change is made, this must be persisted to the database with
                AltAlleleGroupAdaptor->store or ->update
  Example     : $aag->add_member(1040032,$types_hash);
                $aaga->update($aag); # updating the whole group is necessary.
=cut

sub add_member {
    my $self = shift;
    my ($gene_id,$type_hash) = @_;
    
    my $members = $self->{'MEMBERS'};
    push @$members,[$gene_id,$type_hash];
    $self->{'MEMBERS'} = $members;
    return;
}

sub get_all_members_with_type {
    my $self = shift;
    my $type = shift;
    
    my @filtered_members;
    my $members = $self->{'MEMBERS'};
    foreach my $member (@$members) {
        if (exists($member->[1]->{$type})) {
            push @filtered_members,$member;
        }
    }
    return \@filtered_members;
}

sub attribs {
    my $self = shift;
    my $member_id = shift;
    
    
}

=head2 remove_all_members

  Description: Remove members from this object, but NOT the database. See
               AltAlleleGroupAdaptor->remove()
               Use in conjunction with add_member if members need to be altered
=cut

sub remove_all_members {
    my $self = shift;
    $self->{'MEMBERS'} = [];
    return;
}

=head2 rep_Gene_id

  Arg[1]     : Optional - set a new representative Gene id for the group
  Description: Reports or sets the representative Gene for this AltAlleleGroup
               If you wish to remove the representative status of all genes without
               setting a new one, see unset_rep_Gene_id
  Returntype : Integer or undef if none set
=cut

sub rep_Gene_id {
    my $self = shift;
    my $new_id = shift;
    my $list = $self->{'MEMBERS'};
    my $change;
    
    foreach my $allele (@$list) {
        my ($gene_id,$type) = @$allele;
        if (exists($type->{IS_REPRESENTATIVE}) && !defined($new_id) ) {
            return $gene_id;
        }
        
        if ($new_id) {
            unless ($gene_id == $new_id) {delete($allele->[1]->{IS_REPRESENTATIVE})}
            else {
                $allele->[1]->{IS_REPRESENTATIVE} = 1; 
                $change = $new_id;
            }
        }
    }
    
    if ($change) {
        $self->{'MEMBERS'} = $list;
        return $new_id;
    } elsif ($new_id && !$change) {
        my $db_id = $self->dbID() || 'unknown';
        throw("Requested representative gene ID was not set because it is not in this AltAlleleGroup, ID $db_id");
    }
    else {
        warning("No representative allele currently set for this AltAlleleGroup");
        return;
    }
}

=head2 unset_rep_Gene_id

  Description: Removes the representative Gene flag from this AltAlleleGroup.
               This action is not possible through rep_Gene_id due to
               validation of inputs.
  Returntype : 

=cut

sub unset_rep_Gene_id {
    my $self = shift;
    my $list = $self->{'MEMBERS'};
    
    foreach my $allele (@$list) {
        delete($allele->[1]->{IS_REPRESENTATIVE});
    }
    $self->{'MEMBERS'} = $list;
    return;
}

=head2 get_all_Gene_ids

  Arg[1]      : Boolean - Do not include representative gene in list of ids.
  Description : fetches all the Gene dbIDs within the allele group. It can also
                be used to list those ids that are not the representative Gene.
                
  Returntype  : listref of gene dbIDs

=cut

sub get_all_Gene_ids {
    my $self = shift;
    my $all_but_rep = shift;
    my $list = $self->{'MEMBERS'};
    
    my @gene_ids;
    
    foreach my $allele (@$list) {
        my ($gene_id,$type) = @$allele;
        if ($all_but_rep && $type->{IS_REPRESENTATIVE}) {next;} 
        push @gene_ids,$gene_id;
    }
    return \@gene_ids;
}

sub get_representative_Gene {
    my $self = shift;
    my $ga = $self->adaptor->db->get_GeneAdaptor;
    
    return $ga->fetch_by_dbID($self->rep_Gene_id);
}

sub get_all_Genes {
    my $self = shift;
    my $all_but_rep = shift; # falls through to get_all_Gene_ids
    my $gene_ids = $self->get_all_Gene_ids($all_but_rep);
    my $genes;
    my $ga = $self->adaptor->db->get_GeneAdaptor;
    $genes = $ga->fetch_all_by_dbID_list($gene_ids);
    
    return $genes;
}

sub size {
    my $self = shift;
    my $list = $self->{'MEMBERS'};
    return scalar(@$list);
}


=head2 get_all_members
  Description: Retrieves all of the information about all members.
  Returntype : Listref of id and type list: [gene_id,type]
  Caller     : AltAlleleGroupAdaptor->store
=cut

sub get_all_members {
    my $self = shift;
    my $members = $self->{'MEMBERS'};
    return $members;
}


1;