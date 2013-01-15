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

  Bio::EnsEMBL::AltAlleleGroup

=head1 SYNOPSIS

  use Bio::EnsEMBL::AltAlleleGroup;
  use Bio::EnsEMBL::DBSQL::AltAlleleGroupAdaptor;
  
  my $aag_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor("Human","core","AltAlleleGroup");
  
  $aag_adaptor->fetch_all_Groups_by_type('MANUAL');

=head1 DESCRIPTION

    Alt allele groups keep track of which alleles are tied to a particular Gene
    They allow related genes to be located.
=cut

package Bio::EnsEMBL::AltAlleleGroup;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(warning throw);

use base qw/Bio::EnsEMBL::Storable/;

=head2 new

  Arg [-MEMBERS]: A list reference of [gene_id,is_ref,type]
                : gene_id is a dbID for Gene (consistent only within one release)
                : is_ref is a boolean flag denoting the reference allele
                : type is a string for descriptive purposes
  Example    : $aag = Bio::EnsEMBL::AltAlleleGroup->new(...);
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

sub add_member {
    my $self = shift;
    my ($gene_id,$is_ref,$type) = @_;
    
    my $members = $self->{'MEMBERS'};
    push @$members,[$gene_id,$is_ref,$type];
    $self->{'MEMBERS'} = $members;
    return 1;
}

sub remove_all_members {
    my $self = shift;
    $self->{'MEMBERS'} = [];
    return 1;
}

=head2 ref_Gene_id

  Arg[1]     : Optional - set a new reference Gene id for the group
  Description: Reports or sets the reference Gene for this AltAlleleGroup
               If you wish to remove the reference status of all genes without
               setting a new one, see unset_ref_Gene_id
  Returntype : Integer or undef if none set
=cut

sub ref_Gene_id {
    my $self = shift;
    my $new_id = shift;
    my $list = $self->{'MEMBERS'};
    my $change;
    
    foreach my $allele (@$list) {
        my ($gene_id,$is_ref,$type) = @$allele;
        if ($is_ref && !defined($new_id) ) {
            return $gene_id;
        }
        
        if ($new_id) {
            unless ($gene_id == $new_id) {$allele->[1] = 0}
            else {
                $allele->[1] = 1; 
                $change = $new_id;
            }
        }
    }
    
    if ($change) {
        $self->{'MEMBERS'} = $list;
        return $new_id;
    } elsif ($new_id && !$change) {
        throw("Requested reference gene ID was not set because it is not in this AltAlleleGroup");
    }
    else {
        warning("No reference allele currently set for this AltAlleleGroup");
        return;
    }
}

sub unset_ref_Gene_id {
    my $self = shift;
    my $list = $self->{'MEMBERS'};
    
    foreach my $allele (@$list) {
        $allele->[1] = 0;
    }
    $self->{'MEMBERS'} = $list;
    return 1;
}

=head2 get_all_Gene_ids

  Arg[1]      : Boolean - Do not include reference gene in list of ids.
  Description : fetches all the Gene dbIDs within the allele group. It can also
                be used to list those ids that are not reference.
                
  Returntype  : listref of gene dbIDs

=cut

sub get_all_Gene_ids {
    my $self = shift;
    my $all_but_ref = shift;
    my $list = $self->{'MEMBERS'};
    
    my @gene_ids;
    
    foreach my $allele (@$list) {
        my ($gene_id,$is_ref,$type) = @$allele;
        if ($all_but_ref && $is_ref) {next;} 
        push @gene_ids,$gene_id;
    }
    return \@gene_ids;
}

sub get_ref_Gene {
    my $self = shift;
    my $ga = $self->adaptor->db->get_GeneAdaptor;
    
    return $ga->fetch_by_dbID($self->ref_Gene_id);
}

sub get_all_Genes {
    my $self = shift;
    my $all_but_ref = shift; # falls through to get_all_Gene_ids
    my $gene_ids = $self->get_all_Gene_ids($all_but_ref);
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
  Returntype : Listref of triplets: [gene_id,is_ref,type]
  Caller     : AltAlleleGroupAdaptor->store
=cut

sub get_all_members {
    my $self = shift;
    my $members = $self->{'MEMBERS'};
    return $members;
}

1;