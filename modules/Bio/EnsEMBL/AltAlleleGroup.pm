=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

  my %type_flags = ('IS_MOST_COMMON_ALLELE' => '1','AUTOMATICALLY_ASSIGNED' => '1');
  
  $aag = Bio::EnsEMBL::AltAlleleGroup->new(
     -MEMBERS => [ [$gene_id,\%type_flags ] ],
  );
  $aag->remove_all_members;
  $aag->add_member($gene_id,\%type_flags);
  
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
use Bio::EnsEMBL::Utils::Scalar qw(check_ref assert_integer assert_ref);

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
    
    $self->{'MEMBERS'} = $list || [];
    
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
  my ($self, $gene_id,$type_hash) = @_;
  if(!$self->contains_member($gene_id)) {
    push(@{$self->{MEMBERS}}, [$gene_id, {}]);
    $self->set_attribs($gene_id, $type_hash);
  }
  return;
}

=head2 get_all_members_with_type

  Arg [1]     : String The type to search members by
  Description : Loops through the internal members array returning all
                attributes of the same type as what has been specified
  Example     : my $members = $aag->get_all_members_with_type('IS_VALID_ALTERNATE');

=cut

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

=head2 attribs

  Arg [1]     : Int gene id to record attributes against
  Description : Returns all known attributes of the given gene id. Attributes
                are returned as a HashRef but is a copy of the interally
                held attribute list
  Returntype  : HashRef copy of all the given id's attributes
  Example     : $aag->attribs(10, 'IS_VALID_ALTERNATE');
                $aag->attribs(10, [ 'IS_VALID_ALTERNATE' ]);
                $aag->attribs(10, {IS_VALID_ALTERNATE => 1});

=cut

sub attribs {
  my ($self, $gene_id) = @_;
  assert_integer($gene_id, 'gene_id');
  foreach my $member (@{$self->{MEMBERS}}) {
    if($member->[0] == $gene_id) {
      my $attribs = $member->[1];
      return { %{$attribs} }; # make a copy. never leak
    }
  }
  return {};
}

=head2 set_attribs

  Arg [1]     : Int gene id to set attributes against
  Arg [2]     : ArrayRef/HashRef/Scalar The attribute you wish to record
  Description : Adds the given type to the specified gene id in this group. You
                can specify the type using an ArrayRef, HashRef or a single scalar
  Example     : $aag->attribs(10, 'IS_VALID_ALTERNATE');
                $aag->attribs(10, [ 'IS_VALID_ALTERNATE' ]);
                $aag->attribs(10, {IS_VALID_ALTERNATE => 1});

=cut

sub set_attribs {
  my ($self, $gene_id, $attribs) = @_;
  assert_integer($gene_id, 'gene_id');
  my $current_attribs = $self->attribs($gene_id);
  if(check_ref($attribs, 'ARRAY')) {
    #Loop & add
    $current_attribs->{uc($_)} = 1 for @{$attribs};
  }
  elsif(check_ref($attribs, 'HASH')) {
    #loop through the keys adding them in
    foreach my $key (keys %{$attribs}) {
      $current_attribs->{uc($key)} = 1;
    }
  }
  #Simple scalar value so just add it in
  else {
    $current_attribs->{uc($attribs)} = 1;
  }
  foreach my $member (@{$self->{MEMBERS}}) {
    if($member->[0] == $gene_id) {
      $member->[1] = $current_attribs;
    }
  }
  return;
}

=head2 remove_attribs

  Arg [1]     : Int gene id to retrieve attributes against
  Arg [2]     : ArrayRef/HashRef/Scalar The attribute you wish to remove
  Description : Removes the given type from this group against the specified
                gene identifier
  Example     : $aag->remove_attribs(10, 'IS_VALID_ALTERNATE');
                $aag->remove_attribs(10, [ 'IS_VALID_ALTERNATE' ]);
                $aag->remove_attribs(10, {IS_VALID_ALTERNATE => 1});

=cut

sub remove_attribs {
  my ($self, $gene_id, $attribs) = @_;
  assert_integer($gene_id, 'gene_id');
  my @to_remove;
  if(check_ref($attribs, 'ARRAY')) {
    @to_remove = map { uc($_) } @{$attribs};
  }
  elsif(check_ref($attribs, 'HASH')) {
    @to_remove = map { uc($_) } keys %{$attribs};
  }
  #Simple scalar value so just add it in
  else {
    @to_remove = uc($attribs);
  }
  foreach my $member (@{$self->{MEMBERS}}) {
    if($member->[0] == $gene_id) {
      my $current_attribs = $member->[1];
      delete $current_attribs->{$_} for @to_remove;
    }
  }
  return;
}

=head2 remove_member

  Arg [1]     : Int gene id to retrieve attributes against
  Arg [2]     : ArrayRef/HashRef/Scalar The attribute you wish to remove
  Description : Removes the given member from this group. Any changes
                must be persisted back to the database via update() or
                store() methods in Bio::EnsEMBL::DBSQL::AltAlleleGroupAdaptor.
  Example     : $aag->remove_member(10);

=cut

sub remove_member {
  my ($self, $gene_id) = @_;
  assert_integer($gene_id, 'gene_id');
  my $members = $self->{MEMBERS};
  my $size = scalar(@{$members});
  for(my $i = 0; $i < $size; $i++) {
    my $current_id = $members->[$i]->[0];
    #If this was the ID then splice it out of the array and exit
    if($current_id == $gene_id) {
      splice(@{$members}, $i, 1);
      last;
    }
  }
  return;
}

=head2 contains_member

  Arg [1]     : Int gene id to retrieve attributes against
  Description : Searches through the members list looking for the
                specified gene id. Returns true if it was found
                or false if not.
  Returntype  : Boolean indicating if the given gene id is held in this group
  Example     : $aag->contains_member(10);

=cut

sub contains_member {
  my ($self, $gene_id) = @_;
  assert_integer($gene_id, 'gene_id');
  foreach my $member (@{$self->{MEMBERS}}) {
    if($member->[0] == $gene_id) {
      return 1;
    }
  }
  return 0;
}

=head2 remove_all_members

  Description : Remove members from this object, but NOT the database. See
                AltAlleleGroupAdaptor->remove() to remove the group from the
                database
=cut

sub remove_all_members {
  my $self = shift;
  $self->{'MEMBERS'} = [];
  return;
}

=head2 rep_Gene_id

  Arg[1]      : Optional - set a new representative Gene id for the group
  Description : Reports or sets the representative Gene for this AltAlleleGroup
                If you wish to remove the representative status of all genes without
                setting a new one, see unset_rep_Gene_id
  Returntype  : Integer or undef if none set
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

  Description : Removes the representative Gene flag from this AltAlleleGroup.
                This action is not possible through rep_Gene_id due to
                validation of inputs.

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
  Arg[2]      : ArrayRef - Can contain dbIDs or Gene objects to exclude from the returned list
  Description : fetches all the Gene dbIDs within the allele group. It can also
                be used to list those ids that are not the representative Gene.
  Returntype  : ArrayRef of gene dbIDs

=cut

sub get_all_Gene_ids {
    my $self = shift;
    my $all_but_rep = shift;
    my $excluded_genes = shift;
    my $list = $self->{'MEMBERS'};

    my %gene_exclusions;
    if($excluded_genes) {
      assert_ref($excluded_genes, 'ARRAY', 'excluded genes');
      foreach my $gene (@{$excluded_genes}) {
        my $gene_id = (ref($gene)) ? $gene->dbID() : $gene;
        $gene_exclusions{$gene_id} = $gene_id;
      }
    }
    
    my @gene_ids;
    
    foreach my $allele (@$list) {
        my ($gene_id,$type) = @$allele;
        if ($all_but_rep && $type->{IS_REPRESENTATIVE}) {next;} 
        if(exists $gene_exclusions{$gene_id}) { next; }
        push @gene_ids,$gene_id;
    }
    return [sort {$a <=> $b} @gene_ids];
}

=head2 get_representative_Gene

  Description : Used to fetch a Gene object which has been marked as the
                representative Gene for this alt allele group.
  Returntype  : Bio::EnsEMBL::Gene object which is the representative gene

=cut


sub get_representative_Gene {
    my $self = shift;
    my $ga = $self->adaptor->db->get_GeneAdaptor;
    
    return $ga->fetch_by_dbID($self->rep_Gene_id);
}

=head2 get_all_Genes

  Arg[1]      : Boolean - Do not include representative gene in list of ids.
  Arg[2]      : ArrayRef - Can contain dbIDs or Gene objects to exclude from the returned list
  Description : Fetches all the Gene objects within the allele group. It can also
                be used to list those Genes that are not the representative Gene.
  Returntype  : ArrayRef of Bio::EnsEMBL::Gene objects

=cut


sub get_all_Genes {
  my ($self, $all_but_rep, $excluded_genes) = @_;
  my $gene_ids = $self->get_all_Gene_ids($all_but_rep, $excluded_genes);
  return $self->adaptor()->db()->get_GeneAdaptor()->fetch_all_by_dbID_list($gene_ids);
}


=head2 get_all_Genes_types

  Arg[1]      : Boolean - Do not include representative gene in list of ids.
  Arg[2]      : ArrayRef - Can contain dbIDs or Gene objects to exclude from the returned list
  Description : Fetches all the Gene objects within the allele group and their
                associcated attributes. It can also be used to list those 
                Genes that are not the representative Gene.
  Returntype  : ArrayRef. 2 dimensional holding [Bio::EnsEMBL::Gene, {attribute_hash}]

=cut

sub get_all_Genes_types {
  my ($self, $all_but_rep, $excluded_genes) = @_;
  my $gene_ids = $self->get_all_Gene_ids($all_but_rep, $excluded_genes);
  my $ga = $self->adaptor()->db()->get_GeneAdaptor();
  my @output;
  my $members = $self->{MEMBERS};
  foreach my $allele (@{$members}) {
    my ($gene_id,$attribs) = @$allele;
    if ($all_but_rep && $attribs->{IS_REPRESENTATIVE}) {next;} 
    my $gene = $ga->fetch_by_dbID($gene_id);
    my %attribs_copy = %{$attribs};
    push(@output, [$gene, \%attribs_copy])
  }
  return \@output;
}

=head2 size

  Description : Returns the current size of this group in members
  Returntype  : Int the size of the current alt allele group

=cut

sub size {
  my $self = shift;
  my $list = $self->{'MEMBERS'};
  return scalar(@$list);
}


=head2 get_all_members

  Description : Retrieves all of the information about all members. Be aware
                that this emits the interal data structure so direct modification
                should be done with caution.
  Returntype  : ArrayRef of id and type list: [gene_id,type]
  Caller      : AltAlleleGroupAdaptor->store

=cut

sub get_all_members {
  my $self = shift;
  my $members = $self->{'MEMBERS'};
  return $members;
}


1;
