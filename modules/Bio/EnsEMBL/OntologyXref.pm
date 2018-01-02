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

Bio::EnsEMBL::OntologyXref

=head1 DESCRIPTION

This class extends the DBEntry in order to associate Evidence Tags
to the relationship between EnsEMBL objects and ontology accessions
(primarily GO accessions).

The relationship to GO that is stored in the database is actually
derived through the relationship of EnsEMBL peptides to SwissProt
peptides, i.e. the relationship is derived like this:

  ENSP -> SWISSPROT -> GO

And the evidence tag describes the relationship between the SwissProt
Peptide and the GO entry.

In reality, however, we store this in the database like this:

  ENSP -> SWISSPROT
  ENSP -> GO

and the evidence tag hangs off of the relationship between the ENSP and
the GO identifier.  Some ENSPs are associated with multiple closely
related Swissprot entries which may both be associated with the same GO
identifier but with different evidence tags.  For this reason a single
'OntologyXref' can have multiple evidence tags.

=head1 SYNOPSIS

  my $ontology_xref = Bio::EnsEMBL::OntologyXref->new();
  $ontology_xref->add_linkage_type('IEA');

  foreach my $evtag ( @{ $ontology_xref->get_all_linkage_types() } ) {
    print "$evtag\n";
  }

=head1 METHODS

=cut

package Bio::EnsEMBL::OntologyXref;

use strict;
use warnings;

use base qw( Bio::EnsEMBL::DBEntry );

require Bio::EnsEMBL::Registry;

=head2 get_OntologyTerm

  Example    : $ontology_xref->get_OntologyTerm();
  Description: Queries the OntologyTermAdaptor for a term which is the same 
               as the primary id of this object. This method requires a
               OntologyDBAdaptor to be available in the Bio::EnsEMBL::Registry. 
               If you have loaded data from an Ensembl release using
               Bio::EnsEMBL::Registry->load_registry_from_db() then this should
               work.
  Returntype : Bio::EnsEMBL::OntologyTerm
  Exceptions : None
  Caller     : general
  Status     : Experimantal

=cut


sub get_OntologyTerm {
  my ($self) = @_;
  my $dbas = Bio::EnsEMBL::Registry->get_all_DBAdaptors(-GROUP => 'ontology');
  foreach my $ontology_dba (@{$dbas}) {
    my $ota = $ontology_dba->get_OntologyTermAdaptor();
    my $term = $ota->fetch_by_accession($self->primary_id());
    return $term if $term;
  }
  return;
}

=head2 add_linkage_type

  Arg [1]    : string $value
               allowed values:
               'IC', 'IDA', 'IEA', 'IEP', 'IGI', 'IMP', 'IPI',
               'ISS', NAS', 'ND', 'TAS', 'NR', 'RCA'
  Arg [2]    : (optional) Bio::EnsEMBL::DBEntry $source
  Example    : $ontology_xref->add_linkage_type('IGI');
  Description: Associates a linkage type and source DBEntry with
               this ontology_xref
  Returntype : integer; number of linkages
  Exceptions : thrown if $linkage_type argument not supplied or
               the optional DBEntry is not a DBEntry object.
  Caller     : DBEntryAdaptor
  Status     : Experimantal

=cut

sub add_linkage_type {
  my ( $self, $lt, $source_dbentry ) = @_;

  if ( !defined($lt) ) {
    $self->throw("linkage type argument required");
  }

  if ( defined($source_dbentry)
       && !$source_dbentry->isa('Bio::EnsEMBL::DBEntry') )
  {
    $self->throw("source_dbentry must be a Bio::EnsEMBL::DBEntry");
  }

  $self->{'linkage_types'} ||= [];

  push @{ $self->{'linkage_types'} },
    [ $lt, ( $source_dbentry || () ) ];
}


=head2 get_all_linkage_info

  Arg [1]    : none
  Example    :

    foreach ( @{ $ontology_xref->get_all_linkage_info() } ) {
      print "evidence: $_->[0] via $_->[1]->display_id";
    }

  Description: Retrieves a list of evidence-tag/source-DBEntry pairs
               associated with this ontology_xref
  Returntype : listref of listrefs
  Exceptions : none
  Caller     : geneview? general.
  Status     : Experimental

=cut

sub get_all_linkage_info {
  my ($self) = @_;

  return $self->{'linkage_types'} || [];
}


=head2 get_all_linkage_types

  Arg [1]    : none
  Example    :

    print( join( ' ', @{ $ontology_xref->get_all_linkage_types() } ),
           "\n" );

  Description: Retrieves a unique list of evidence tags associated with
               this ontology_xref
  Returntype : none
  Exceptions : none
  Caller     : geneview? general
  Status     : Stable

=cut

sub get_all_linkage_types {
  my ($self) = @_;

  my %seen;
  return [ grep  { !$seen{$_}++ }
             map { $_->[0] } @{ $self->{'linkage_types'} } ];

  #return [ map{ $_->[0]} @{ $self->{'linkage_types'} || [] } ];
}


=head2 flush_linkage_types

  Arg [1]    : none
  Example    : $ontology_xref->flush_linkage_types();
  Description: Removes any associated evidence tags
  Returntype : none
  Exceptions : none
  Caller     : general 
  Status     : Stable

=cut

sub flush_linkage_types {
  my ($self) = @_;

  $self->{'linkage_types'} = [];
}


=head2 add_associated_xref

  Arg [1]    : Bio::EnsEMBL::DBEntry $associated_xref
               or an Array of Bio::EnsEMBL::DBEntry for compound annotations
  Arg [2]    : Bio::EnsEMBL::DBEntry $source_dbentry
  Arg [3]    : string $condition_type or an Array of string $condition_types
               matching the order of Arg[1] for compound queries.
  Arg [4]    : (optional) Integer $group id for compound annotations.
  Arg [5]    : (optional) Integer $rank order for a term within a compound annotation.
  Example    : $ontology_xref->add_associated_xref(
                                  $associated_xref,
                                  $source_dbentry,
                                  'with',
                                  42,
                                  5);
  Description: Associates a linkage type and source DBEntry with
               this ontology_xref
  Returntype : none
  Exceptions : thrown if $linkage_type argument not supplied or
               the optional DBEntry is not a DBEntry object.
  Caller     : DBEntryAdaptor
  Status     : Experimantal

=cut

sub add_associated_xref {
  my ( $self, $associated_xref, $source_dbentry, $condition_type, $group, $rank ) = @_;
  if ( ref($associated_xref) eq 'ARRAY' )
  {
    foreach my $e ( @{ $associated_xref } ) {
      if ( !$e->isa('Bio::EnsEMBL::DBEntry') ) {
        $self->throw("associated_xref must be a Bio::EnsEMBL::DBEntry or an 
                  Array of Bio::EnsEMBL::DBEntry objects.");
      }
    }
  } elsif ( defined($associated_xref)
       && !$associated_xref->isa('Bio::EnsEMBL::DBEntry') )
  {
    $self->throw("associated_xref must be a Bio::EnsEMBL::DBEntry or an Array
                  of Bio::EnsEMBL::DBEntry objects.");
  }
  
  if ( defined($source_dbentry)
       && !$source_dbentry->isa('Bio::EnsEMBL::DBEntry') )
  {
    $self->throw("source_dbentry must be a Bio::EnsEMBL::DBEntry");
  }
  
  if ( !defined $condition_type ) {
      $self->throw("condition must be a string");
  }
  
  if (!defined $group) {
    $group = 1 + scalar keys %{ $self->{'associated_xref'} };
  }
  
  if (!defined $rank) {
    #$rank = 0;
    $rank = 1 + scalar keys %{ $self->{'associated_xref'}->{ $group } };
  }

  $self->{'associated_xref'} ||= {};

  $self->{'associated_xref'}->{ $group }->{$rank} = 
    [ $associated_xref, $source_dbentry, $condition_type ];
}


=head2 add_linked_associated_xref

  Arg [1]    : Bio::EnsEMBL::DBEntry $associated_xref
               or an Array of Bio::EnsEMBL::DBEntry for compound annotations
  Arg [2]    : Bio::EnsEMBL::DBEntry $source_dbentry
  Arg [3]    : string $condition_type or an Array of string $condition_types
               matching the order of Arg[1] for compound queries.
  Arg [4]    : Integer $group id.
  Arg [5]    : Integer $rank id.
  Example    : $ontology_xref->add_associated_xref(
                                  $associated_xref,
                                  $source_dbentry,
                                  'with',
                                  42,
                                  5);
  Description: Associates a linkage type and source DBEntry with this
               ontology_xref that have come from the same annotation source
  Returntype : none
  Exceptions : thrown if $linkage_type argument not supplied or
               the optional DBEntry is not a DBEntry object.
  Caller     : DBEntryAdaptor
  Status     : Experimantal

=cut

sub add_linked_associated_xref {
  my ( $self, $associated_xref, $source_dbentry, $condition_type, $associate_group_id, $associate_group_rank ) = @_;
  
  if ( !defined($associated_xref) )
  {
    $self->throw("associated_xref must be a Bio::EnsEMBL::DBEntry or an Array
                  of Bio::EnsEMBL::DBEntry objects.");
  }
  
  if ( defined($associated_xref)
       && !$associated_xref->isa('Bio::EnsEMBL::DBEntry') )
  {
    $self->throw("associated_xref must be a Bio::EnsEMBL::DBEntry or an Array
                  of Bio::EnsEMBL::DBEntry objects.");
  }
  
  if ( defined($source_dbentry)
       && !$source_dbentry->isa('Bio::EnsEMBL::DBEntry') )
  {
    $self->throw("source_dbentry must be a Bio::EnsEMBL::DBEntry");
  }
  
  if ( !defined $condition_type ) {
      $self->throw("condition must be a string");
  }
  
  if ( !defined $associate_group_id )
  {
    $self->throw("$associate_group_id must be an integer");
  }
  
  if ( !defined $associate_group_rank )
  {
    $self->throw("$associate_group_rank must be an integer");
  }
  
#  print "\t" . $associate_group_id;
#  print "\t" . $associate_group_rank;
#  print "\t|" . defined($associated_xref) . '|';
#  #print "\t" . defined($associated_xref->primary_id);
#  print "\t" . $associated_xref->primary_id;# . ' (' . $associated_xref->display_id . ')';
#  print "\t" . $source_dbentry->primary_id;
#  print "\t" . $condition_type . "\n"; 
  
  my $associated_xref_array = {};
  my $matching_link = 0;
  
  my $load_postion = 0;
  if ( !defined $self->{'associated_xref'} ) {
    $associated_xref_array->{$associate_group_id}->{$associate_group_rank} = [
      $associated_xref,
      $source_dbentry,
      $condition_type
    ];
    $load_postion = 1;
  } else {
    $associated_xref_array = $self->{'associated_xref'};
    
    if (!defined $associated_xref_array->{$associate_group_id} ) {
      $associated_xref_array->{$associate_group_id}->{$associate_group_rank} = [
        $associated_xref,
        $source_dbentry,
        $condition_type
      ];
      $load_postion = 2;
    } else {
      my $already_loaded = 0;
      foreach my $rank ( keys %{$associated_xref_array->{$associate_group_id}} ) {
        my @ax_gr_set = @{ $associated_xref_array->{$associate_group_id}->{$rank} };
        if (
          $ax_gr_set[0]->primary_id eq $associated_xref->primary_id && 
          $ax_gr_set[1]->primary_id eq $source_dbentry->primary_id &&
          $ax_gr_set[2] eq $condition_type
        ) {
          $already_loaded = 1;
          $load_postion = 4;
          #last;
        }
      }
      if ( !$already_loaded ) {
        if (!defined $associated_xref_array->{$associate_group_id}->{$associate_group_rank} ) {
          $associated_xref_array->{$associate_group_id}->{$associate_group_rank} = [
            $associated_xref,
            $source_dbentry,
            $condition_type
          ];
          $load_postion = 5;
        } else {
          $associated_xref_array->{$associate_group_id}->{scalar keys %{ $associated_xref_array->{$associate_group_id} }} = [
            $associated_xref,
            $source_dbentry,
            $condition_type
          ];
          $load_postion = 3;
        }
      }
    }
  }
  $self->{'associated_xref'} = $associated_xref_array;
#  print "\t\tLoaded at " . $load_postion . "\n";
#  if ($associate_group_id == 28792 || $associate_group_id == 28793) {
#    print Data::Dumper->Dumper([$self->{'associated_xref'}]);
#  }
}


=head2 get_all_associated_xrefs

  Arg [1]    : none
  Example    :

    foreach ( @{ $ontology_xref->get_all_associated_xref() } ) {
      print "evidence: $_->[0] via $_->[1]->display_id";
    }

  Description: Retrieves a list of associated-DBEntry/source-DBEntry/condition
               sets associated with this ontology_xref
  Returntype : listref of listrefs
  Exceptions : none
  Caller     : geneview? general.
  Status     : Experimental

=cut

sub get_all_associated_xrefs {
  my ($self) = @_;

  return $self->{'associated_xref'} || {};
}

=head2 get_extensions
  Arg [1]    : none
  Example    :

    use Data::Dumper;
    print Dumper @{ $ontology_xref->get_extensions_for_web() }
    
  Returns    :
    $VAR1 = {
              'source' => '<a href="<Link to CiteXplore>">11937031</a>',
              'evidence' => 'IDA',
              'description' => '<strong>has_direct_input</strong> 
                   <a href="http://www.pombase.org/spombe/result/SPBC32F12.09">
                     SPBC32F12.09
                   </a>'
            };

  Description: Retrieves a list of associated-DBEntry/source-DBEntry/condition
               sets associated with this ontology_xref and formats them ready
               for web display in a group per row fashion.
               The accessions for ontologies are linkable. Extra links need to
               be added for each distinct database that is reference. 
  Returntype : listref of hashrefs
  Exceptions : none
  Caller     : 
  Status     : Experimental
=cut
sub get_extensions {
  my ($self) = @_;
  
  if ( !defined $self->{'associated_xref'} ) {
    return [];
  }
  
  my @annotExtRows = ();
  
  my %external_urls = (
    'GO'     => 'http://www.ebi.ac.uk/ego/GTerm?id=',
    'GO_REF' => 'http://www.geneontology.org/cgi-bin/references.cgi#',
    'SO'     => 'http://www.sequenceontology.org/miso/current_cvs/term/',
    #'MOD'    => 'href=mod',
    'PomBase' => 'http://www.pombase.org/spombe/result/',
    'PomBase_Systematic_ID' => 'http://www.pombase.org/spombe/result/',
    'PUBMED' => 'http://europepmc.org/abstract/MED/'
  );
  
  foreach my $groupId (keys %{ $self->{'associated_xref'} } ) {
    my $description = '';
    my $evidence    = '';
    my $source      = '';
    
    foreach my $rank ( keys %{ $self->{'associated_xref'}->{$groupId} } ) {
      if ( $self->{'associated_xref'}->{$groupId}->{$rank}->[2]
           eq 'evidence' ) {
        if ( exists $external_urls{$self->{'associated_xref'}->{$groupId}->{$rank}->[0]->dbname} ) {
          $evidence = '<a href="' . $external_urls{$self->{'associated_xref'}->{$groupId}->{$rank}->[0]->dbname}
                    . $self->{'associated_xref'}->{$groupId}->{$rank}->[0]->primary_id . '">'
                    . $self->{'associated_xref'}->{$groupId}->{$rank}->[0]->display_id
                    . '</a>';
        } else {
          $evidence = $self->{'associated_xref'}->{$groupId}->{$rank}->[0]->display_id;
        }
      } else {
        if (length $description > 0) {
          $description .= ', ';
        }
        $description .= '<strong>' . $self->{'associated_xref'}->{$groupId}->{$rank}->[2] . '</strong> ';
        
        if (exists $external_urls{$self->{'associated_xref'}->{$groupId}->{$rank}->[0]->dbname} ) {
          $description .= '<a href="' . $external_urls{$self->{'associated_xref'}->{$groupId}->{$rank}->[0]->dbname}
                        . $self->{'associated_xref'}->{$groupId}->{$rank}->[0]->primary_id . '">'
                        . $self->{'associated_xref'}->{$groupId}->{$rank}->[0]->display_id
                        . '</a>';
        } else {
          $description .= $self->{'associated_xref'}->{$groupId}->{$rank}->[0]->display_id;
        }
        
        
      }
      
      if ( !undef $external_urls{$self->{'associated_xref'}->{$groupId}->{$rank}->[1]} ) {
        if ( exists $external_urls{$self->{'associated_xref'}->{$groupId}->{$rank}->[1]->dbname} ) {
          $source = '<a href="' . $external_urls{$self->{'associated_xref'}->{$groupId}->{$rank}->[1]->dbname}
                  . $self->{'associated_xref'}->{$groupId}->{$rank}->[1]->primary_id . '">';
          if ( $self->{'associated_xref'}->{$groupId}->{$rank}->[1]->dbname ne 'GO_REF' ) {
            $source .= 'PMC:';
          }
          $source .= $self->{'associated_xref'}->{$groupId}->{$rank}->[1]->display_id . '</a>';
        } elsif ($self->{'associated_xref'}->{$groupId}->{$rank}->[1]->display_id eq 'PMPB:0') {
          $source = '';
        } else {
          $source = $self->{'associated_xref'}->{$groupId}->{$rank}->[1]->display_id;
        }
      }
    }
    
    if ($evidence ne '' and $description ne '') {
      my %row =  ('description' => $description,
                  'evidence'    => $evidence,
                  'source'      => $source);
      push @annotExtRows, (\%row);
    }
  }
  
  return \@annotExtRows;
}

1;
