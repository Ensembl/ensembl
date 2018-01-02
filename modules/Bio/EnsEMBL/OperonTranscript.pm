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

Bio::EnsEMBL::OperonTranscript - Object representing a polycistronic transcript that is part of an operon

=head1 SYNOPSIS

my $operon_transcript =
  Bio::EnsEMBL::OperonTranscript->new( -START  => $start,
									   -END    => $end,
									   -STRAND => $strand,
									   -SLICE  => $slice );
$operon->add_OperonTranscript($operon_transcript);

=head1 DESCRIPTION

A representation of a polycistronic transcript from an operon within the Ensembl system. An operon is a collection of one or more polycistronic transcripts, which contain one or more genes.

=head1 METHODS

=cut

package Bio::EnsEMBL::OperonTranscript;

use strict;
use warnings;

use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning deprecate);
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Feature);

=head2 new

  Arg [-START]  : 
       int - start postion of the operon transcript
  Arg [-END]    : 
       int - end position of the operon transcript
  Arg [-STRAND] : 
       int - 1,-1 tehe strand the operon transcript is on
  Arg [-SLICE]  : 
       Bio::EnsEMBL::Slice - the slice the operon transcript is on
  Arg [-STABLE_ID] :
        string - the stable identifier of this operon transcript
  Arg [-VERSION] :
        int - the version of the stable identifier of this operon transcript
  Arg [-CREATED_DATE]:
        string - the date the operon transcript was created
  Arg [-MODIFIED_DATE]:
        string - the date the operon transcript was last modified

  Example    : $gene = Bio::EnsEMBL::OperonTranscript->new(...);
  Description: Creates a new operon transcript object
  Returntype : Bio::EnsEMBL::OperonTranscript
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub new {
	my $caller = shift;

	my $class = ref($caller) || $caller;
	my $self = $class->SUPER::new(@_);
	my ( $stable_id, $version, $created_date, $modified_date, $display_label ) =
	  rearrange( [  'STABLE_ID',    'VERSION',
					'CREATED_DATE', 'MODIFIED_DATE', 'DISPLAY_LABEL' ],
				 @_ );

	$self->stable_id($stable_id);
	$self->version($version);
	$self->{'created_date'}  = $created_date;
	$self->{'modified_date'} = $modified_date;
	$self->display_label($display_label);
	return $self;
}

=head2 created_date

  Arg [1]    : (optional) String - created date to set (as a UNIX time int)
  Example    : $gene->created_date('1141948800');
  Description: Getter/setter for attribute created_date
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub created_date {
  my $self = shift;
  $self->{'created_date'} = shift if ( @_ );
  return $self->{'created_date'};
}


=head2 modified_date

  Arg [1]    : (optional) String - modified date to set (as a UNIX time int)
  Example    : $gene->modified_date('1141948800');
  Description: Getter/setter for attribute modified_date
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub modified_date {
  my $self = shift;
  $self->{'modified_date'} = shift if ( @_ );
  return $self->{'modified_date'};
}


=head2 display_label

  Arg [1]    : (optional) String - the name/label to set
  Example    : $operon->name('accBCD');
  Description: Getter/setter for attribute name.
  Returntype : String or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub display_label {
	my $self = shift;
	$self->{'display_label'} = shift if (@_);
	return $self->{'display_label'};
}

=head2 stable_id

  Arg [1]    : (optional) String - the stable ID to set
  Example    : $operon->stable_id("accR2A");
  Description: Getter/setter for stable id for this operon transcript.
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub stable_id {
	my $self = shift;
	$self->{'stable_id'} = shift if (@_);
	return $self->{'stable_id'};
}

=head2 version

  Arg [1]    : (optional) Int - the stable ID version to set
  Example    : $operon->version(1);
  Description: Getter/setter for stable id version for this operon transcript.
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
sub version {
  my $self = shift;
  $self->{'version'} = shift if(@_);
  return $self->{'version'};
}

=head2 stable_id_version

  Arg [1]    : (optional) String - the stable ID with version to set
  Example    : $operon->stable_id("accR2A.3");
  Description: Getter/setter for stable id with version.
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub stable_id_version {
    my $self = shift;
    if(my $stable_id = shift) {
	# See if there's an embedded period, assume that's a
	# version, might not work for some species but you
	# should use ->stable_id() and version() if you're worried
	# about ambiguity
	my $vindex = rindex($stable_id, '.');
	# Set the stable_id and version pair depending on if
	# we found a version delimiter in the stable_id
	($self->{stable_id}, $self->{version}) = ($vindex > 0 ?
						  (substr($stable_id,0,$vindex), substr($stable_id,$vindex+1)) :
						  $stable_id, undef);
    }
    return $self->{stable_id} . ($self->{version} ? ".$self->{version}" : '');
}

=head2 operon

  Example    : $operon = $ot->operon();
  Description: getter for the operon to which this transcript belongs
  Returntype : Bio::EnsEMBL::Operon
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
sub operon {
	my $self = shift;
	if ( !exists $self->{'operon'} ) {
		if ( defined $self->adaptor() ) {
			my $ta = $self->adaptor()->db()->get_OperonAdaptor();
			my $operon = $ta->fetch_by_operon_transcript($self);
			$self->{'operon'} = $operon;
		}
	}
	return $self->{'operon'};
}

=head2 get_all_Genes

  Example    : $genes = $ot->get_all_Genes();
  Description: get all the genes that are attached to this operon transcript
  Returntype : Arrayref of Bio::EnsEMBL::Gene
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
sub get_all_Genes {
	my $self = shift;
	if(! defined $self->{_gene_array}) {
	  if(defined $self->dbID() && defined $self->adaptor()) {
	    my $ta = $self->adaptor()->db()->get_OperonTranscriptAdaptor();
      my $transcripts = $ta->fetch_genes_by_operon_transcript($self);
      $self->{_gene_array} = $transcripts;
	  }
	  else {
	    $self->{_gene_array} = [];
	  }
	}
	return $self->{_gene_array};
}
=head2 add_Gene

  Arg [1]    : Bio::EnsEMBL::Gene - gene to attach to this polycistronic transcript
  Example    : $operon->add_gene($gene);
  Description: Attach a gene to this polycistronic transcript
  Exceptions : if argument is not Bio::EnsEMBL::Gene
  Caller     : general
  Status     : Stable

=cut
sub add_Gene {
	my ($self,$gene) = @_;
	assert_ref($gene,'Bio::EnsEMBL::Gene');
	push @{$self->get_all_Genes()},$gene;
	return;
}

=head2 add_DBEntry

  Arg [1]    : Bio::EnsEMBL::DBEntry $dbe
               The dbEntry to be added
  Example    : my $dbe = Bio::EnsEMBL::DBEntery->new(...);
               $operon->add_DBEntry($dbe);
  Description: Associates a DBEntry with this operon. Note that adding DBEntries
               will prevent future lazy-loading of DBEntries for this operon
               (see get_all_DBEntries).
  Returntype : none
  Exceptions : thrown on incorrect argument type
  Caller     : general
  Status     : Stable

=cut

sub add_DBEntry {
  my $self = shift;
  my $dbe = shift;

  unless($dbe && ref($dbe) && $dbe->isa('Bio::EnsEMBL::DBEntry')) {
    throw('Expected DBEntry argument');
  }

  $self->{'dbentries'} ||= [];
  push @{$self->{'dbentries'}}, $dbe;
}


=head2 get_all_Attributes

  Arg [1]    : (optional) String $attrib_code
               The code of the attribute type to retrieve values for
  Example    : my ($author) = @{ $ot->get_all_Attributes('author') };
               my @ot_attributes = @{ $ot->get_all_Attributes };
  Description: Gets a list of Attributes of this operon transcript.
               Optionally just get Attributes for given code.
  Returntype : Listref of Bio::EnsEMBL::Attribute
  Exceptions : warning if gene does not have attached adaptor and attempts lazy
               load.
  Caller     : general
  Status     : Stable

=cut

sub get_all_Attributes {
	my $self        = shift;
	my $attrib_code = shift;

	if ( !exists $self->{'attributes'} ) {
		if ( !$self->adaptor() ) {
			return [];
		}

		my $attribute_adaptor = $self->adaptor->db->get_AttributeAdaptor();
		$self->{'attributes'} = $attribute_adaptor->fetch_all_by_OperonTranscript($self);
	}

	if ( defined $attrib_code ) {
		my @results =
		  grep { uc( $_->code() ) eq uc($attrib_code) }
		  @{ $self->{'attributes'} };
		return \@results;
	} else {
		return $self->{'attributes'};
	}
}

=head2 get_all_DBEntries

  Arg [1]    : (optional) String, external database name

  Arg [2]    : (optional) String, external_db type

  Example    : @dbentries = @{ $gene->get_all_DBEntries() };

  Description: Retrieves DBEntries (xrefs) for this operon transcript.  This does
               *not* include DBEntries that are associated with the
               transcripts and corresponding translations of this
               gene (see get_all_DBLinks()).

               This method will attempt to lazy-load DBEntries
               from a database if an adaptor is available and no
               DBEntries are present on the gene (i.e. they have not
               already been added or loaded).

  Return type: Listref of Bio::EnsEMBL::DBEntry objects
  Exceptions : none
  Caller     : get_all_DBLinks, OperontTranscriptAdaptor::store
  Status     : Stable

=cut

sub get_all_DBEntries {
	my ( $self, $db_name_exp, $ex_db_type ) = @_;

	my $cache_name = 'dbentries';

	if ( defined($db_name_exp) ) {
		$cache_name .= $db_name_exp;
	}

	if ( defined($ex_db_type) ) {
		$cache_name .= $ex_db_type;
	}

	# if not cached, retrieve all of the xrefs for this gene
	if ( !defined( $self->{$cache_name} ) && defined( $self->adaptor() ) ) {
		$self->{$cache_name} =
		  $self->adaptor()->db()->get_DBEntryAdaptor()
		  ->fetch_all_by_Operon( $self->operon(), $db_name_exp, $ex_db_type );
	}

	$self->{$cache_name} ||= [];

	return $self->{$cache_name};
}

1;

