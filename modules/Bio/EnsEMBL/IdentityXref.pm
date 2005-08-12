#
#
# EnsEMBL module for IdentityXref
#
# Copyright Emmanuel Mongin
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::IdentiyXref

=head1 SYNOPSIS

my $xref = Bio::EnsEMBL::IdentityXref->new
                               (-QUERY_IDENTITY  => 80.4,
                                -TARGET_IDENTITY => 90.1,
                                -SCORE           => 90,
                                -EVALUE          => 12,
                                -CIGAR_LINE      => '23MD3M2I40M',
                                -QUERY_START     => 1,
                                -QUERY_END       => 68,
                                -TRANSLATION_START => 10,
                                -TRANSLATION_END   => 77,
                                -ANALYSIS        => $analysis,
                                -ADAPTOR         => $adaptor,
                                -PRIMARY_ID      => $primary_id,
                                -DBNAME          => 'SwissProt');
                                
                                         

=head1 CONTACT

Post questions to the ensembl development list: <ensembl-dev@ebi.ac.uk>

=head1 METHODS

=cut

package Bio::EnsEMBL::IdentityXref;
use vars qw(@ISA $AUTOLOAD);
use strict;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

@ISA = qw( Bio::EnsEMBL::DBEntry );


=head2 new

  Arg [...]  : QUERY_IDENTITY TARGET_IDENTITY SCORE EVALUE CIGAR_LINE
             : QUERY_START QUERY_END TRANSLATION_START TRANSLATION_END
             : ANALYSIS   pairs
  Example    : see synopsis 
  Description: Create a new Bio::EnsEMBL::IdentityXref object
  Returntype : Bio::EnsEMBL::IdentityXref
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    my ($query_identity, $target_identity, $score, $evalue, 
        $cigar_line, $query_start, $query_end, $translation_start,
        $translation_end, $analysis ) = rearrange(
        [qw(QUERY_IDENTITY TARGET_IDENTITY SCORE EVALUE CIGAR_LINE 
            QUERY_START QUERY_END TRANSLATION_START TRANSLATION_END
            ANALYSIS)], @args);

    $self->{'query_identity'} = $query_identity;
    $self->{'target_identity'} = $target_identity;
    $self->{'score'} = $score;
    $self->{'evalue'} = $evalue;
    $self->{'cigar_line'} = $cigar_line;
    $self->{'query_start'} = $query_start;
    $self->{'query_end'} = $query_end;
    $self->{'translation_start'} = $translation_start;
    $self->{'translation_end'} = $translation_end;
    $self->{'analysis'} = $analysis;

    return $self;
}

=head2 query_identity

  Arg [1]    : (optional) string $value
  Example    : $query_identity = $id_xref->query_identity;
  Description: Getter/Setter for query identity
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub query_identity{
   my $obj = shift;
   if( @_ ) {
       my $value = shift;
       $obj->{'query_identity'} = $value;
   }
   return $obj->{'query_identity'};

}


=head2 target_identity

  Arg [1]    : (optional) string $value
  Example    : $target_identity = $id_xref->target_identity;
  Description: Getter/Setter for query identity
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub target_identity{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'target_identity'} = $value;
    }
    return $obj->{'target_identity'};

}



=head2 cigar_line

  Arg [1]    : string $cigar_line
  Example    : none
  Description: get/set for attribute cigar_line
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub cigar_line {
  my $self = shift;
  $self->{'cigar_line'} = shift if( @_ );
  return $self->{'cigar_line'};
}


=head2 translation_start

  Arg [1]    : string $translation_start
  Example    : none
  Description: get/set for attribute translation_start
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub translation_start {
   my $self = shift;
  $self->{'translation_start'} = shift if( @_ );
  return $self->{'translation_start'};
}


=head2 translation_end

  Arg [1]    : string $translation_end
  Example    : none
  Description: get/set for attribute translation_end
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub translation_end {
   my $self = shift;
  $self->{'translation_end'} = shift if( @_ );
  return $self->{'translation_end'};
}


=head2 query_start

  Arg [1]    : string $query_start
  Example    : none
  Description: get/set for attribute query_start
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub query_start {
   my $self = shift;
  $self->{'query_start'} = shift if( @_ );
  return $self->{'query_start'};
}


=head2 query_end

  Arg [1]    : string $query_end
  Example    : none
  Description: get/set for attribute query_end
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub query_end {
   my $self = shift;
  $self->{'query_end'} = shift if( @_ );
  return $self->{'query_end'};
}


=head2 score

  Arg [1]    : string $score
  Example    : none
  Description: get/set for attribute score
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub score {
   my $self = shift;
  $self->{'score'} = shift if( @_ );
  return $self->{'score'};
}


=head2 evalue

  Arg [1]    : string $evalue
  Example    : none
  Description: get/set for attribute evalue
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub evalue {
   my $self = shift;
  $self->{'evalue'} = shift if( @_ );
  return $self->{'evalue'};
}


=head2 analysis

  Arg [1]    : Bio::EnsEMBL::Analysis $analysis
  Example    : none
  Description: get/set for attribute analysis
  Returntype : Bio::EnsEMBL::Analysis
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub analysis {
   my $self = shift;
  $self->{'analysis'} = shift if( @_ );
  return $self->{'analysis'};
}




=head2 get_mapper

  Args       : none
  Example    : none
  Description: produces a mapper object that takes coordinates from one side 
               of the alignment to the other side. "ensembl" and "external" 
               are the two coordinate systems contained.
  Returntype : Bio::EnsEMBL::Mapper
  Exceptions : none
  Caller     : general, ProteinDAS subsystem
  Status     : Stable

=cut


sub get_mapper {
  my ( $self ) = @_;
  # parse the cigar_line and create a mapper ...
  if( exists $self->{'_cached_mapper'} ) {
    return $self->{'_cached_mapper'};
  }

  my ( @lens, @chars );

  # if there is no cigar line, nothing is going to be loaded
  if( $self->cigar_line() ) {
    my @pre_lens = split( '[DMI]', $self->cigar_line() );
    @lens = map { if( ! $_ ) { 1 } else { $_ }} @pre_lens;
    @chars = grep { /[DMI]/ } split( //, $self->cigar_line() );
  }
  my $translation_start = $self->translation_start();
  my $translation_end = $self->translation_end();
  my $query_start = $self->query_start();
  my $query_end = $self->query_end();

  #  my $hit_id = $self->display_id();
  my $ensembl_id = "ensembl_id";
  my $external_id = "external_id";
  # now build the mapper
  my $mapper = Bio::EnsEMBL::Mapper->new( "external", "ensembl" );


  for( my $i=0; $i<=$#lens; $i++ ) {
    my $length = $lens[$i];
    my $char = $chars[$i];
    if( $char eq "M" ) {
      $mapper->add_map_coordinates( $external_id, $query_start,
                                    $query_start + $length - 1, 1,
                                    $ensembl_id, $translation_start,
                                    $translation_start + $length - 1);
      $query_start += $length;
      $translation_start += $length;

    } elsif( $char eq "D" ) {
      $translation_start += $length;
    } elsif( $char eq "I" ) {
      $query_start += $length;
    }
  }

  $self->{'_cached_mapper'} = $mapper;

  return $mapper;
}



=head2 transform_feature

  Arg [1]    : a feature type with start and end $feature
               This doesnt need to be a Bio::EnsEMBL::Feature as it doesnt
               need an attached slice. We may have to introduce an appropriate
               object type.
  Example    : my $ens_prot_feature_list = 
                    $ident_xref->transform_feature( $swiss_prot_feature );
  Description: a list of potential partial features which represent all 
               mappable places
               of the original feature in ensembl translation coordinates.
  Returntype : listref of whatever was put in
  Exceptions : none
  Caller     : general, ProteinDAS subsystem
  Status     : Stable

=cut


sub transform_feature {
  my $self= shift;
  my $feature = shift;

  my $fstart = $feature->start();
  my $fend = $feature->end();

  my $mapper = $self->get_mapper();
  my @result;

  my @coords = $mapper->map_coordinates( "external_id", $fstart, $fend, 
                                           1, "external" );

  for my $coord ( @coords ) {
    if( $coord->isa( "Bio::EnsEMBL::Mapper::Coordinate" )) {
      my $new_feature;
      %{$new_feature} = %$feature;
      bless $new_feature, ref( $feature );
      $new_feature->start( $coord->start() );
      $new_feature->end( $coord->end() );

      push( @result, $new_feature );
    }
  }

  return \@result;
}



=head2 map_feature

  Arg [1]    : a start,end capable feature object
  Example    : none
  Description: 
  Returntype : list of Coordinates/Gaps which represents the mapping
  Exceptions : none
  Caller     : another way of doing ProteinDAS
  Status     : Stable

=cut

sub map_feature {
  my $self = shift;
  my $feature = shift;


  my $fstart = $feature->start();
  my $fend = $feature->end();

  my $mapper = $self->get_mapper();
  my @coords = $mapper->map_coordinates( "external_id", $fstart, $fend, 1,
                                         "external" );

  return @coords;
}



1;
