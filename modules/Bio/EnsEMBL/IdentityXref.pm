#
#
# BioPerl module for SimilarityXref.pl
#
# Cared for by Emmanuel Mongin <mongin@ebi.ac.uk>
#
# Copyright Emmanuel Mongin
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::IdentiyXref

=head1 SYNOPSIS

my $xref = Bio::EnsEMBL::IdentityXref->new;

=head1 CONTACT

Post questions to the ensembl development list: <ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...
package Bio::EnsEMBL::IdentityXref;
use vars qw(@ISA $AUTOLOAD);
use strict;


@ISA = qw( Bio::EnsEMBL::DBEntry );

=head2 new
  
  See Bio::EnsEMBL::DBEntry::new

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    my ($query_identity, $target_identity, $score, $evalue, 
        $cigar_line, $query_start, $query_end, $translation_start,
        $translation_end, $analysis ) = $self->_rearrange(
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

=cut

sub analysis {
   my $self = shift;
  $self->{'analysis'} = shift if( @_ );
  return $self->{'analysis'};
}




=head2 get_mapper

  Args       : none
  Example    : none
  Description: produces a mapper object that takes coordinates from one side of 
               the alignment to the other side. "ensembl" and "external" are the 
               two coordinate systems contained. 
  Returntype : Bio::EnsEMBL::Mapper
  Exceptions : none
  Caller     : general, ProteinDAS subsystem

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
  my $hit_id = "external_id";
  my $translation_id = "translation_id";
  # now build the mapper
  my $mapper = Bio::EnsEMBL::Mapper->new( "ensembl", "external" );


  for( my $i=0; $i<=$#lens; $i++ ) {
    my $length = $lens[$i];
    my $char = $chars[$i];
    if( $char eq "M" ) {
      $mapper->add_map_coordinates( $translation_id, $translation_start,
				    $translation_start + $length - 1,
				    1, $hit_id, $query_start, $query_start + $length - 1 );
      $query_start += $length;
      $translation_start += $length;

    } elsif( $char eq "D" ) {
      $query_start += $length;
    } elsif( $char eq "I" ) {
      $translation_start += $length;
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
  Example    : my $ens_prot_feature_list = $ident_xref->
                         transform_feature( $swiss_prot_feature ); 
  Description: a list of potential partial features which represent all mappable places
               of the original feature in ensembl translation coordinates.
  Returntype : listref of whatever was put in
  Exceptions : none
  Caller     : general, ProteinDAS subsystem

=cut


sub transform_feature {
  my $self= shift;
  my $feature = shift;

  my $fstart = $feature->start();
  my $fend = $feature->end();
  
  my $mapper = $self->get_mapper();
  my @result;

  my @coords = $mapper->map_coordinates( "external_id", $fstart, $fend, 1, "external" );
  
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

=cut

sub map_feature {
  my $self = shift;
  my $feature = shift;


  my $fstart = $feature->start();
  my $fend = $feature->end();
  
  my $mapper = $self->get_mapper();
  my @coords = $mapper->map_coordinates( "external_id", $fstart, $fend, 1, "external" );
  
  return @coords;
}







1;
