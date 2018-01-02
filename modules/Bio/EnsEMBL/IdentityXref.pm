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

Bio::EnsEMBL::IdentityXref

=head1 SYNOPSIS

  my $xref = Bio::EnsEMBL::IdentityXref->new(
    -XREF_IDENTITY    => 80.4,
    -ENSEMBL_IDENTITY => 90.1,
    -SCORE            => 90,
    -EVALUE           => 12,
    -CIGAR_LINE       => '23MD3M2I40M',
    -XREF_START       => 1,
    -XREF_END         => 68,
    -ENSEMBL_START    => 10,
    -ENSEMBL_END      => 77,
    -ADAPTOR          => $adaptor,
    -PRIMARY_ID       => $primary_id,
    -DBNAME           => 'SwissProt'
  );

=head1 METHODS

=cut

package Bio::EnsEMBL::IdentityXref;
use vars qw(@ISA $AUTOLOAD);
use strict;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( deprecate );

@ISA = qw( Bio::EnsEMBL::DBEntry );


=head2 new

  Arg [...]  : XREF_IDENTITY ENSEMBL_IDENTITY SCORE EVALUE CIGAR_LINE
             : XREF_START XREF_END ENSEMBL_START ENSEMBL_END
             : ANALYSIS pairs
  Example    : see synopsis 
  Description: Create a new Bio::EnsEMBL::IdentityXref object
  Returntype : Bio::EnsEMBL::IdentityXref
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

my $error_shown = 0;

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    my ($query_identity, $target_identity, $score, $evalue, 
        $cigar_line, $query_start, $query_end, $translation_start,
        $translation_end, $analysis, $xref_identity, $ensembl_identity, 
        $xref_start, $xref_end, $ensembl_start, $ensembl_end) = rearrange(
        [qw(QUERY_IDENTITY TARGET_IDENTITY SCORE EVALUE CIGAR_LINE 
            QUERY_START QUERY_END TRANSLATION_START TRANSLATION_END
            ANALYSIS XREF_IDENTITY ENSEMBL_IDENTITY XREF_START XREF_END ENSEMBL_START ENSEMBL_END)], @args);

    if((defined($query_identity) or defined($target_identity) or defined($query_start) or defined ($query_end) or
       defined($translation_start) or defined($translation_end)) and !$error_shown){
      print STDERR "Arguments have now been changed to stop confusion so please replace the following\n";
      print STDERR "\tQUERY_IDENTITY\t->\tXREF_IDENTITY\n";
      print STDERR "\tTARGET_IDENTITY\t->\tENSEMBL_IDENTITY\n";
      print STDERR "\tQUERY_START\t->\tXREF_START\n";
      print STDERR "\tQUERY_END\t->\tXREF_END\n";
      print STDERR "\tTRANSLATION_START\t->\tENSEMBL_START\n";
      print STDERR "\tTRANSLATION_END\t->\tENSEMBL_END\n";
      print STDERR "The old arguments will be removed in a futute release so please change your code to the new names\n";
      $error_shown = 1;
    }
    $self->{'xref_identity'} = $query_identity || $xref_identity;
    $self->{'ensembl_identity'} = $target_identity || $ensembl_identity;
    $self->{'score'} = $score;
    $self->{'evalue'} = $evalue;
    $self->{'cigar_line'} = $cigar_line;
    $self->{'xref_start'} = $query_start || $xref_start;
    $self->{'xref_end'} = $query_end || $xref_end;
    $self->{'ensembl_start'} = $translation_start || $ensembl_start;
    $self->{'ensembl_end'} = $translation_end || $ensembl_end;
    $self->{'analysis'} = $analysis;

    return $self;
}

=head2 xref_identity

  Arg [1]    : (optional) string $value
  Example    : $xref_identity = $id_xref->xref_identity;
  Description: Getter/Setter for xref identity
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub xref_identity{
   my $obj = shift;
   if( @_ ) {
       my $value = shift;
       $obj->{'xref_identity'} = $value;
   }
   return $obj->{'xref_identity'};

}


=head2 ensembl_identity

  Arg [1]    : (optional) string $value
  Example    : $ensembl_identity = $id_xref->ensembl_identity;
  Description: Getter/Setter for ensembl identity
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub ensembl_identity{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'ensembl_identity'} = $value;
    }
    return $obj->{'ensembl_identity'};

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


=head2 ensembl_start

  Arg [1]    : string $ensembl_start
  Example    : none
  Description: get/set for attribute ensembl_start
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub ensembl_start {
   my $self = shift;
  $self->{'ensembl_start'} = shift if( @_ );
  return $self->{'ensembl_start'};
}


=head2 ensembl_end

  Arg [1]    : string $ensembl_end
  Example    : none
  Description: get/set for attribute ensembl_end
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub ensembl_end {
   my $self = shift;
  $self->{'ensembl_end'} = shift if( @_ );
  return $self->{'ensembl_end'};
}


=head2 xref_start

  Arg [1]    : string $xref_start
  Example    : none
  Description: get/set for attribute xref_start
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub xref_start {
   my $self = shift;
  $self->{'xref_start'} = shift if( @_ );
  return $self->{'xref_start'};
}


=head2 xref_end

  Arg [1]    : string $xref_end
  Example    : none
  Description: get/set for attribute xref_end
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub xref_end {
   my $self = shift;
  $self->{'xref_end'} = shift if( @_ );
  return $self->{'xref_end'};
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
  my $translation_start = $self->ensembl_start();
  my $translation_end = $self->ensembl_end();
  my $query_start = $self->xref_start();
  my $query_end = $self->xref_end();

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
