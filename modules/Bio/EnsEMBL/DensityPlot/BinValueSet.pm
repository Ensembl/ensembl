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

Bio::EnsEMBL::DensityPlot::BinValueSet

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 METHODS

=cut

package Bio::EnsEMBL::DensityPlot::BinValueSet;

use vars qw($AUTOLOAD @ISA);
use strict;
use Bio::EnsEMBL::DensityPlot::BinValue;

# Object preamble - inheriets from Bio::Root::Object

use Bio::EnsEMBL::Root;

@ISA = qw(Bio::EnsEMBL::Root Exporter);
#@EXPORT_OK = qw();
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub new {
  my ($class,@args) = @_;
    
    my $self = {};
    bless $self,$class;
    $self->{'_bin_array'} = [];
    return $self;
}



=head2 add_binvalue

 Title   : add_binValue
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_binvalue{
    my ($self,$value) = @_;

    defined ($value->chromosomestart)   || $self->throw( "Bin Value object does not contain a ChromosomeStart method" );
    defined ($value->chromosomeend)     || $self->throw( "Bin Value object does not contain a ChromosomeEnd method"   );
    defined ($value->value)             || $self->throw( "Bin Value object does not contain a Value method"           );
    $self->_store_biggest($value->value);
    $self->_store_smallest($value->value);

    push(@{$self->{'_bin_array'}},$value);
}

=head2 get_binvalues

 Title   : get_binvalues
 Usage   : my @binvalue_objects = $BVSet->get_binvalues
 Function: scales all the binvalues by the scale_factor and returns them.
 Example :
 Returns : array of BinValue objects 
 Args    : none


=cut

sub get_binvalues{
    my $self = shift;
    my $biggest_value   = $self->{'_biggest_value'} || 0;
    my $smallest_value  = $self->{'_smallest_value'} || 0;
   
    if (!defined ($biggest_value)||!defined($smallest_value)){
    	$self->throw("Cannot scale - no values to scale against");  
    }

    my $width = $self->scale_to_fit();

    if ($self->stretch && ($biggest_value-$smallest_value) ){
    	foreach my $bv (@{ $self->{'_bin_array'}}){
	        my $scaledval = (($bv->value - $smallest_value) / ($biggest_value-$smallest_value) )* $width;
	        $bv->scaledvalue($scaledval);
    	}
    } elsif($biggest_value) {
        foreach my $bv (@{ $self->{'_bin_array'}}){
	        my $scaledval = ($bv->value / $biggest_value) * $width;
	        $bv->scaledvalue($scaledval);
	    }
    } else {
        foreach my $bv (@{ $self->{'_bin_array'}}){
	        $bv->scaledvalue(0);
	    }
    }

   return ( @{ $self->{'_bin_array'}}      );  

}

sub size {
    my $self = shift;
    return scalar @{$self->{'_bin_array'}};
}

=head2 position

 Title   : position
 Usage   : $obj->position($newval)
 Function: 
 Returns : value of position
 Args    : newvalue (optional)


=cut

sub position{
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'position'} = $value;
    }
    return $self->{'position'};

}


=head2 label

 Title   : label
 Usage   : $obj->label($newval)
 Function: 
 Returns : value of label
 Args    : newvalue (optional)


=cut

sub label{
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'label'} = $value;
    }
    return $self->{'label'};

}


=head2 label2

 Title   : label2
 Usage   : $obj->label2($newval)
 Function: 
 Returns : value of label2
 Args    : newvalue (optional)


=cut

sub label2{
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'label2'} = $value;
    }
    return $self->{'label2'};

}



=head2 color

 Title   : color
 Usage   : $obj->color($newval)
 Function: 
 Returns : value of color
 Args    : newvalue (optional)


=cut

sub color{
   my $self = shift;


   if( @_ ) {
      my $value = shift;
      $self->{'color'} = $value;
    }
    return $self->{'color'};

}

=head2 shape

 Title   : shape
 Usage   : $obj->shape($newval)
 Function: 
 Returns : value of shape
 Args    : newvalue (optional)


=cut

sub shape{
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'shape'} = $value;
    }
    return $self->{'shape'};

}



=head2 stretch

 Title   : stretch
 Usage   : $obj->stretch($newval)
 Function: gets/sets a boolean for whether we should stretch the data over the
 range (i.e. from min to max rather than absolute numbers).
 Returns : value of _stretch
 Args    : newvalue (optional)


=cut

sub stretch{
   my ($self,$value) = @_;
   if( defined $value ) {
      $self->{'_stretch'} = $value;
    }
    return $self->{'_stretch'};
}


=head2 scale_to_fit

 Title   : scale_to_fit
 Usage   : $obj->scale_to_fit($newval)
 Function: gets/sets the number that the BinValues are to be scaled against -
 i.e. the greatest BinValue->value will be scaled to this number, and the rest
 scaled in proportion.
 Returns : scale_to_fit value
 Args    : newvalue (optional)


=cut

sub scale_to_fit{
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'_scale_to_fit'} = $value;
    }
    return $self->{'_scale_to_fit'};

}


=head2 _store_biggest

 Title   : _store_biggest
 Usage   : $self->_store_biggest($newval)
 Function: internal method for storing the largest BinValue->value in this set.
 Returns : biggest value seen so far
 Args    : value 


=cut

sub _store_biggest {
    my ($self,$val) = @_;
    
    if (!defined $self->{'_biggest_value'} ||
	$val > $self->{'_biggest_value'}){
	$self->{'_biggest_value'}=$val;
    }

    return $self->{'_biggest_value'};
}



=head2 _store_smallest

 Title   : _store_smallest
 Usage   : $self->_store_smallest($newval)
 Function: internal method for storing the smallest BinValue->value in this set.
 Returns : smallest value seen so far
 Args    : value 

=cut

sub _store_smallest {
    my ($self,$val) = @_;
  
    if (!defined($self->{'_smallest_value'})){
	$self->{'_smallest_value'}=$val;
    }

    if (!defined($self->{'_smallest_value'}) ||
	$val < $self->{'_smallest_value'}){
	$self->{'_smallest_value'}=$val;
    }
    return $self->{'_smallest_value'};
}



1;
