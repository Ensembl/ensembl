
#
# BioPerl module for BinValue
#
# Cared for by Emmanuel Mongin <mongin@ebi.ac.uk>
#
# Copyright Emmanuel Mongin
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

BinValue

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

This object deals with the raw data to built the density plots

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DensityPlot::BinValue;
use vars qw($AUTOLOAD @ISA);
use strict;

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
    return $self;
}

=head2 chromosomestart

 Title   : ChromosomeStart
 Usage   : $obj->ChromosomeStart($newval)
 Function: 
 Returns : value of ChromosomeStart
 Args    : newvalue (optional)


=cut

sub chromosomestart{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'chromosomestart'} = $value;
    }
    return $obj->{'chromosomestart'};

}

=head2 chromosomeend

 Title   : chromosomesnd
 Usage   : $obj->chromosomeend($newval)
 Function: 
 Returns : value of chromosomeend
 Args    : newvalue (optional)


=cut

sub chromosomeend{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'chromosomeend'} = $value;
    }
    return $obj->{'chromosomeend'};

}


=head2 value

 Title   : value
 Usage   : $obj->value($newval)
 Function: 
 Returns : value of value
 Args    : newvalue (optional)


=cut

sub value{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'value'} = $value;
    }
    return $obj->{'value'};

}



=head2 scaledvalue

 Title   : scaledvalue
 Usage   : $obj->scaledvalue($newval)
 Function: 
 Returns : this object's scaled value
 Args    : newvalue (optional)


=cut

sub scaledvalue{
   my $obj = shift;
   if( @_ ) {
      my $scaledvalue = shift;
      $obj->{'scaledvalue'} = $scaledvalue;
    }
    return $obj->{'scaledvalue'};

}



=head2 url

 Title   : url
 Usage   : $obj->url($newval)
 Function: 
 Returns : this object's url
 Args    : newvalue (optional)


=cut

sub url{
   my $obj = shift;
   if( @_ ) {
      my $url = shift;
      $obj->{'url'} = $url;
    }
    return $obj->{'url'};

}

1;
