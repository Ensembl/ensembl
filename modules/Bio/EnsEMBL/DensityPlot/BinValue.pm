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

Bio::EnsEMBL::DensityPlot::BinValue

=head1 SYNOPSIS

=head1 DESCRIPTION

This object deals with the raw data to built the density plots

=head1 METHODS

=cut

package Bio::EnsEMBL::DensityPlot::BinValue;
use vars qw($AUTOLOAD @ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

@ISA = qw(Exporter);
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
