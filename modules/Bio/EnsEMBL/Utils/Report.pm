=head1 NAME - Bio::EnsEMBL::Utils::Report.pm

=head1 DESCRIPTION

=head1 SYNOPSIS

=head1 CONTACT

eae@sanger.ac.uk

=cut

# Let the code begin ...

package Bio::EnsEMBL::Utils::Report;

use vars qw(@ISA);
use strict;
use diagnostics;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::GeneCluster;
use Bio::EnsEMBL::Utils::TranscriptCluster;
use Bio::EnsEMBL::Utils::GeneComparison;
use Bio::EnsEMBL::Root;

@ISA = qw(Bio::EnsEMBL::Root);

####################################################################################

=head2 new()

=cut

sub new {
  
  my ($class) = @_;
  
  if (ref($class)){
    $class = ref($class);
  }
  my $self = {};
  bless($self,$class);
  return $self;
}

####################################################################################

=head2 header()

get/set method for a string header

=cut

sub header{
  my ($self,$header) = @_;
  if ($header){
    $self->{'_header'} = $header;
  }
  return $self->{'_header'};
}


####################################################################################

=head2 histogram()

get/set method for a hash, aimed to be used to create histograms

=cut

sub histogram{
  my ($self,$histogram) = @_;
  if ($histogram){
    $self->{'_histogram'} = $histogram;
  }
  return $self->{'_histogram'} = $histogram;
}

####################################################################################

=head2 list()

get/set method for an array of arrays, aimed to be used to print out lists

=cut

sub list{
  my ($self,$list) = @_;

  # $list is a reference to an array, whose elementes are arrays
  if ($list){
    $self->{'_list'} = $list;
  }
  return $self->{'_list'};
}

####################################################################################

=head2 to_plainTXT()

sends the info in the containers {'_list'}, {'_histogram'}, ... to plain text output
to file passed as argument, or as default, to SDTOUT

=cut

sub to_plainTXT(){
  my ($self,$filename) = @_;
  
  my $title     =    $self->title;
  my %histogram = %{ $self->histogram };
  my @list      = @{ $self->list      };

  # here we should be able to choose where to send the output
  # define OUTPUT handle
  #if ($filename){
  #  open OUT,">$filename";
  #}
  #else{
  #  OUT = STDOUT;
  #}
  print STDOUT $title."\n";
  foreach my $sublist ( @list ){
    foreach my $item ( @$sublist ){
      print STDOUT $item."\t";
    }
    print STDOUT "\n";
  }

  foreach my $key ( keys(%histogram) ){
    print STDOUT $key." ---> ".$histogram{$key}."\n";
  }
}

####################################################################################

=head2 histogram_to_R()

prints out the {'_histogram'} in format to be read by R and make a histogram

=cut

sub histogram_to_R(){
  ($self, $filename) = @_;
  %histogram = $self->histogram;
  
  print STDERR "printing histogram to file $filename...\n";

  # open filehandle
  my $opening = open OUT, ">$filename";;
  unless ($opening){
    $self->throw("Can't open file $filename");
  }

  # here we process the contents of %histogram and put it into $filename
  foreach my $key ( keys(%histogram) ){
    
    # ...
    # put here the right format
    print OUT $key."\t".$histogram{$key}."\n";

  }

  # close filehandle
  my $closing = close OUT;
  unless ($closing){
    $self->throw("Can't close file $filename");
  }
}

####################################################################################

=head2 to_HTML

sends the info in the containers {'_list'}, {'_histogram'}, ... to plain html format output

=cut

sub to_HTML{
  my ($self,$filename) = @_;
}

1;
