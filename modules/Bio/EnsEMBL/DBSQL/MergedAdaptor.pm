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

Bio::EnsEMBL::DBSQL::MergedAdaptor

=head1 SYNOPSIS

  $merged_adaptor = new Bio::EnsEMBL::DBSQL::MergedAdaptor(
    -species => "human",
    -type    => "Population"
  );

=head1 DESCRIPTION

The MergedAdaptor object is merely a list of adaptors. AUTOLOAD is used
to call a subroutine on each adaptor and merge the results.

=head1 METHODS

=cut


package Bio::EnsEMBL::DBSQL::MergedAdaptor;


use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Utils::Exception qw(throw warning deprecate);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Registry;
my $reg = "Bio::EnsEMBL::Registry";


=head2 new

  Arg [SPECIES]: (optional) string 
                  species name to get adaptors for
  Arg [TYPE]   : (optional) string 
                  type to get adaptors for
  Arg [GROUPS] : (optional) ref to list
  Example      : $MergedAdaptor = new 
               : Bio::EnsEMBL::DBSQL::MergedAdaptor(-species=> 'human', -type =>'Population', -groups => ['Sanger','Ensembl']);

  Description: Creates a new MergedAdaptor
  Returntype : Bio::EnsEMBL::DBSQL::MergedAdaptor
  Exceptions : throws if species or type not specified
  Caller     : general
  Status     : At Risk
             : Under development

=cut

sub new {
  my ($class,@args) = @_;

  my $self ={};
  bless $self,$class;

  my ($species, $type, $groups) =
    rearrange([qw(SPECIES TYPE GROUPS)], @args);

  if(!defined($species)|| !defined($type)){
    die "Species and Type must be specified\n";
  }

  my @adaps;
  if (!defined ($groups)){
      #get all adaptors for that species and type
      @adaps = @{$reg->get_all_adaptors(-species => $species, -type => $type)};
  }
  else{
      #get only specified adaptors for the particular groups
      foreach my $group (@{$groups}){
	  push @adaps, $reg->get_adaptor($species,$group,$type);
      }
  }
 
  my @list =();
  push(@list,@adaps);
  $self->{'list'}= \@list;

  return $self;
}

=head2 add_list

  Example    : $MergedAdaptor->add_list(@adaptors);
  Description: adds a list of adaptors to the Merged adaptor list.
  Returntype : none
  Exceptions : none
  Status     : At Risk
             : Under development

=cut

sub add_list{
  my ($self, @arr) = @_;

  foreach my $adap (@arr){
    $self->add_adaptor($adap);
  }
}

=head2 add_adaptor

  Example    : $MergedAdaptor->add_adaptor(@adaptors);
  Description: adds an adaptor to the Merged adaptor list.
  Returntype : none
  Exceptions : none
  Status     : At Risk
             : Under development

=cut

sub add_adaptor{
  my ($self,$adaptor)=@_;

  if(!defined ($self->{'list'})){
    my @list =();
    push(@list,$adaptor);
    $self->{'list'}= \@list;
  }
  else{
    push(@{$self->{'list'}},$adaptor);
  }
}


sub printit{
  my ($self)=@_;

  foreach my $adaptor (@{$self->{'list'}}){
    print "printit $adaptor\t".$adaptor->db->group()."\n";
  }
}


use vars '$AUTOLOAD';

sub AUTOLOAD {
  my ($self,@args) = @_;
  my @array_return=();
  my $ref_return = undef;
  $AUTOLOAD =~ /^.*::(\w+)+$/ ;

  my $sub = $1;

  foreach my $adaptor (@{$self->{'list'}}) {
    my $ref;
    if($adaptor->can($sub)){
      $ref = $adaptor->$sub(@args);
      if( ref($ref) eq 'ARRAY' ) {
        push @array_return, @{$ref};
      } else {
        push @array_return, $ref;
      }
    }
    else{ # end of can
      warn("In Merged Adaptor $adaptor cannot call sub $sub");
    }
  }
  return \@array_return;
}

sub DESTROY{
}

1;
