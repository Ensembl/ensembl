#
# BioPerl module for Bio::EnsEMBL::Utils::Eprof
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright EBI and GRL
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Utils::Eprof - Bespoke Ensembl profiler

=head1 SYNOPSIS

    use Bio::EnsEMBL::Utils::Eprof('eprof_start','eprof_end','eprof_dump');

    &eprof_start('function-a');
    ... do something
    &eprof_end('function-a');


    &eprof_dump(\*STDERR);

    # there is an object based set for above as well, for running
    # multiple concurrent profilers


=head1 DESCRIPTION

This is an Ensembl profiler as we broke the Perl profilers.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Utils::Eprof;
use vars qw(@ISA @EXPORT_OK);
use strict;
use Exporter;
use Bio::EnsEMBL::Utils::EprofStack;

# Object preamble - inheriets from Bio::Root::Object

use Bio::EnsEMBL::Root;

@ISA = qw(Bio::EnsEMBL::Root Exporter);
@EXPORT_OK = qw(eprof_start eprof_end eprof_dump eprof_reset);

my $global;

sub new { 
    my ($class) = shift;
    my $self = {};
    $self->{'_tags'} = {};
    bless $self,$class;

    return $self;
}

=head2 eprof_start

 Title   : eprof_start
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub eprof_start{
   my ($tag) = @_;
   $global = Bio::EnsEMBL::Utils::Eprof->new() unless defined $global;
   $global->start($tag);
}

=head2 eprof_end

 Title   : eprof_end
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub eprof_end {
   my ($tag) = @_;
   $global = Bio::EnsEMBL::Utils::Eprof->new() unless defined $global;
   $global->end($tag);
}

sub eprof_dump {
    my $fh = shift;

    if( !defined $global ) {
	return;
    }

    $global->dump($fh);
}

=head2 eprof_reset

 Title   : eprof_reset
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub eprof_reset{
  undef($global);
}

=head2 dump

 Title   : dump
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub dump{
  my ($self,$fh) = @_;

  my @tags = sort {  $self->_tags->{$a}->total_time <=> $self->_tags->{$b}->total_time } keys %{$self->_tags};
   
  foreach my $tag ( @tags ) {
    my $st = $self->_tags->{$tag};
    next if $st->number == 0;
    my $STD = '---';
    if($st->number>1) {
      my $SS = $st->total_time_time - $st->total_time*$st->total_time/$st->number;
      $STD = sprintf "%6f", sqrt( $SS/$st->number/($st->number-1) ) if $SS>0;
    }
    print $fh sprintf("Eprof: %20s  %6f  %6f  %d  %s  [%6f,%6f]\n",$st->tag,$st->total_time,$st->total_time/$st->number,$st->number,$STD,$st->min_time,$st->max_time);
  }
}


=head2 start

 Title   : start
 Usage   : $eprof->start('this_tag');
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub start{
  my ($self,$tag) = @_;
  $self->throw("Must start on tag")                                 unless defined $tag;
  $self->_tags->{$tag} = Bio::EnsEMBL::Utils::EprofStack->new($tag) unless defined $self->_tags->{$tag};
  $self->_tags->{$tag}->push_stack();
}

=head2 end

 Title   : end
 Usage   : $eprof->end('this_tag');
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub end{
  my ($self,$tag) = @_;
  $self->throw("Must end on tag")               unless defined $tag;
  $self->throw("Ending with a nonexistant tag") unless defined $self->_tags->{$tag};
  $self->_tags->{$tag}->pop_stack();
}


=head2 _tags

 Title   : _tags
 Usage   : $obj->_tags($newval)
 Function: 
 Returns : value of _tags
 Args    : newvalue (optional)


=cut

sub _tags{
  my $obj = shift;
  return $obj->{'_tags'};
}

1;

