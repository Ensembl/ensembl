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

Bio::EnsEMBL::Utils::Eprof - Bespoke Ensembl profiler

=head1 SYNOPSIS

  use Bio::EnsEMBL::Utils::Eprof( 'eprof_start', 'eprof_end',
    'eprof_dump' );

  &eprof_start('function-a');
  # ... do something
  &eprof_end('function-a');

  &eprof_dump( \*STDERR );

  # there is an object based set for above as well, for running
  # multiple concurrent profilers

=head1 DESCRIPTION

This is an Ensembl profiler as we broke the Perl profilers.

=head1 METHODS

=cut

package Bio::EnsEMBL::Utils::Eprof;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception ('throw');
use Bio::EnsEMBL::Utils::EprofStack;

use base qw( Exporter );

our @EXPORT_OK =
  ( 'eprof_start', 'eprof_end', 'eprof_dump', 'eprof_reset' );

my $global;

sub new {
  my ($proto) = @_;

  my $class = ref($proto) || $proto;
  my $self = bless( { '_tags' => {} }, $class );

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

sub eprof_start {
  my ($tag) = @_;

  if ( !defined($global) ) {
    $global = Bio::EnsEMBL::Utils::Eprof->new();
  }

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

  if ( !defined($global) ) {
    $global = Bio::EnsEMBL::Utils::Eprof->new();
  }

  $global->end($tag);
}

sub eprof_dump {
  my ($fh) = @_;

  if ( !defined($global) ) { return }

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

sub eprof_reset { undef($global) }

=head2 dump

 Title   : dump
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub dump {
  my ( $self, $fh ) = @_;

  my @tags = sort {
    $self->_tags()->{$a}->total_time()
      <=> $self->_tags()->{$b}->total_time()
  } keys %{ $self->_tags() };

  foreach my $tag (@tags) {
    my $st = $self->_tags->{$tag};

    if ( $st->number() == 0 ) { next }

    my $STD = '---';

    if ( $st->number() > 1 ) {
      my $SS =
        $st->total_time_time() -
        $st->total_time()*$st->total_time()/$st->number();

      if ( $SS > 0 ) {
        $STD = sprintf( "%6f",
                        sqrt( $SS/$st->number()/( $st->number() - 1 ) )
        );
      }
    }

    print( $fh sprintf( "Eprof: %20s  %6f  %6f  %d  %s  [%6f,%6f]\n",
                        $st->tag(), $st->total_time(),
                        $st->total_time()/$st->number(), $st->number(),
                        $STD, $st->min_time(),
                        $st->max_time() ) );
  } ## end foreach my $tag (@tags)
} ## end sub dump

=head2 start

 Title   : start
 Usage   : $eprof->start('this_tag');
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub start {
  my ( $self, $tag ) = @_;

  if ( !defined($tag) ) {
    $self->throw("No tag, can't start.");
  }

  if ( !defined( $self->_tags()->{$tag} ) ) {
    $self->_tags()->{$tag} = Bio::EnsEMBL::Utils::EprofStack->new($tag);
  }

  $self->_tags()->{$tag}->push_stack();
}

=head2 end

 Title   : end
 Usage   : $eprof->end('this_tag');
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub end {
  my ( $self, $tag ) = @_;

  if ( !defined($tag) ) {
    $self->throw("No tag, can't end.");
  }

  if ( !defined( $self->_tags()->{$tag} ) ) {
    $self->throw(
                sprintf( "Ending with a nonexistant tag '%s'", $tag ) );
  }

  $self->_tags->{$tag}->pop_stack();
}

=head2 _tags

 Title   : _tags
 Usage   : $obj->_tags($newval)
 Function: 
 Returns : value of _tags
 Args    : newvalue (optional)


=cut

sub _tags {
  my ($obj) = @_;
  return $obj->{'_tags'};
}

1;

