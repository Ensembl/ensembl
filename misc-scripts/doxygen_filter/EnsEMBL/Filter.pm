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


# =======================================================================
# Doxygen Pre-Processor for Perl
# Copyright (C) 2002  Bart Schuller
# Copyright (C) 2006  Phinex Informatik AG
# All Rights Reserved
# 
# Doxygen Filter is free software; you can redistribute it and/or modify
# it under the same terms as Perl itself.
# 
# Larry Wall's 'Artistic License' for perl can be found in
# http://www.perl.com/pub/a/language/misc/Artistic.html
# 
# =======================================================================
# 
# Author: Aeby Thomas, Phinex Informatik AG,
# 	  Based on DoxygenFilter from Bart Schuller
# E-Mail: tom.aeby@phinex.ch
# 
# Phinex Informatik AG
# Thomas Aeby
# Kirchweg 52
# 1735 Giffers
# 
# =======================================================================
# 
# @(#) $Id$
# 
# Revision History prior to import into EnsEMBL:
# 
# Revision 1.1  2011-10-18 13:46:34  kt7
# Transferring Doxygen filter into EnsEMBL CVS.
#
# Revision 1.1  2011-05-16 15:31:44  kt7
# Preservative backup of all code required to run Doxygen to
# generate API reference.
#
# Revision 1.2  2006/01/31 16:53:52  aeby
# added copyright info
#
#  
# =======================================================================

## @file
# implementation of Filter as derrived from Doxygen::Filter.


## @class
# Filter from non-C++ syntax API docs to Doxygen-compatible syntax.
# This class is meant to be used as a filter for the
# <a href="http://www.doxygen.org/">Doxygen</a> documentation tool.
package EnsEMBL::Filter;

use warnings;
use strict;

my $debug_home_location = '';
if($ENV{PERL_DOXYGEN_DEBUG_LOCATION}) {
  $debug_home_location = $ENV{PERL_DOXYGEN_DEBUG_LOCATION};
}

## @cmethod object new($outfh)
# create a filter object.
# @param outfh optional output filehandle; defaults to STDOUT
# @return filter object
sub new {
    my $class = shift;
    my $outfh = shift || \*STDOUT;
    
    my $self = bless {outfh => $outfh}, $class;
    
    if($debug_home_location) {
      open(my $debug, '>', $debug_home_location);
      $self->{debug_fh} = $debug;
    }
    
    return $self;
}

sub DESTROY {
  my ($self) = @_;
  close $self->{debug_fh} if $self->{debug_fh};
  return;
}

## @method virtual void filter($infh)=0
# do the filtering.
# @param infh input filehandle, normally STDIN
sub filter {
    die "subclass responsibility";
}

## @method protected string protection($sig)
# Return the protection of a method/function signature.
# @param sig the method signature
# @return Either "Public" or "Private".
sub protection {
    my($self, $sig) = @_;
    return $sig =~ /^(private|protected)/ ? "\u$1" : 'Public';
}

## @method void start($command)
# start a doc comment.
# Outputs the start of a javadoc comment.
# @param command the javadoc command
sub start {
    my $self = shift;
    my $command = shift;
    $self->print("/** $command\n");
    return $self;
}

## @method void end()
# end a doc comment.
# Outputs the end of a javadoc comment.
sub end {
    my ($self, $tag) = @_;
    if (defined $tag)
    {
	$self->print("$tag */\n");
    }
    else
    {
	$self->print("*/\n");
    }
    return $self;
}

## @method void push($section)
# Start a diversion to a section.
# @param section The name of the section to divert all output to.
# @see pop(), print(), flush()
sub push {
    my($self, $section) = @_;
    $self->{current_section} = $section;
    return $self;
}

## @method void pop()
# End a diversion to a section.
# @see push(), flush()
sub pop {
    my($self) = @_;
    delete $self->{current_section};
    return $self;
}

## @method void print(@args)
# print a string to the output handle.
# If a diversion to a specific section is in effect: saves the text under
# that section.
# @param args the strings to be printed
# @see push(), flush()
sub print {
    my $self = shift;
    return unless @_;
    if (my $section = $self->{current_section}) {
        CORE::push @{$self->{sections}{$section}}, @_;
    } else {
        my $outfh = $self->{outfh};
        print $outfh @_;
        my $debug_fh = $self->{debug_fh};
        print $debug_fh @_ if $debug_fh;
    }
    return $self;
}

sub echo {
    my $self = shift;
    return unless @_;
    print STDERR @_;
    return $self;
}

## @method void more(@args)
# process the follow-up lines after the initial apidoc line.
# @param args the lines to be processed
sub more {
    my $self = shift;
    $self->print(@_);
    return $self;
}

my @order = (
    'Public Class Methods',
    'Public Object Methods',
    'Public Functions',
    'Protected Class Methods',
    'Protected Object Methods',
    'Protected Functions',
    'Private Class Methods',
    'Private Object Methods',
    'Private Functions',
    );
## @method void flush()
# Flush the saved sections. Should be called at the end of a class.
# @see push(), print()
sub flush {
    my $self = shift;
    my $sections = $self->{sections};
    foreach (@order) {
        next unless $sections->{$_};
        $self->start("\@name $_\n")->end;
        $self->start("\@{")->end;
	$self->print("\n");
        $self->print(@{$sections->{$_}});
	$self->print("\n");
        $self->start("\@}")->end;
	$self->print("\n");
    }
    delete $self->{sections};
    return $self;
}

1;

__END__

=head1 NAME

EnsEMBL::Filter - use DoxyGen with Perl and other languages.

=head1 DESCRIPTION

Filter from non-C++ syntax API docs to Doxygen-compatible syntax.
This class is meant to be used as a filter for the
Doxygen (http://www.doxygen.org/) documentation tool

