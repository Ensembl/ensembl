# EnsEMBL RawContig object
#
# Copyright EMBL-EBI 2001
#
# 
# Date : 22.11.2001
#


=head1 NAME

Bio::EnsEMBL::RawContig
  Contig object which represents part of an EMBL Clone.Mostly for database usage

=head1 SYNOPSIS

=head1 CONTACT

  Arne Stabenau: stabenau@ebi.ac.uk
  Ewan Birney  : birney@ebi.ac.uk

=head1 APPENDIX

=cut



package Bio::EnsEMBL::RawContig;

use vars qw( @ISA );
use strict;


use Bio::EnsEMBL::DBSQL::DBAdaptor;

@ISA = qw( Bio::Root::RootI );


sub new {
  my ( $class, @args ) = @_;

  my $self = {};
  bless $self, $class;
  


}


1;
