# EnsEMBL Gene reading writing adaptor for mySQL
#
# Copyright EMBL-EBI 2001
#
# Author: Graham McVicker
# 
# Date : 22.07.2002
#

=head1 NAME

Bio::EnsEMBL::DBSQL::ProxyGeneAdaptor

=head1 SYNOPSIS

Designed as an abstraction over the database specific GeneAdaptors.  The proxy
gene adaptor normally behaves just as a normal core GeneAdaptor, however, for
certain requests it may decide to instead forward the request to another
database (such as the lite database if it is available).

=head1 CONTACT

  Arne Stabenau: stabenau@ebi.ac.uk
  Ewan Birney  : birney@ebi.ac.uk
  Graham McVicker : mcvicker@ebi.ac.uk

=head1 APPENDIX

=cut

use strict;

package Bio::EnsEMBL::DBSQL::ProxyGeneAdaptor;

use Bio::EnsEMBL::DBSQL::GeneAdaptorI;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;

use vars '@ISA';

@ISA = qw(Bio::EnsEMBL::DBSQL::GeneAdaptorI Bio::EnsEMBL::DBSQL::BaseAdaptor);

#implement the interface GeneAdaptorI
use implements qw(Bio::EnsEMBL::DBSQL::GeneAdaptorI);

sub new {
  my($class, $db, $core_adaptor) = @_;

  #call superclass constructor
  my $self = $class->SUPER::new($db);
  $self->{'_core_adaptor'} = $core_adaptor;
  return $self;
}

sub fetch_by_Slice {
  my ($self, @args) = @_;

  my $lite_db = $self->db()->lite_DBAdaptor();
 
  if(defined $lite_db) {
    #use the Lite database if it is available
    return $lite_db->get_GeneAdaptor()->fetch_by_Slice(@args);
  }

  #otherwise use the core database
  return $self->{'_core_adaptor'}->fetch_by_Slice(@args);
}


sub fetch_by_transcript_stable_id {
  my ($self, @args) = @_;

  my $lite_db = $self->db()->lite_DBAdaptor();
  
  if(defined $lite_db) {
    #use the Lite database if it is available
    return $lite_db->get_GeneAdaptor()->fetch_by_transcript_stable_id(@args);
  }
  
  #otherwise use the core database
  return $self->{'_core_adaptor'}->fetch_by_transcript_stable_id(@args);
}


sub list_geneIds {
   my ($self, @args) = @_;

   #use core db
   return $self->{'_core_adaptor'}->list_geneIds(@args);
}

sub list_stable_geneIds {
   my ($self, @args) = @_;

   #use core db
   return $self->{'_core_adaptor'}->list_stable_geneIds(@args);
}

sub fetch_by_dbID {
  my ( $self, @args) = @_;

  #use core db
  return $self->{'_core_adaptor'}->fetch_by_dbID(@args);
}

sub fetch_by_stable_id{
  my ($self, @args) = @_;
  print STDERR ( "ProxyGeneAdaptor is called.\n" );

  my $lite_db = $self->db()->lite_DBAdaptor();
  
  if(defined $lite_db) {
    #use the Lite database if it is available
    print STDERR "Lite Database used for fetch_by_stable_id.\n" ; 
    return $lite_db->get_GeneAdaptor()->fetch_by_stable_id(@args);
  }
  
  #otherwise use the core database
  return $self->{'_core_adaptor'}->fetch_by_stable_id(@args);
}

sub fetch_by_contig_list{
  my ($self, @args) = @_;

  #use core db
  return $self->{'_core_adaptor'}->fetch_by_contig_list(@args);
}


sub fetch_by_Transcript_id {
    my ($self, @args) = @_;

    #use core db
    return $self->{'_core_adaptor'}->fetch_by_Transcript_id(@args);
}

sub fetch_by_Peptide_id {
    my ($self, @args) = @_;
    
    #use core db
    return $self->{'_core_adaptor'}->fetch_by_Peptise_id(@args);
}

sub fetch_by_maximum_DBLink {
    my ($self, @args) = @_;

    #use core db
    return $self->{'_core_adaptor'}->fetch_by_maximum_DBLink(@args);
}


sub get_description {
    my ($self, @args) = @_;

    #use core db
    return $self->{'_core_adaptor'}->get_description(@args);
}


sub get_stable_entry_info {
    my ($self, @args) = @_;

    #use core db
    return $self->{'_core_adaptor'}->get_stable_entry_info(@args);
}


sub fetch_by_DBEntry {
    my ($self, @args) = @_;

    #use core db
    return $self->{'_core_adaptor'}->fetch_by_DBEntry(@args);
}


sub store {
    my ($self, @args) = @_;

    return $self->{'_core_adaptor'}->store(@args);
}


sub remove {
    my ($self, @args) = @_;

    return $self->{'_core_adaptor'}->remove(@args);
}


=head2 get_Interpro_by_geneid

 Title   : get_Interpro_by_geneid
 Usage   : @interproid = $geneAdaptor->get_Interpro_by_geneid($gene->id);
 Function: gets interpro accession numbers by geneid. A hack really -
           we should have a much more structured system than this
 Example :
 Returns : 
 Args    :


=cut

sub get_Interpro_by_geneid {
    my ($self, @args) = @_;

    return $self->{'_core_adaptor'}->get_Interpro_by_geneid(@args);
}


sub create_tables {
  my ($self, @args) = @_;

  return $self->{'_core_adaptor'}->create_tables(@args);
}


