
#
# BioPerl module for Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor - DESCRIPTION of Object

=head1 SYNOPSIS

    # $dba is an Bio::EnsEMBL::DBSQL::DBAdaptor
    my $ass_adptr = $dba->get_AssemblyMapperAdaptor;

=head1 DESCRIPTION

Adaptor for handling Assembly mappers.  This is a
I<Singleton> class.  ie: There is only one per
database (C<DBAdaptor>).

=head1 AUTHOR - Ewan Birney

Email birney@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::RootI

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::AssemblyMapper;


@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);
# new() can be inherited from Bio::Root::RootI


sub new {
  my($class, $dbadaptor) = @_;

  my $self = $class->SUPER::new($dbadaptor);

  $self->{'_type_cache'} = {};

  return $self;
}
  

=head2 fetch_by_type

    my $ass_mapr = $ass_adptr->fetch_by_type('UCSC');

Fetches a C<Bio::EnsEMBL::AssemblyMapper> object
from the adaptor for a particular assembly
(golden path) type (eg: B<USCS> or B<SANGER>).

The answer is cached for a particular assembly type.

=cut

sub fetch_by_type{
   my ($self,$type) = @_;

   if( !defined $self->{'_type_cache'}->{$type} ) {
       $self->{'_type_cache'}->{$type} = Bio::EnsEMBL::AssemblyMapper->new($self,$type);
   }
   
   return $self->{'_type_cache'}->{$type};
}



=head2 register_region

    $ass_adptr->register_region($ass_mapr, 'UCSC', '20', 1_000_000, 2_000_000);



=cut

sub register_region{
   my ($self,$assmapper,$type,$chr_name,$start,$end) = @_;

   # This is to reassure people that I've got this piece of
   # SQL correct ;) EB
   #
   # JGRG: The NOT clauses are slightly counterintuitive,
   # but fetch every contig that overlaps start and end.
   # We search for every contig that doesn't begin after
   # our end (doesn't overlap) and doesn't end before
   # our start (doesn't overlap), and therefore get every
   # contig that overlaps.  This takes care of the condition
   # where our start and end lie within a contig.



   my $sth = $self->prepare("select ass.contig_start,ass.contig_end,ass.contig_id,ass.contig_ori,chr.name,ass.chr_start,ass.chr_end from assembly ass,chromosome chr where chr.name = '$chr_name' AND ass.chromosome_id = chr.chromosome_id and NOT (ass.chr_start > $end) and NOT (ass.chr_end < $start) and ass.type = '$type'");

   $sth->execute();

   while( my $arrayref = $sth->fetchrow_arrayref ) {
       my ($contig_start,$contig_end,$contig_id,$contig_ori,$chr_name,$chr_start,$chr_end) = @$arrayref;
       if( $assmapper->_have_loaded_contig($contig_id) == 0 ) {
	   $assmapper->_register_contig($contig_id);
	   $assmapper->_mapper->add_map_coordinates($contig_start,$contig_end,$contig_id,$chr_start,$chr_end,$chr_name,$contig_ori);
       }
   }

}




1;

