
#
# Ensembl module for Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor
#
# Written by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright GRL/EBI
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

=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::EnsEMBL::Root

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::AssemblyMapper;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

=head2 new

  Arg [1]    : Bio::EnsEMBL::DBAdaptor $dbadaptor the adaptor for the database
               this assembly mapper is using.
  Example    : my $asma = new Bio::EnsEMBL::AssemblyMapperAdaptor($dbadaptor);
  Description: Creates a new AssemblyMapperAdaptor object
  Returntype : Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBSQL::DBAdaptor

=cut

sub new {
  my($class, $dbadaptor) = @_;

  my $self = $class->SUPER::new($dbadaptor);

  $self->{'_type_cache'} = {};

  return $self;
}



=head2 fetch_by_type

  Arg [1]    : char $type 
  Description: Fetches a Bio::EnsEMBL::AssemblyMapper object from the adaptor 
               for a particular assembly (golden path) type (e.g. NCBI_xx). 
               The result is cached for a particular assembly type. 
  Returntype : Bio::EnsEMBL::AssemblyMapper 
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor 

=cut

sub fetch_by_type{
  my ($self,$type) = @_;
  
  if( !defined $self->{'_type_cache'}->{$type} ) {
    $self->{'_type_cache'}->{$type} = 
      Bio::EnsEMBL::AssemblyMapper->new($self,$type);
  }
  
  return $self->{'_type_cache'}->{$type};
}


=head2 register_region

  Arg [1]    : Bio::EnsEMBL::AssemblyMapper $assmapper
	       A valid AssemblyMapper object
  Arg [2]    : char $type
	       golden path type (e.g. NCBI_xx)
  Arg [3]    : char $chr_name
	       chromosome name (e.g. '2', 'X')
  Arg [4]    : int $start
	       chromosome start coordinate
  Arg [5]    : int $end
	       chromosome end coordinate
  Description: Declares a chromosomal region to the AssemblyMapper.
               This extracts the relevant data from the assembly
               table and stores it in a Bio::EnsEMBL::Mapper.
               It therefore must be called before any mapping is
               attempted on that region. Otherwise only gaps will
               be returned!
  Returntype : none
  Exceptions : thrown if $assmapper arg is not a Bio::EnsEMBL::AssemblyMapper
  Caller     : Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor

=cut

sub register_region{
  my ($self, $assmapper, $type, $chr_name, $start, $end) = @_;

  $self->throw("$assmapper is not a Bio::EnsEMBL::AssemblyMapper")
    unless $assmapper->isa("Bio::EnsEMBL::AssemblyMapper");

  my $chr = $self->db->get_ChromosomeAdaptor()->fetch_by_chr_name( $chr_name );
  my $chr_id = $chr->dbID();
  my $max_assembly_contig = $self->db()->get_MetaContainer->get_max_assembly_contig();

  if (! defined $max_assembly_contig) {
      $max_assembly_contig = 1_000_000_000;
  }

  my $select = qq{
      select
         ass.contig_start,
         ass.contig_end,
         ass.contig_id,
         ass.contig_ori,
         ass.chr_start,
         ass.chr_end
      from
         assembly ass
      where
         ass.chromosome_id = $chr_id and
         ? <= ass.chr_end  and
         ? >= ass.chr_start  and
	 ? <= ass.chr_start and
         ass.type = '$type'
   };

  my $sth = $self->prepare($select);
   
  $sth->execute( $start, $end, $start-$max_assembly_contig );
   
  while( my $arrayref = $sth->fetchrow_arrayref ) {
    my ($contig_start, $contig_end, $contig_id, $contig_ori,
	$chr_start, $chr_end) = @$arrayref;
    if( $assmapper->_have_registered_contig($contig_id) == 0 ) {
      $assmapper->_register_contig($contig_id);
      $assmapper->_mapper->add_map_coordinates(
					       $contig_id, $contig_start, 
					       $contig_end, $contig_ori,
					       $chr_name, $chr_start, $chr_end
					      );
    }
  } 
}


=head2 register_contig

  Arg  [1]   : Bio::EnsEMBL::AssemblyMapper $assmapper
	       A valid AssemblyMapper object
  Arg  [2]   : char $type
	       golden path type (e.g. NCBI_xx)
  Arg  [3]   : int $contig_id
	       RawContig internal ID

  Description: Declares a chromosomal region to the AssemblyMapper.
               This extracts the relevant data from the assembly
               table and stores it in a Bio::EnsEMBL::Mapper.
               It therefore must be called before any mapping is
               attempted on that region. Otherwise only gaps will be returned!
  Returntype : 1 if the contig is present in assembly, 
               0 if the contig is absent
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor

=cut

sub register_contig {
   my ($self, $assmapper, $type, $contig_id ) = @_;

   my $sth = $self->prepare(qq{
      select
	c.name,
	a.chr_start,
	a.chr_end
      from
	assembly a,
	chromosome c
      where
	 contig_id = '$contig_id' and
	 type = '$type' and
	 c.chromosome_id = a.chromosome_id
   });

   $sth->execute();

   my @ctg_list;

   while (my $row = $sth->fetchrow_arrayref) {
       push @ctg_list, $row;
   }

   if (@ctg_list == 0) {
     #$self->warn("Not found contig $contig_id");
     return ();
   }

   if (@ctg_list > 1) {
     $self->warn("Contig $contig_id is ambiguous in assembly type $type");
     return ();
   }

   my $chr_name = $ctg_list[0]->[0];
   my $chr_start = $ctg_list[0]->[1];
   my $chr_end = $ctg_list[0]->[2];

   return ($chr_name, $chr_start, $chr_end);
}




1;
