
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

# Object preamble - inherits from Bio::Root::RootI

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::AssemblyMapper;


@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);


sub new {
  my($class, $dbadaptor) = @_;

  my $self = $class->SUPER::new($dbadaptor);

  $self->{'_type_cache'} = {};

  return $self;
}


=head2 fetch_by_type

  Arg  1      char $type
  Function    Fetches a Bio::EnsEMBL::AssemblyMapper object
              from the adaptor for a particular assembly
              (golden path) type (e.g. NCBI_xx). The result
              is cached for a particular assembly type.
  Returntype  Bio::EnsEMBL::AssemblyMapper
  Exceptions  none
  Caller      Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor

=cut

sub fetch_by_type{
   my ($self,$type) = @_;

   if( !defined $self->{'_type_cache'}->{$type} ) {
       $self->{'_type_cache'}->{$type} = Bio::EnsEMBL::AssemblyMapper->new($self,$type);
   }

   return $self->{'_type_cache'}->{$type};
}


=head2 register_region

  Arg  1      Bio::EnsEMBL::AssemblyMapper $assmapper
	      A valid AssemblyMapper object
  Arg  2      char $type
	      golden path type (e.g. NCBI_xx)
  Arg  3      char $chr_name
	      chromosome name (e.g. '2', 'X')
  Arg  4      int $start
	      chromosome start coordinate
  Arg  5      int $end
	      chromosome end coordinate
  Function    Declares a chromosomal region to the AssemblyMapper.
              This extracts the relevant data from the assembly
              table and stores it in a Bio::EnsEMBL::Mapper.
              It therefore must be called before any mapping is
              attempted on that region. Otherwise only gaps will
              be returned!
  Returntype  none
  Exceptions  none
  Caller      Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor

=cut

sub register_region{
   my ($self, $assmapper, $type, $chr_name, $start, $end) = @_;

   $self->throw("$assmapper is not a Bio::EnsEMBL::AssemblyMapper")
      unless $assmapper->isa("Bio::EnsEMBL::AssemblyMapper");


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
   #
   # SCP: Those not clauses are a bit confusing, so I've
   # ditched them :-)  "NOT a > b" equiv to "a <= b"

   # No. SCP - not right! The NOT has a bracket'd AND and so
   # one can't just expand it like this EB. Replaced it and
   # if *anyone* wants to change this, talk to me or James G first!


   my $select = qq{
      select
         ass.contig_start,
         ass.contig_end,
         ass.contig_id,
         ass.contig_ori,
         chr.name,
         ass.chr_start,
         ass.chr_end
      from
         assembly ass,
         chromosome chr
      where
         chr.name = '$chr_name' and
         ass.chromosome_id = chr.chromosome_id and
         NOT (
         $start >= ass.chr_end  and
         $end <= ass.chr_start)  and
         ass.type = '$type'
   };

   #print STDERR "Going to select $select\n";
   my $sth = $self->prepare($select);
   
   $sth->execute();
   
   while( my $arrayref = $sth->fetchrow_arrayref ) {
     my ($contig_start,$contig_end,$contig_id,$contig_ori,$chr_name,$chr_start,$chr_end) = @$arrayref;
     if( $assmapper->_have_registered_contig($contig_id) == 0 ) {
       $assmapper->_register_contig($contig_id);
       $assmapper->_mapper->add_map_coordinates(
						$contig_id, $contig_start, $contig_end, $contig_ori,
						$chr_name, $chr_start, $chr_end
					       );
     }
   }
   
}


=head2 register_region

  Arg  1      Bio::EnsEMBL::AssemblyMapper $assmapper
	      A valid AssemblyMapper object
  Arg  2      char $type
	      golden path type (e.g. NCBI_xx)
  Arg  3      int $contig_id
	      RawContig internal ID
  Arg  4      int $left
	      5 prime (chromosomal) extension
  Arg [5]     int $right
	      optional 3 prime extension
	      (same as 5 prime if not defined)
  Function    Declares a chromosomal region to the AssemblyMapper.
              This extracts the relevant data from the assembly
              table and stores it in a Bio::EnsEMBL::Mapper.
              It therefore must be called before any mapping is
              attempted on that region. Otherwise
              only gaps will
              be returned!
  Returntype  none
  Exceptions  none
  Caller      Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor

=cut

sub register_region_around_contig {
   my ($self, $assmapper, $type, $contig_id, $left, $right) = @_;

   $right = $left unless defined $right;

   my $sth = $self->prepare(qq{
      select
         chr_start,
         chr_end
      from
         assembly
      where
	 contig_id = '$contig_id' and
	 type = '$type'
   });

   $sth->execute();

   my @ctg_list;

   while (my $row = $sth->fetchrow_arrayref) {
       push @ctg_list, $row;
   }

   if (@ctg_list == 0) {
     $self->warn("Not found contig $contig_id");
     return;
   }

   if (@ctg_list > 1) {
     $self->warn("Contig $contig_id is ambiguous in assembly type $type");
     return;
   }

   my $start = $ctg_list[0]->[0] - $left;
   my $end   = $ctg_list[0]->[1] + $right;

   $sth = $self->prepare(qq{
      select
         ass.contig_start,
         ass.contig_end,
         ass.contig_id,
         ass.contig_ori,
         chr.name,
         ass.chr_start,
         ass.chr_end
      from
         assembly ass,
         chromosome chr
      where
         ass.chromosome_id = chr.chromosome_id and
         ass.chr_start <= $end and
         ass.chr_end   >= $start and
         ass.type = '$type'
   });

   $sth->execute();

   while( my $arrayref = $sth->fetchrow_arrayref ) {
       my ($contig_start,$contig_end,$contig_id,$contig_ori,$chr_name,$chr_start,$chr_end) = @$arrayref;
       if( $assmapper->_have_registered_contig($contig_id) == 0 ) {
	   $assmapper->_register_contig($contig_id);
	   $assmapper->_mapper->add_map_coordinates(
               $contig_id, $contig_start, $contig_end, $contig_ori,
	       $chr_name, $chr_start, $chr_end
           );
       }
   }

}

1;
