
#
# BioPerl module for Bio::EnsEMBL::PerlDB::Clone
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::PerlDB::Clone - Pure Perl implementation of Clone

=head1 SYNOPSIS

    $clone = Bio::EnsEMBL::PerlDB::Clone->new();
 
    $clone->add_Contig($contig);
    

=head1 DESCRIPTION

Describe the object here

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::PerlDB::Clone;
use vars qw($AUTOLOAD @ISA);
use strict;
use Bio::EnsEMBL::DB::CloneI;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::RootI;


@ISA = qw(Bio::Root::RootI Bio::EnsEMBL::DB::CloneI);

sub new {
  my($class,@args) = @_;

  my $self = bless {
      _contig_hash     => {},
      _id              => undef,
      _embl_id         => undef,
      _embl_version    => undef,
      _htg_phase       => 1,
      _created         => time,
      _modified        => time,
      _seq_date        => 0,
      _contig_overlaps => [],
  }, $class;

 return $self; # success - we hope!
}

=head2 get_Contig

 Title   : get_Contig
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Contig{
   my ($self,$id) = @_;

   return $self->{'_contig_hash'}->{$id};
}


=head2 get_all_Contigs

 Title   : get_all_Contigs
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_Contigs{
   my ($self) = @_;

   return values %{$self->{'_contig_hash'}};
}

=head2 add_Contig

 Title   : add_Contig
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_Contig{
   my ($self,$contig) = @_;


   if( ! $contig->isa('Bio::EnsEMBL::DB::ContigI') ) {
       $self->warn("$contig is not a contigI object...");
   }

   $self->{'_contig_hash'}->{$contig->id()} = $contig;
}


=head2 get_all_Genes

 Title   : get_all_Genes
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_Genes{
   my ($self) = @_;
   my %h;

   # read into a hash to make unique
   foreach my $contig ( $self->get_all_Contigs ) {
       foreach my $gene ( $contig->get_all_Genes ) {
	   $h{$gene->id()} = $gene;
       }

   }

   return values %h;

}

=head2 id

 Title   : id
 Usage   : $obj->id($newval)
 Function: 
 Example : 
 Returns : value of id
 Args    : newvalue (optional)


=cut

sub id{
   my ($obj,$value) = @_;
   if( defined $value) {
       $obj->{'_id'} = $value;
   }
   return $obj->{'_id'};
   
}


=head2 embl_id

 Title   : embl_id
 Usage   : $obj->embl_id($newval)
 Function: 
 Returns : value of embl_id
 Args    : newvalue (optional)


=cut

sub embl_id{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_embl_id'} = $value;
    }
    return $obj->{'_embl_id'};

}

=head2 embl_version

 Title   : embl_version
 Usage   : $obj->embl_version($newval)
 Function: 
 Returns : value of embl_version
 Args    : newvalue (optional)


=cut

sub embl_version{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_embl_version'} = $value;
    }
    return $obj->{'_embl_version'};

}

=head2 sv

 Title   : sv
 Function: returns the version number (not the acc.version, just verision).
 Example :
 Returns : 
 Args    :


=cut

sub sv {
   my ($self,@args) = @_;
   return $self->version(@args);
}

=head2 seq_date

 Title   : seq_date
 Usage   : $obj->seq_date($newval)
 Function: 
 Returns : value of seq_date
 Args    : newvalue (optional)


=cut

sub seq_date{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_seq_date'} = $value;
    }
    return $obj->{'_seq_date'};

}

=head2 version

 Title   : version
 Usage   : $obj->version($newval)
 Function: 
 Returns : value of version
 Args    : newvalue (optional)


=cut

sub version{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_version'} = $value;
    }
    return $obj->{'_version'};

}
=head2 htg_phase

 Title   : htg_phase
 Usage   : $obj->htg_phase($newval)
 Function: 
 Returns : value of htg_phase
 Args    : newvalue (optional)


=cut

sub htg_phase{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_htg_phase'} = $value;
    }
    return $obj->{'_htg_phase'};

}

=head2 modified

 Title   : modified
 Usage   : $clone->modified()
 Function: Gives the unix time value of the modified 
           datetime field, which indicates
           the last time this clone was modified in ensembl
 Example : $clone->modified()
 Returns : unix time
 Args    : none


=cut

sub modified{
    my ($self,$value) = @_;
    if( defined($value) && $value ne '' ) {
      $self->{'_modified'} = $value;
    }
    return $self->{'_modified'};
}

=head2 created

 Title   : created
 Usage   : $clone->created()
 Function: Gives the unix time value of the created 
           datetime field, which indicates
           the first time this clone was put in ensembl
 Example : $clone->created()
 Returns : unix time
 Args    : none


=cut
sub created{
    my ($self,$value) = @_;
    if( defined($value) && $value ne '' ) {
      $self->{'_created'} = $value;
    }
    return $self->{'_created'};
}

=head2 get_all_ContigOverlaps 

 Title   : get_all_ContigOverlaps
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_ContigOverlaps {
    my ($self) = @_;
    
    return @{$self->{'_contig_overlaps'}};
}

1;









