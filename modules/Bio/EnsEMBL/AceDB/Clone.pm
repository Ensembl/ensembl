
#
# BioPerl module for DB::Clone
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::AceDB::Clone - Object representing one clone

=head1 SYNOPSIS

    # $db is Bio::EnsEMBL::AceDB::Obj 

    $clone = $db->get_Clone();

    @contig = $clone->get_Contigs();

    @genes  = $clone->get_all_Genes();

=head1 DESCRIPTION

Represents information on one Clone

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::AceDB::Clone;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;


@ISA = qw(Bio::Root::Object Bio::EnsEMBL::DB::CloneI );
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize(@args);

  # set stuff in self from @args
  my ($dbobj,$id) = $self->_rearrange([qw(DBOBJ
					  ID
					  )],@args);

  $id || $self->throw("Cannot make contig db object without id");
  $dbobj || $self->throw("Cannot make contig db object without db object");
  $dbobj->isa('Bio::EnsEMBL::AceDB::Obj') || 
      $self->throw("Cannot make contig db object with a $dbobj object");

  $self->id($id);
  $self->_dbobj($dbobj);

  return $make; # success - we hope!
}


=head2 seq

 Title   : seq
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=cut

sub seq {
   my ($self) = @_;

   my ($contig) = $self->get_Contig($self->id());
   return $contig->seq();
}


=head2 created

 Title   : created
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=cut

sub created {
    my ($self) = @_;
    my ($contig) = $self->get_Contig($self->id()); 
    if (my $date = $contig->ace_seq->at('Properties.Status.Finished[1]')) {
        return $date;
    }
    return 0;   
}


=head2 embl_version

 Title   : embl_version
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=cut

sub embl_version {
   my ($self) = @_;
   my ($contig) = $self->get_Contig($self->id());
   if (my $version = $contig->ace_seq->at('DB_info.Sequence_version[1]')) {
        return $version->name;
   }
   return;
}


=head2 embl_id

 Title   : embl_id
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=cut

sub embl_id {
   my ($self) = @_;

   my ($contig) = $self->get_Contig($self->id());
   if (my $database = $contig->ace_seq->at('DB_info.Database[1]')) {

        if ($database eq "EMBL") {
            if (my $embl_id = $contig->ace_seq->at('DB_info.Database[2]')) {
          
                return $embl_id->name;
            }
        }       
   }

   return $self->id;

}


=head2 htg_phase

 Title   : htg_phase
 Usage   : $obj->id($newval)
 Function: 
 Example : 
 Returns : value of id
 Args    : newvalue (optional)

=cut

sub htg_phase {
   my ($obj) = @_;
    return 3;
}


=head2 id

 Title   : id
 Usage   : $obj->id($newval)
 Function: 
 Example : 
 Returns : value of id
 Args    : newvalue (optional)

=cut

sub id {
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'_clone_id'} = $value;
    }
    return $obj->{'_clone_id'};

}


=head2 modified

 Title   : modified
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=cut

sub modified {
    my ($self) = @_;
    my ($contig) = $self->get_Contig($self->id()); 
    if (my $date = $contig->ace_seq->at('Properties.Status.Finished[1]')) {
        return $date->name;
    }
    return 0;   
}


=head2 seq_date

 Title   : seq_date
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=cut

sub seq_date {
    my ($self) = @_;
    my ($contig) = $self->get_Contig($self->id()); 
    return $contig->seq_date();   
}


=head2 sv

 Title   : sv
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=cut

sub sv {
    my ($self) = @_;
    return 1;  
}


=head2 version

 Title   : version
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=cut

sub version {
   my ($self) = @_;
 
   my ($contig) = $self->get_Contig($self->id()); 
   if (my $version = $contig->ace_seq->at('DB_info.Sequence_version[1]')) {   
        return $version->name;
   }
   # If the version isn't defined just return 1.     
   return 1;
}


=head2 get_all_Contigs

 Title   : get_Contigs
 Usage   : foreach $contig ( $clone->get_Contigs ) 
 Function:
 Example :
 Returns : 
 Args    :

=cut

sub get_all_Contigs {
   my ($self) = @_;
   my $contig = new Bio::EnsEMBL::AceDB::Contig ( -dbobj => $self->_dbobj,
					   '-id' => $self->id(), 
                                           '-clone' => $self);
                                           
   my @contigs;
   push(@contigs, $contig);
   return @contigs;   
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
    
    my @overlaps;

    foreach my $contig ($self->get_all_Contigs) {
	if (defined($contig->get_left_overlap)) {
	    
	    my $overlap    = $contig->get_left_overlap;
	    my $type;
	    
	    if ($overlap->sister_polarity == 1) {
		$type = 'left2right';
	    } 
            elsif ($overlap->sister_polarity == -1) {
		$type = 'left2left';
	    } 
            else {
		$self->throw("Invalid value [" .$overlap->sister_polarity . "] for polarity");
	    }
	    
	    my $tmpoverlap = new Bio::EnsEMBL::ContigOverlap(-contiga => $contig,
							     -contigb => $overlap->sister,
							     -positiona => $overlap->self_position,
							     -positionb => $overlap->sister_position,
							     -source    => $overlap->source,
							     -distance  => $overlap->distance,
							     -overlap_type => $type);
	    
	    push(@overlaps,$tmpoverlap);
	}

	if (defined($contig->get_right_overlap)) {
	    
	    my $overlap    = $contig->get_right_overlap;
	    my $type;
	    
	    if ($overlap->sister_polarity == 1) {
		$type = 'right2left';
	    } 
            elsif ($overlap->sister_polarity == -1) {
		$type = 'right2right';
	    } 
            else {
		$self->throw("Invalid value [" .$overlap->sister_polarity . "] for polarity");
	    }
	    
	    my $tmpoverlap = new Bio::EnsEMBL::ContigOverlap(-contiga => $contig,
							     -contigb => $overlap->sister,
							     -positiona => $overlap->self_position,
							     -positionb => $overlap->sister_position,
							     -source    => $overlap->source,
							     -distance  => $overlap->distance,
							     -overlap_type => $type);
	    
	    push(@overlaps,$tmpoverlap);
	    
	}
    }

    return (@overlaps);
}


=head2 get_Contig

 Title   : get_Contig
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=cut

sub get_Contig {
   my ($self,$contigid) = @_;

   if( $contigid ne $self->id() ) {
       $self->throw("In an Acedb database, trying to get a contigid $contigid not on the clone. Indicates an error!");
   }

   my ($c) = $self->get_all_Contigs(); 
   return $c;
}


=head2 get_all_Genes

 Title   : get_all_Genes
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=cut

sub get_all_Genes {
   my ($self,@args) = @_;
   my (@genes);

   foreach my $contig ( $self->get_all_Contigs ) {
       push(@genes,$contig->get_all_Genes());
   }
   return @genes;
}


=head2 _dbobj

 Title   : _dbobj
 Usage   : $obj->_dbobj($newval)
 Function: 
 Example : 
 Returns : value of _dbobj
 Args    : newvalue (optional)

=cut

sub _dbobj {
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'_dbobj'} = $value;
    }
    return $obj->{'_dbobj'};

}

1;
