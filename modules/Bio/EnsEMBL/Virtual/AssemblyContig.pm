package Bio::EnsEMBL::Virtual::AssemblyContig;

use vars qw(@ISA);
use strict;

use Bio::Root::RootI;
use Bio::SeqFeatureI;

@ISA = qw(Bio::SeqFeatureI Bio::Root::RootI);

sub new {
    my ($class,@args) = @_;
    
    my $self = {};
    bless $self,$class;

    my ($dbID,$display_id,$chr_name,$chr_start,$chr_end,$start,$end,$orientation,$assembly_type) = 
	$self->_rearrange([qw( DBID
			       DISPLAY_ID
			       CHR_NAME
			       CHR_START
			       CHR_END
			       START
			       END
			       ORIENTATION
			       ASSEMBLY_TYPE
			       )],@args);

    if( 
	!defined $display_id ||
	!defined $chr_name ||
	!defined $chr_start ||
	!defined $chr_end ||
	!defined $start ||
	!defined $end ||
	!defined $orientation ||
	!defined $assembly_type ) {
      $self->throw("Did not pass all arguments into AssemblyContig");
    }

    $self->dbID($dbID);
    $self->display_id($display_id);
    $self->chr_name($chr_name);
    $self->chr_start($chr_start);
    $self->chr_end($chr_end);
    $self->start($start);
    $self->end($end);
    $self->orientation($orientation);
    $self->assembly_type($assembly_type);

    if ($start > $end) {
      $self->throw("Contig start must be less than end. start: $start > end: $end");
    }
    return $self;
}

=head2 display_id

 Title   : display_id
 Usage   : $obj->display_id($newval)
 Function: 
 Example : 
 Returns : value of contig display iud
 Args    : newvalue (optional)


=cut

sub display_id {
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'display_id'} = $value;
    }
    return $obj->{'display_id'};

}


=head2 dbID

 Title   : dbIDd
 Usage   : $obj->id($newval)
 Function: 
 Example : 
 Returns : value of contig internal db id
 Args    : newvalue (optional)


=cut

sub dbID {
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'id'} = $value;
    }
    return $obj->{'id'};

}

=head2 start

 Title   : start
 Usage   : $obj->start($newval)
 Function: 
 Example : 
 Returns : value of contig start
 Args    : newvalue (optional)


=cut

sub start {
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'start'} = $value;
    }
    return $obj->{'start'};

}


=head2 end

 Title   : end
 Usage   : $obj->end($newval)
 Function: 
 Example : 
 Returns : value of contig end
 Args    : newvalue (optional)


=cut

sub end {
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'end'} = $value;
    }
    return $obj->{'end'};

}


=head2 chr_start

 Title   : chr_start
 Usage   : $obj->chr_start($newval)
 Function: 
 Example : 
 Returns : value of chr start
 Args    : newvalue (optional)


=cut

sub chr_start {
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'chr_start'} = $value;
    }
    return $obj->{'chr_start'};

}


=head2 chr_end

 Title   : chr_end
 Usage   : $obj->chr_end($newval)
 Function: 
 Example : 
 Returns : value of chr end
 Args    : newvalue (optional)


=cut

sub chr_end {
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'chr_end'} = $value;
    }
    return $obj->{'chr_end'};

}


=head2 chr_name

 Title   : chr_name
 Usage   : $obj->chr_name($newbal)
 Function: 
 Example : 
 Returns : value of chr_name
 Args    : newvalue (optional)


=cut

sub chr_name {
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'chr_name'} = $value;
    }
    return $obj->{'chr_name'};

}


=head2 orientation

 Title   : orientation
 Usage   : $obj->orientation($newval)
 Function: 
 Example : 
 Returns : value of orientation
 Args    : newvalue (optional)


=cut

sub orientation {
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'orientation'} = $value;
    }
    return $obj->{'orientation'};

}

=head2 assembly_type

 Title   : assembly_type
 Usage   : $obj->assembly_type($newval)
 Function: 
 Example : 
 Returns : value of assembly_type
 Args    : newvalue (optional)


=cut

sub assembly_type {
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'assembly_type'} = $value;
    }
    return $obj->{'assembly_type'};

}

=head2 SeqFeature methods

These methods are here for seqfeature I compliance

=cut


=head2 strand

 Title   : strand
 Usage   : this is for seqfeatureI compliance
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub strand{
   my ($self) = @_;

   return $self->orientation;
}

=head2 source_tag

 Title   : source_tag
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub source_tag {
   my ($self,@args) = @_;

   return 'ensembl';
}
    

=head2 primary_tag

 Title   : primary_tag
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub primary_tag{
   my ($self) = @_;

   return 'fragment';
}

=head2 has_tag_value

 Title   : has_tag_value
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub has_tag_value{
   my ($self,@args) = @_;

   return 0;
}

=head2 each_tag_value

 Title   : each_tag_value
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub each_tag_value{
   my ($self,@args) = @_;

   return ();

}


=head2 all_tags

 Title   : all_tags
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub all_tags{
   my ($self) = @_;

   return ();
}


=head2 to_FTHelper

 Title   : to_FTHelper
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub to_FTHelper{
   my ($self,@args) = @_;

   my $ft = Bio::SeqIO::FTHelper->new();

   my $loc = $self->start."..".$self->end;
   $ft->loc($loc);
   $ft->key('misc');
   $ft->add_field('note','Component DNA fragment');
   # grrr. So frustrating. Have to get clone version
   my $clone = $self->contig->dbobj->get_Clone($self->contig->cloneid);
   $ft->add_field('note',"accession=".$self->contig->cloneid.".".$clone->embl_version);
   $ft->add_field('note',"start=".($self->rawcontig_start+$self->contig->embl_offset));
   $ft->add_field('note',"end=".($self->rawcontig_end+$self->contig->embl_offset));
   $ft->add_field('note',"orientation=".$self->orientation);
   
   return $ft;
}


1;

