
# BioPerl module for Bio::EnsEMBL::TimDB::Clone
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::TimDB::Clone - Perl wrapper over Tim's directories for Clones

=head1 SYNOPSIS

    $clone = Bio::EnsEMBL::TimDB::Clone->new();
 
    $clone->add_Contig($contig);
    

=head1 DESCRIPTION

Describe the object here

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let's begin the code:

package Bio::EnsEMBL::TimDB::Clone;
use vars qw($AUTOLOAD @ISA);
use strict;
use Bio::EnsEMBL::DB::CloneI;
use Bio::EnsEMBL::TimDB::Contig;
use Bio::EnsEMBL::ContigOverlap;

use Bio::SeqIO;

use NDBM_File;
use Fcntl qw( O_RDONLY );

# Object preamble - inheriets from Bio::Root::Object
use Bio::Root::Object;

@ISA = qw(Bio::Root::Object Bio::EnsEMBL::DB::CloneI);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called
sub _initialize {
  my($self,@args) = @_;
  my $make = $self->SUPER::_initialize(@args);

  # set stuff in self from @args
  my ($dbobj,$id,$cgp,$disk_id,$sv,$emblid,$htgsp,$byacc,$chr,$species)=
      $self->_rearrange([qw(DBOBJ
			    ID
			    )],@args);

  $id      || $self->throw("Cannot make contig db object without id");
  $dbobj   || $self->throw("Cannot make contig db object without db object");

  $dbobj->isa('Bio::EnsEMBL::TimDB::Obj') || 
      $self->throw("Cannot make contig db object with a $dbobj object");

  $self->_dbobj      ($dbobj);  
  $self->id          ($id);
  
  return $make; # success - we hope!
}

=head2 fetch

 Title   : fetch
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=cut

sub fetch {
    my ($self) = @_;

    my $id=$self->id;
    
    # translate incoming id to ensembl_id, taking into account nacc flag
    my($disk_id,$cgp,$sv,$emblid,$htgsp,$chr,$species);
    ($id,$disk_id,$cgp,$sv,$emblid,$htgsp,$chr,$species)=
	$self->_dbobj->get_id_acc($id);
    if($id eq 'unk'){
	$self->throw("Cannot get accession for $disk_id");
    }
    
    # if its already been created, get it from the hash
    if($self->_dbobj->{'_clone_array'}->{$id}){
	return $self->_dbobj->{'_clone_array'}->{$id};
    }
    
    # else, have to build it
    
    # can only build it, if it was 'loaded' in the initial call to timdb
    # (i.e. that it is in the active list)
    unless($self->_dbobj->{'_active_clones'}->{$id}){
	$self->throw("$id fetched - not loaded in original object call - locked?");
    }

    # test if clone is not locked (for safety); don't check for valid SV's
    # (probably overkill, as locking at Obj.pm level)
    my($flock,$fsv,$facc)=$self->_dbobj->_check_clone_entry($disk_id);
    if($flock){
	$self->throw("$id is locked by TimDB");
    }
    if($facc){
	$self->throw("$id does not have an accession number");
    }

    my $byacc= $self->{'_byacc'};
    
    # Fill clone object
    $self->disk_id     ($disk_id);
    $self->embl_version($sv);
    $self->embl_id     ($emblid);
    $self->chromosome  ($chr);
    $self->species     ($species);
    $self->htg_phase   ($htgsp);
    $self->byacc       ($byacc);

    # construct and test the directory of the clone
    # fast (direct)
    my $cgp_dir         = $self->_dbobj->{'_unfin_data_root_cgp'}->{$cgp};
    my $clone_dir       = "$cgp_dir/data/$disk_id";
    my $contig_dbm_file = "$cgp_dir/unfinished_ana.dbm";
    
    $self->clone_dir     ($clone_dir,$disk_id);
    $self->build_contigs ($contig_dbm_file);
    
    
    # save it to hash
    $self->_dbobj->{'_clone_array'}->{$id}=$self;
    
    return $self;
}

sub build_contigs {
    my ($self,$contig_dbm_file) = @_;

    # build list of contigs for clone
    # (methods get_all_Contigs and get_Contig look at this list of objects, 
    # rather than build it)
    
    my %unfin_contig;

    unless(tie(%unfin_contig,'NDBM_File',$contig_dbm_file,O_RDONLY,0644)){
	$self->throw("Error opening contig dbm file");
    }

    my ($key,$val);

    my $spacing = $Bio::EnsEMBL::DB::CloneI::CONTIG_SPACING;
    my $disk_id = $self->disk_id;
    my $id      = $self->id;

    my @contigs;

    while (($key,$val) = each %unfin_contig) {

#	print(STDERR "[$key][$val]\n");

	if($key=~/^$disk_id/){
	    
	    my($len,$checksum,$embl_offset,$embl_order) = split(/,/,$val);
	  
	    # all of these values should be positive, non zero
	    $self->throw("Error: invalid length [$len] for contig $key") 		if !$len         || $len<1;
	    $self->throw("Error: invalid checksum [$checksum] for contig $key") 	if !$checksum    || $checksum<1;
	    $self->throw("Error: invalid embl_order [$embl_order] for contig $key")	if !$embl_order  || $embl_order<1;
	    $self->throw("Error: invalid embl_order [$embl_offset] for contig $key")	if !$embl_offset  || $embl_offset<1;
	    
	    my $disk_key = $key;
	    $key =~ s/^$disk_id/$id/;

	    print STDERR "Attempting to retrieve contig with $disk_key [$key]\n";
	    
	    my $tmpcontig = new Bio::EnsEMBL::TimDB::Contig( -dbobj => $self->_dbobj,
							     -id    => $key);

	    $tmpcontig->length     ($len);
	    $tmpcontig->disk_id    ($disk_key);
	    $tmpcontig->_clone_dir ($self->clone_dir);
	    $tmpcontig->chromosome ($self->chromosome,$self->species);
	    $tmpcontig->checksum   ($checksum);
	    $tmpcontig->embl_order($embl_order);
	    $tmpcontig->embl_offset($embl_offset);

	    push(@contigs,$tmpcontig);

	    $self->add_Contig($tmpcontig);
	}

    }

    $self->_make_ContigOverlaps;
    
    print STDERR scalar($self->get_all_Contigs) . " contigs found in clone\n";
}

sub _make_ContigOverlaps {
    my ($self) = @_;

    my $clone_order = $self->_dbobj->{'_contig_order_hash'}->{$self->id};
    my $spacing     = $Bio::EnsEMBL::DB::CloneI::CONTIG_SPACING;
    
    return unless defined($clone_order);
    
    my @pieces = split(/:/,$clone_order);

    # Each contig is in the $clone_order string separated by either a colon or
    # a semi-colon.  Contigs separated by semi-colons are ordered. The ones
    # separated by colons are unordered

    my $finalcontig;
    my $finalorient;
    my @unordered;

    foreach my $piece (@pieces) {
	# Do we have an ordered piece (otherwise put into @unordered )
	if ($piece =~ /;/) {

	    my @joins      = split(/;/,$piece);
	    my $numcontigs = scalar(@joins);
	    
	    # First of all attach to previous piece

	    if (defined($finalcontig)) {
		my $newtype;
		my $positiona;
		my $positionb;

		my ($contigid,$fr) = ($joins[0] =~ /(.*)\.([FR])$/);
		my $contiga = $self->get_Contig($contigid);

		$self->throw("No contig [$contigid] existis in clone") unless defined($contiga);

		if ($finalorient == 1) {
		    $positiona = $finalcontig->length;
		    
		    if ($fr eq 'R') {
			$newtype = 'right2right';
			$positionb = $contiga->length;
		    } elsif ($fr eq 'F') {
			$newtype = 'right2left';
			$positionb = 1;
		    } else {
			$self->throw("Wrong frame [$fr]");
		    }
		} elsif ($finalorient == -1) {
		    $positiona = 1;
		    if ($fr eq 'R') {
			$newtype = 'left2right';
			$positionb = $contiga->length;
		    } elsif ($fr eq 'F') {
			$newtype = 'left2left';
			$positionb = 1;
		    } else {
			$self->throw("Wrong frame [$fr]");
		    }
		} else {
		    $self->throw("Wrong finalorient [$finalorient]");
		}
		print(STDERR "Joinging " .$finalcontig->id . "\t" . $contiga->id . "\n");
		my $tmpoverlap = new Bio::EnsEMBL::ContigOverlap(-contiga   => $finalcontig,
								 -contigb   => $contiga,
								 -positiona => $positiona,
								 -positionb => $positionb,
								 -source    => 'UNORDERED',
								 -distance  => $spacing,
								 -overlap_type => $newtype);
		
		$self->add_ContigOverlap($tmpoverlap);
	    }
	    
	    # Now join the ordered contigs together
	    for (my $i = 0; $i < ($numcontigs-1); $i++) {
		# Split the pieces into id and orientation (D89885.00001.F)
		my ($contig1,$fr1)= ($joins[$i] =~ /(.*)\.([FR])$/);
		my ($contig2,$fr2)= ($joins[$i+1] =~ /(.*)\.([FR])$/);

		my $contiga = $self->get_Contig($contig1) || $self->throw("No contig with id [$contig1]");
		my $contigb = $self->get_Contig($contig2) || $self->throw("No contig with id [$contig2]");

		my $type;
		my $positiona;
		my $positionb;

		if ($fr1 eq 'F' && $fr2 eq 'F') {
		    $type = 'right2left';
		    $positiona = $contiga->length;
		    $positionb = 1;

		} elsif ($fr1 eq 'F' && $fr2 eq 'R') {

		    $type = 'right2right';
		    $positiona = $contiga->length;
		    $positionb = $contigb->length;

		} elsif ($fr1 eq 'R' && $fr2 eq 'F') {
		    $type = 'left2left';
		    $positiona = 1;
		    $positionb = 1;

		} elsif ($fr1 eq 'R' && $fr2 eq 'R') {
		    $type = 'left2right';
		    $positiona = 1;
		    $positionb = $contigb->length;
		}

		my $overlap = new Bio::EnsEMBL::ContigOverlap(-contiga   => $contiga,
							      -contigb   => $contigb,
							      -positiona => $positiona,
							      -positionb => $positionb,
							      -source    => 'CLONE',
							      -distance  => $spacing,
							      -overlap_type => $type);
		
		$self->add_ContigOverlap($overlap);

		$finalcontig = $contigb;
		$finalorient = 1  if $type =~ /2left/;
		$finalorient = -1 if $type =~ /2right/;

	    }
	    
	    
	} else {
	    push(@unordered,$piece);
	}
    }

    # Now process the unorderedcontigs
    my @newcontigs;
    
    foreach my $contigid (@unordered) {
	$contigid =~ s/\.[FR]$//;
	my $contig = $self->get_Contig($contigid) || $self->throw("No contig with id [$contigid]");
	push(@newcontigs,$contig);
    }
    
    # Sort longest to shortest
    @newcontigs = sort { $b->length <=> $a->length } @newcontigs;
    
    my $numcontigs = scalar(@newcontigs);
    
    # Join the final ordered contig onto the first unordered contig
    if (defined($finalcontig) && $numcontigs >  0) {
	my $type;
	my $positiona;
	my $positionb;
	
	if ($finalorient == 1) {
	    $type      = 'right2left';
	    $positiona = $finalcontig->length;
	} elsif ($finalorient == -1) {
	    $type      = 'left2left';
	    $positiona = 1;
	} else {
	    $self->throw("Final orientation not set");
	}
	
	my $overlap = new Bio::EnsEMBL::ContigOverlap(-contiga   => $finalcontig,
						      -contigb   => $newcontigs[0],
						      -positiona => $positiona,
						      -positionb => 1,
						      -source    => 'UNORDERED',
						      -distance  => $spacing,
						      -overlap_type => $type);
	
	$self->add_ContigOverlap($overlap);
    }
    
    
    # Create overlaps between the unordered contigs
    for (my $i = 0; $i < $numcontigs-1; $i++) {
	my $overlap = new Bio::EnsEMBL::ContigOverlap(-contiga   => $newcontigs[$i],
						      -contigb   => $newcontigs[$i+1],
						      -positiona => $newcontigs[$i]->length,
						      -positionb => 1,
						      -source    => 'UNORDERED',
						      -distance  => $spacing,
						      -overlap_type => 'right2left');
	
	$self->add_ContigOverlap($overlap);
    }
}

						      
=head2 add_ContigOverlap

 Title   : add_ContigOverlap
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=cut

sub add_ContigOverlap {
    my ($self,$overlap) = @_;

    if (!defined($overlap) || !($overlap->isa("Bio::EnsEMBL::ContigOverlap"))) {
	$self->throw("Argument [$overlap] is not a Bio::EnsEMBL::ContigOverlap");
    }

    if (!defined($self->{_overlaps})) {
	$self->{_overlaps} = [];
    }

    push(@{$self->{_overlaps}},$overlap);
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

    if (!defined($self->{_overlaps})) {
	$self->{_overlaps} = [];
    }

    return @{$self->{_overlaps}};
}

=head2 get_all_Contigs

 Title   : get_all_Contigs
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=cut

sub get_all_Contigs {
   my ($self) = @_;

   if (!defined($self->{_contigs})) {
       $self->{_contigs} = [];
   }

   return @{$self->{_contigs}};
}

=pod

=head2 get_all_Genes

 Title   : get_all_Genes
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=cut


sub get_all_Genes{
    my ($self,$evidence) = @_;
    my %h;
    
    # map if not already mapped
    $self->_dbobj->map_etg unless $self->_dbobj->{'_mapped'};

    my @features;

    foreach my $contig ($self->get_all_Contigs) {

	if ($evidence eq 'evidence') {
	    push(@features,$contig->get_all_SimilarityFeatures);
	}

	foreach my $gene ($contig->get_all_Genes){
	    # read into a hash to make unique
	    $h{$gene->id()} = $gene;
	}
    }

    # Now attach the evidence if necessary
    if ($evidence eq 'evidence') {
	foreach my $gene (values %h) {
	    foreach my $exon ($gene->each_unique_Exon) {
		$exon ->find_supporting_evidence (\@features);
	    }
	}
    }
    # DEBUG
    print STDERR "Clone contains ".scalar(keys %h)." genes\n";
    return values %h;
}

#
# Seq method from the cloneI object now
#


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

=head2 disk_id

 Title   : disk_id
 Usage   : $obj->disk_id($newval)
  Function:
 Example :
 Returns : value of disk_id
 Args    : newvalue (optional)

=cut

sub disk_id {
    my ($obj,$value) = @_;
    if( defined $value) {
	$obj->{'_clone_disk_id'} = $value;
    }
    return $obj->{'_clone_disk_id'};
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

=head2 embl_id

 Title   : embl_id
 Usage   : this is the embl_id for this clone, to generate nice looking files
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub embl_id {
    my ($obj,$value) = @_;
    if( defined $value) {

	# FIXME
	# may be '' - what to do in this case?

	$obj->{'_clone_embl_id'} = $value;
    }
    return $obj->{'_clone_embl_id'};
}


=head2 sv

 Title   : sv
 Function: returns the version number (not the acc.version, just verision).
 Example :
 Returns : 
 Args    :


=cut

sub sv {
    my($self)=@_;
    $self->warn("DEPRECATING METHOD replace \$self->sv with \$self->embl_version");
    return $self->embl_version;
}


=head2 htg_phase

 Title   : htg_phase
 Usage   : this is the phase being 0,1,2,3,4 (4 being finished).
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub htg_phase {
    my ($obj,$value) = @_;
    if( defined $value) {
	if($value=~/^[01234]$/){
	    $obj->{'_clone_htgsp'} = $value;
	}else{
	#    $obj->throw("Invalid value for htg_phase $value");
	     $obj->warn("Invalid value for htg_phase $value. Storing undef for the moment");
	}
    }
    return $obj->{'_clone_htgsp'};
}


=head2 byacc

 Title   : byacc
 Usage   : $obj->byacc($newval)
 Function: 
 Returns : value of byacc
 Args    : newvalue (optional)


=cut

sub byacc{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'byacc'} = $value;
    }
    return $obj->{'byacc'};

}


=head2 compare_dna

 Title   : compare_dna
 Usage   : $obj->compare_dna($file)
 Function: 
 Returns : 1 if dna in $file is different (checksum comparision) to dna in clone object
 Args    : newvalue (optional)


=cut

sub compare_dna {
    my $self = shift;
    my $file = shift;

    print "Comparing with $file\n";

    my $seqio=Bio::SeqIO->new( '-format' => 'Fasta', -file => $file);

    my $fseq;
    my %contigs;
    my $eflag=0;

    my $id      = $self->id;
    my $disk_id = $self->disk_id;

    while ($fseq=$seqio->next_seq()){
	my $seqid = $fseq->id; 
	$seqid =~ s/^$disk_id/$id/;
	my $seqlen   = $fseq->seq_len;
	my $seq      = $fseq->seq;
	my $checksum = unpack("%32C*",$seq) % 32767;

	my $contig = $self->get_Contig($seqid);

	if (defined($contig)) {
	    my $checksum2 = $contig->checksum;
	    $contigs{$seqid} = 1;
	    
	    if ($checksum == $checksum2){
		print STDERR "Contig $seqid same in Ensembl: $checksum $checksum2\n";
	    } else {
		print STDERR "ERROR: Contig $seqid different in Ensembl: $checksum $checksum2\n";
		$eflag=1;
	    }
	} else {
	    print STDERR "ERROR: Contig $seqid not in Ensembl\n";
	    $eflag=1;
	}
    }
    
    foreach my $contig ($self->get_all_ContigIds){
	if(!$contigs{$contig}){
	    print STDERR "ERROR: Contig $contig only in Ensembl\n";
	}
    }
    return $eflag;
}


=head2 embl_version

 Title   : embl_version
 Usage   : $clone->embl_version()
 Function: Gives the value of the EMBL version, i.e. the sequence data version [SV]
 Example : $clone->embl_version()
 Returns : version number
 Args    : none


=cut

sub embl_version {
    my ($obj,$value) = @_;
    if (defined $value) {
        $value ||= -1;
	if ($value =~ /^-?\d+$/) {
	    $obj->{'_clone_sv'} = $value;
	} else {
	    $obj->throw("Invalid value for embl_version [SV] '$value'");
	}
    }
    return $obj->{'_clone_sv'};
}

=head2 seq_date

 Title   : seq_date
 Usage   : $clone->seq_date()
 Function: In TimDB there is one DNA file per clone, so just checks this date
 Example : $clone->seq_date()
 Returns : unix time
 Args    : none


=cut

sub seq_date {
    my ($self) = @_;
    my $dnafile=$self->clone_dir . "/" . $self->disk_id . ".seq";
    my $dnafiledate=(stat($dnafile))[9];
    return $dnafiledate;
}


=head2 version

 Title   : version
 Function: Schema translation
 Example :
 Returns : 
 Args    :


=cut

sub version {
    my ($self,@args) = @_;
    # this value is incremented each time a clone is unlocked after processing in TimDB
    return $self->_clone_status(6);
}


=head2 created

 Title   : created
 Usage   : $clone->created()
 Function: Gives the unix time value of the created datetime field, which indicates
           the first time this clone was put in ensembl
 Example : $clone->created()
 Returns : unix time
 Args    : none


=cut

sub created {
    my ($self) = @_;
    # this value is the time set when a clone is unlocked after being first created in TimDB
    return $self->_clone_status(5);
}


=head2 modified

 Title   : modified
 Usage   : $clone->modified()
 Function: Gives the unix time value of the modified datetime field, which indicates
           the last time this clone was modified in ensembl
 Example : $clone->modified()
 Returns : unix time
 Args    : none


=cut

sub modified{
    my ($self) = @_;
    # this value is the time set when a clone is unlocked in TimDB
    return $self->_clone_status(0);
}

sub _clone_status{
    my($self,$field)=@_;
    my $disk_id=$self->disk_id;
    my $val=$self->_dbobj->{'_clone_update_dbm'}->{$disk_id};
    if(!$val){
	my $id=$self->id;
	$self->throw("No Status entry for $disk_id [$id]");
    }
    my @fields=split(',',$val);
    return $fields[$field];
}

sub chromosome {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_chromosome} = $arg;
    }
    return $self->{_chromosome};
}

sub species {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_species} = $arg;
    }

    return $self->{_species};
}

sub clone_dir {
    my ($self,$arg,$disk_id) = @_;

    if (defined($arg)) {

	# Check the argument is a directory
	unless(-d $arg){
	    $self->throw("Cannot find directory for $disk_id");
	}

	# check for sequence file
	if (!-e "$arg/$disk_id.seq"){
	    $self->throw("Error: no sequence file for entry " . $self->id . " ($disk_id)");
	}

	$self->{_clone_dir} = $arg;
    } 
    
    return $self->{_clone_dir};
}


=head2 add_Contig

 Title   : add_Contig
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=cut

sub add_Contig {
    my ($self,$contig) = @_;
    
    if (!defined($contig) || !$contig->isa("Bio::EnsEMBL::TimDB::Contig")) {
	$self->throw("[$contig] is not a Bio::EnsEMBL::TimDB::Contig");
    }

    if (!defined($self->{_contigs})) {
	$self->{_contigs} = [];
    }

    $contig->validate();

    # check for gs file
    my $gsfile = $self->clone_dir . "/" . $contig->disk_id . ".gs";
    print(STDERR "Genscan file [$gsfile]\n");

    if(!-e $gsfile){
	$self->throw("Error: no gs file [$gsfile] for contig " . $contig->id);
    }

    push(@{$self->{_contigs}},$contig); 

    # Also add to the contig hash so we can retrieve it by name

    my $id = $contig->id;
    $self->{_contighash}{$id} = $contig;
    
}


=head2 get_Contig

 Title   : get_Contig
 Usage   :
 Function:
 Example :
 Returns : contig object
 Args    : contig_id

=cut

sub get_Contig {
    my ($self,$arg) = @_;

    $self->throw("No id input") unless defined($arg);

    return $self->{_contighash}{$arg};
}


=head2 get_all_ContigIds

 Title   : get_all_ContigIds
 Usage   :
 Function:
 Example :
 Returns : Array of strings
 Args    : none

=cut

sub get_all_ContigIds {
    my ($self) = @_;

    return (keys %{$self->{_contighash}});
}

1;



