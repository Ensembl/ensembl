
# BioPerl module for Bio::EnsEMBL::TimDB::Clone
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::SlimTim::Clone - Perl wrapper over Tim's directories for Clones

=head1 SYNOPSIS


=head1 DESCRIPTION

Describe the object here

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let's begin the code:

package Bio::EnsEMBL::SlimTim::Clone;
use vars qw($AUTOLOAD @ISA);
use strict;
use Bio::EnsEMBL::DB::CloneI;
use Bio::EnsEMBL::SlimTim::Contig;
use Bio::EnsEMBL::Analysis::ensConf qw(
				       EXON_ID_SUBSCRIPT
				       );
use Bio::EnsEMBL::Analysis::LegacyParser;
use Bio::SeqIO;

use Fcntl qw( O_RDONLY );
use FileHandle;

# Object preamble - inheriets from Bio::Root::Object
use Bio::Root::Object;

@ISA = qw(Bio::Root::Object Bio::EnsEMBL::DB::CloneI);
# new() is inherited from Bio::Root::Object



# _initialize is where the heavy stuff will happen when new is called
sub _initialize {
  my($self,@args) = @_;
  my $make = $self->SUPER::_initialize(@args);

  # set stuff in self from @args
  my ($dbobj,$id)=
      $self->_rearrange([qw(DBOBJ
			    ID
			    )],@args);

  $id      || $self->throw("Cannot make contig db object without id");
  $dbobj   || $self->throw("Cannot make contig db object without db object");

  $dbobj->isa('Bio::EnsEMBL::SlimTim::Obj') || 
      $self->throw("Cannot make contig db object with a $dbobj object");

  my $dir = $dbobj->_dir . "$id";
  if( ! -e $dir ) {
      $self->throw("Directory $dir does not exist - cannot load clones from it");
  }
  print STDERR "Using $dir\n";

  # we have to load the contigs first to get out clone 
  # information - so we need somewhere to put them.
  $self->{_contig_hash} = {};

  # the sequences are in one file (doh!) so we need to get them
  # out first and store them
  $self->{_sequence_hash} = {};

  $self->_dir($dir);
  $self->id          ($id);
  $self->fetch_sequences();
  $self->fetch_contigs();
  $self->fetch();

  return $make; # success - we hope!
}

=head2 fetch_sequences

 Title   : fetch_sequences
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub fetch_sequences{
   my ($self) = @_;
   my $file = $self->_dir . "/". $self->id.".seq"; 

   if( !-e  $file ) {
       $self->throw("$file does not exist!");
   }

   my $seqio = Bio::SeqIO->new(-file => $file, -format => 'Fasta' );
   my $seq;
   while( $seq = $seqio->next_seq ) {
       $self->{_sequence_hash}->{$seq->id} = $seq;
   }

}

=head2 fetch_contigs

 Title   : fetch_contigs
 Usage   :
 Function:
 Example :
 Returns :  
 Args    :


=cut

sub fetch_contigs{
   my ($self,@args) = @_;

   # fetch sequences is called before fetch contigs
   foreach my $id ( keys %{$self->{_sequence_hash}} ) {
       my $contig = Bio::EnsEMBL::SlimTim::Contig->new( -id => $id,
						-disk_id => $id,
						-dir => $self->_dir,
						-seq => $self->{_sequence_hash}->{$id}
						);
       $contig->seq_date(0);
       $self->{_contig_hash}->{$id} = $contig;
   }

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
    my @contigs = $self->get_all_Contigs();

    my $contig = shift @contigs;

    # Fill clone object
    my $accsv = $contig->clone_accession;
    if ( !($accsv =~ /(\S+)\.(\S+)/)) {
	$self->throw("Bad accsv $accsv");
    }

    my $acc = $1;
    my $sv = $2;

    $self->embl_version($sv);
    $self->embl_id     ($acc);
    #$self->chromosome  ($chr);
    #$self->species     ($species);
    my $clonephase = $contig->clone_phase;
    if( $clonephase =~ /HTGS_PHASE(\d)/ ) {
	$clonephase = $1;
    } else {
	$clonephase = 4;
    }
    $self->version(100);
    $self->created(0);
    $self->modified(0);
    $self->htg_phase   ($clonephase);

    return $self;
}

=head2 chromosome

 Title   : chromosome
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub chromosome{
   my ($self,@args) = @_;
   return Bio::EnsEMBL::Species->chromosome_by_name('human','unknown');
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

   return values %{$self->{_contig_hash}};
}


=head2 _parse_exon

 Title   : _parse_exon
    (cut down version of routine in Legacyparser, using MappedExon.
     Missing pepide, phase, etc.]
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=cut

sub _parse_exon{
    my ($self,$file,$contig_id) = @_;

    my @exons;

    my $fh = new FileHandle;
    $fh->open($file) || 
	$self->throw("Could not open file [", $file, "]");

    my $id=$self->id;
    my $disk_id=$self->disk_id;
    
    my $is = $fh->input_record_separator('>');
    my $dis = <$fh>; # skip first record (dull!)
    while( <$fh> ) {
	if ( /^(\S+)\s(\S+\.\S+):(\d+)-(\d+):*(\d*)\s.*(\d{4}-\d\d-\d\d_[\d:]+)\s+(\d{4}-\d\d-\d\d_[\d:]+).*\n(.*)/  ) {
	    my $e = $1;
	    my $contigid = $2;
	    my $start = $3;
	    my $end = $4;
	    my $phase = $5;
	    my $created = &Bio::EnsEMBL::Analysis::LegacyParser::acetime($6,1);
	    my $modified = &Bio::EnsEMBL::Analysis::LegacyParser::acetime($7,1);
	    my $pep = $8;
	    
	    $pep =~ s/\s+//g;
	    
	    my($eid,$ever);
	    if($e=~/^$EXON_ID_SUBSCRIPT(\d+)\.(\d+)/){
		$eid=$EXON_ID_SUBSCRIPT.$1;
		$ever=$2;
	    }else{
		$self->throw("Exon identifier could not be parsed: $eid");
	    }

	    # convert from disk_id to whatever internal clone_id is
	    $contigid=~s/^$disk_id/$id/;
	    my $cloneid=$contigid;
	    $cloneid=~s/\.\d+$//;

	    # if calling by contig, skip those that are wrong
	    next if($contig_id && $contigid ne $contig_id);

	    my $exon = Bio::EnsEMBL::MappedExon->new();
	    $exon->id($eid);
	    $exon->version($ever);

	    $exon->contig_id($contigid);
	    $exon->clone_id($cloneid);
	    
	    if( $end < $start ) {
		my $s = $end;
		$end = $start;
		$start = $s;
		$exon->strand(-1);
	    } else {
		$exon->strand(1);
	    }

	    $exon->start($start);
	    $exon->end($end);
	    
	    $exon->created($created);
	    $exon->modified($modified);
	    push(@exons,$exon);
	} else {
	    chomp;
	    $self->throw("Yikes. Line with > but not parsed! [$_]");
	}
    }
    $fh->input_record_separator($is);
    return @exons;
}

=head2 get_all_ContigOverlaps

 Title   : get_all_ContigOverlaps
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_ContigOverlaps{
   my ($self,@args) = @_;

   return ();
}

sub get_all_Genes {
    return ();
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
    my $dnafile=$self->_dir . "/" . $self->id . ".seq";
    my $dnafiledate=(stat($dnafile))[9];
    return $dnafiledate;
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
   $obj->warn(" *** WARNING *** Ewan has not figured out how to get clone version yet");
   if( @_ ) {
      my $value = shift;
      $obj->{'version'} = $value;
    }
    return $obj->{'version'};

}

=head2 created

 Title   : created
 Usage   : $obj->created($newval)
 Function: 
 Returns : value of created
 Args    : newvalue (optional)


=cut

sub created{
   my $obj = shift;
   $obj->warn("*** WARNING *** Ewan has not figured out how to get created value for a clone");
   if( @_ ) {
      my $value = shift;
      $obj->{'created'} = $value;
    }
    return $obj->{'created'};

}

=head2 modified

 Title   : modified
 Usage   : $obj->modified($newval)
 Function: 
 Returns : value of modified
 Args    : newvalue (optional)


=cut

sub modified{
   my $obj = shift;

   $obj->warn("*** WARNING *** Ewan has not figured out how to get created value for a clone");

   if( @_ ) {
      my $value = shift;
      $obj->{'modified'} = $value;
    }
    return $obj->{'modified'};

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

sub freeze {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_freeze} = $arg;
    }

    return $self->{_freeze};
}

=head2 _dir

 Title   : _dir
 Usage   : $obj->_dir($newval)
 Function: 
 Returns : value of _dir
 Args    : newvalue (optional)


=cut

sub _dir{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_dir'} = $value;
    }
    return $obj->{'_dir'};

}

1;




