
#
# BioPerl module for DB/ContigI.pm
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DB::ContigI.pm - Abstract Interface for Contig

=head1 SYNOPSIS

This is the abstract definition of a Contig, along with 'decorator'
functions

=head1 DESCRIPTION

Describe the object here

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::EMBLLOAD::Contig;
use vars qw(@ISA);
use strict;
use Bio::Root::RootI;

@ISA = qw(Bio::Root::RootI Bio::EnsEMBL::DB::ContigI);
use Bio::EnsEMBL::EMBLLOAD::Obj;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::DBEntry;


sub new {
    my($class,@args) = @_;
    
    my $self = {};
    bless $self,$class;

    my ($annseq,$id)=$self->_rearrange([qw(ANNSEQ)],@args);

    $self->_get_Seq($annseq);

    # HACK by th, for ensembl100:
    # inherit id from clone NOT annseq
    $id="$id.00001";
    $self->{'_id'} = $self->id($id);

    return $self; 
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

   my ($self,$value) = @_;
   if($value){
       $self->{'_id'}=$value;
   }
   return $self->{'_id'};

}


=head2 internal_id

 Title   : internal_id
 Usage   : $obj->internal_id($newval)
 Function: 
 Example : 
 Returns : value of internal_id
 Args    : newvalue (optional)


=cut

sub internal_id{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'internal_id'} = $value;
    }
    return $obj->{'internal_id'};

}



=head2 seq

 Title   : seq
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub seq{
   my ($self,@args) = @_;

   return $self->primary_seq->seq;
}

=head2 primary_seq

 Title   : primary_seq
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub primary_seq{
   my ($self,@args) = @_;

   my $seq=$self->_get_Seq;
   return $seq;
   
}



=head2 get_all_SeqFeatures

 Title   : get_all_SeqFeatures
 Usage   : foreach my $sf ( $contig->get_all_SeqFeatures ) 
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_SeqFeatures {
    my ($self) = @_;
    return ();
}



=head2 get_all_SimilarityFeatures

 Title   : get_all_SimilarityFeatures
 Usage   : foreach my $sf ( $contig->get_all_SimilarityFeatures ) 
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_SimilarityFeatures{
   my ($self) = @_;

   return ();
}




=head2 get_all_Genes

 Title   : get_all_Genes
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

Transcripts are not automatically grouped based on labels in EMBL
files, since these could create illegal ensembl entries, if entries to
not share exons.  If entry with same gene name, but without overlapping
exons exists, gene name is modified to make it unique.

=cut

sub get_all_Genes {
	
    my ($self)=@_;

    my $id = $self->id;
    my $exoncounter = 1;
    my $transcounter = 1;
    my $time = time();

    # first loop over all mRNA features, then loop over all CDS
    # features if CDS fits into an existing transcript, it becomes a
    # translation of that transcript.  If a transcript fits into an
    # existing gene, it becomes part of that gene.

    # exons hash is for fast lookup of exact duplicate exons (coordinates)
    # genes and transcripts hash are for name space clash checks
    my %genes;
    my %transcripts;
    my %exons;

    # WARN: merging genes across clone boundaries not done
    
    # turn on/off mRNA parsing
    my $parse_mrna=0;
    if($parse_mrna){
	foreach my $ft ( $self->_get_Seq->top_SeqFeatures ) {
	    my $ptag = $ft->primary_tag;
	    if($ptag eq 'mRNA') {
		$self->_build_gene($ft,$ptag,$id,\$exoncounter,\$transcounter,$time,
				   \%genes,\%transcripts,\%exons);
	    }
	}
    }
    foreach my $ft ( $self->_get_Seq->top_SeqFeatures ) {
	my $ptag=$ft->primary_tag;
	if($ptag eq 'CDS') {
	    $self->_build_gene($ft,$ptag,$id,\$exoncounter,\$transcounter,$time,
			       \%genes,\%transcripts,\%exons);
	}
    }

    return (values %genes);
}

sub _build_gene{

    my($self,$ft,$ptag,$id,$rexoncounter,$rtranscounter,$time,
       $rhgenes,$rhtranscripts,$rhexons)=@_;

    # get clone_id (accession) from contig id
    my $clone_id=$id;
    $clone_id=~s/\.\d+$//;


    # need to extract and process text from a number of tags

    # this could be accession.number or HUGO -> destined for dbentry table if HUGO
    my $gene_tag;
    if($ft->has_tag('gene')){
	($gene_tag)=$ft->each_tag_value('gene');
    }	

    # contains transcript name and discription.
    # extract transcript name (expect to be accession.number or accession.number.x)
    # and gene name (accession.number)
    my($product_tag,$description,$gene_id,$hack);
    if($ft->has_tag('product')){
	($product_tag)=$ft->each_tag_value('product');
	if($product_tag=~/^((\w+\.\d+)\S*)\s+\((.*)\)/){
	    $product_tag=$1;
	    $gene_id=$2;
	    $description=$3;
	}elsif($product_tag=~/^((\w+\.\d+)\S+)$/){
	    # no description
	    $product_tag=$1;
	    $gene_id=$2;
	}else{
	    # cannot be parsed
	    print "productID could not be extracted: \"$product_tag\"\n";
	    $description=$product_tag;
	    $hack=1;
	}
    }else{
	print "product tag missing\n";
	$hack=1;
    }
    if($hack){
	# HACK - not going to group transcripts into genes in this case
	$product_tag="$clone_id.$$rtranscounter";
	$gene_id=$product_tag;
	$$rtranscounter++;
    }

    # set gene_tag to gene_id if not set
    if(!$gene_tag){
	$gene_tag=$gene_id;
    }

    # destined for dbentry tables
    my $evidence_tag;
    if($ft->has_tag('evidence')){
	($evidence_tag)=$ft->each_tag_value('evidence');
    }
    my @dbentry;
    if($ft->has_tag('db_xref')){
	foreach my $db_xref ($ft->each_tag_value('db_xref')){
	    if($db_xref=~/^(\w+):(\w+)$/){
		my $id=$2;
		my $dbid;
		if($1 eq 'SPTREMBL'){
		    $dbid='SPTREMBL';
		}elsif($1 eq 'SWISS-PROT'){
		    $dbid='SP';
		}
		if($dbid){
		    my $dbentry=Bio::EnsEMBL::DBEntry->new(
							   -primary_id=>$id,
							   -display_id=>$id,
							   -version=>1,
							   -release=>1,
							   -dbname=>$dbid,
							   );
		    push(@dbentry,$dbentry);
		    print "Created DBENTRY: $dbid->$id\n";
		}else{
		    print "UNRECOGNISED: $db_xref\n";
		}
	    }
	}
    }
    if($ft->has_tag('protein_id')){
	my $dbid='protein_id';
	my($id)=($ft->each_tag_value($dbid));
	my $dbentry=Bio::EnsEMBL::DBEntry->new(
					       -primary_id=>$id,
					       -display_id=>$id,
					       -version=>1,
					       -release=>1,
					       -dbname=>$dbid,
					       );
	push(@dbentry,$dbentry);
	print "Created DBENTRY: $dbid->$id\n";
    }
    # accession is another valid dbentry for this transcript
    # add description in here too
    {
	my $dbid='EMBL';
	my $id=$clone_id;
	my $dbentry=Bio::EnsEMBL::DBEntry->new(
					       -primary_id=>$id,
					       -display_id=>$id,
					       -version=>1,
					       -release=>1,
					       -dbname=>$dbid,
					       );
	push(@dbentry,$dbentry);
	print "Created DBENTRY: $dbid->$id\n";
    }

    # exons - now are locations: same single, multi mess as before...
    my $loc=$ft->location;
    my @exon_features;
    if($loc->isa('Bio::Location::Split')){
	print "split loc\n";
	@exon_features=($loc->sub_Location);
    }elsif($loc->isa('Bio::Location::Simple')){
	print "simple loc\n";
	@exon_features=($loc);
    }else{
	$self->throw("unknown location type");
    }

    # in case of transcripts, create any new exons. In case of CDSs,
    # only create exons if cannot found exact fit of CDS inside
    # existing transcript (i.e. where there is no mRNA feature)

    # phases are meaningless except for 'translated bits', so they
    # are added when CDSs are parsed

    my $gene;
    my $trans;

    my $first;
    my $first_start;
    my $last;
    my $last_end;

    if($ptag!~/^mRNA/){
	# if CDS, see if we can identify a fitting transcript
    }
    if(!$trans){
	my $flag_existing_exon;

	# create new transcript, checking for name clash with existing one
	$trans=Bio::EnsEMBL::Transcript->new();
	# get unique id
	if($$rhtranscripts{$product_tag}){
	    $self->throw("transcriptID not unique: \"$product_tag\"");
	}
	$trans->id($product_tag);
	$$rhtranscripts{$product_tag}=$trans;
	$trans->version(1);

	# create exons, avoiding duplicates
	foreach my $sub (@exon_features){

	    my $exon;
	    my $st=$sub->start;
	    my $ed=$sub->end;

	    # if same start/end exists, get it
	    $exon=$$rhexons{"$st:$ed"};

	    # if exon exists, must be a gene associated with it
	    if($exon && !$gene){
		$flag_existing_exon=1;
	    }

	    # if new exon, create it
	    if(!$exon){
		$exon=Bio::EnsEMBL::Exon->new();
		$$rhexons{"$st:$ed"}=$exon;
		$exon->start($st);
		$exon->end($ed);
		$exon->strand($sub->strand);
		$exon->contig_id($self->id);
		$exon->seqname($self->id);
		$exon->version(1);
		$exon->created($time);
		$exon->modified($time);
		my $exon_id=$id.".exon.".$$rexoncounter++;
		$exon->id($exon_id);
		print STDERR "created exon $exon_id [".$trans->id."]\n";
	    }else{
		print STDERR "reused exon ".$exon->id." [".$trans->id."]\n";
	    }
	    $trans->add_Exon($exon);
	}

	# add dbentry objects to transcripts
	foreach my $dbentry (@dbentry){
	    $trans->add_DBLink($dbentry);
	}

	# add this transcript to gene, or create new gene
	if(!$flag_existing_exon){
	    if($$rhgenes{$gene_id}){
		$self->throw("gene $gene_id exists, but no exon overlap");
	    }
	    $gene=Bio::EnsEMBL::Gene->new();
	    $gene->id($gene_id);
	    $gene->version(1);
	    $$rhgenes{$gene_id}=$gene;
	    # add type tag to gene
	    if( $ft->has_tag('pseudo') ) {
		$gene->type('pseudo');
	    } else {
		$gene->type('standard');
	    }
	}else{
	    $gene=$$rhgenes{$gene_id};
	    if(!$gene){
		$self->throw("gene $gene_id does not exist, but exon overlap");
	    }
	}
	$gene->add_Transcript($trans);
    }


    # add Translation if CDS
    if($ptag=~/^CDS/){

	# if range of CDS in transcript not defined, must be entire
	# transcript (i.e. no mRNA record, or reading turned off)
	if(!$first){
	    my @exons = $trans->each_Exon;
	    $first = shift @exons;
	    if( $#exons == -1 ) {
		$last = $first;
	    } else {
		$last = pop @exons;
	    }
	    $first_start=1;
	    $last_end=$last->length;

	    # HACK BELOW
	    my $phase = 0;
	    foreach my $exon ($trans->each_Exon){
		$exon->phase($phase);
		$phase = $exon->end_phase();
	    }

	}

	# add phase for exons in range
	# HACK ABOVE

	my $tranl = Bio::EnsEMBL::Translation->new();
	$tranl->id($trans->id.".transl");
	$tranl->start_exon_id($first->id);
	$tranl->end_exon_id($last->id);
	$tranl->start($first_start);
	$tranl->end($last_end);
	$tranl->version(1);
	$trans->translation($tranl);
    }

}

=head2 length

 Title   : length
 Usage   : 
 Function: Provides the length of the contig
 Example :
 Returns : 
 Args    :


=cut

sub length {
   my ($self,@args) = @_;
   my $length=$self->_get_Seq->length;
   return $length;

}



sub _get_Seq {
    my ($self,$value) = @_;
    if (defined $value){$self->{'annseq'}=$value;}
    return $self->{'annseq'};
}







=head2 offset

 Title   : offset
 Usage   : $offset = $contig->offset()
 Function: Provides the offset of the contig in the clone
         : somehow. 1 means it is the first contig
 Example :
 Returns : 
 Args    :


=cut



sub offset{
   my ($self,@args) = @_;

   my $offset=2;
   #$self->throw("Object did not provide the offset method on Contig interface!");
   return $offset;
}



=head2 created

 Title   : created
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub created{
   my ($self) = @_;
   #$self->throw("Class [$self] has not implemented the created method");
   my $created=4;
   return $created;

}




=head2 seq_date

 Title   : seq_date
 Usage   : $contig->seq_date()
 Function: Gives the unix time value of the dna table created datetime field, which indicates
           the original time of the dna sequence data
 Example : $contig->seq_date()
 Returns : unix time
 Args    : none


=cut

sub seq_date{
    my ($self) = @_;

    my $seq_date=22;
#    $self->throw("Object did not provide the seq_date method on Contig interface!");
    return $seq_date;
}


=head2 embl_offset

 Title   : embl_offset
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub embl_offset{
   my ($self,@args) = @_;

   return 1;
}

=head2 embl_order

 Title   : embl_order
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub embl_order{
   my ($self,@args) = @_;

   return 1;
}


=head2 orientation

 Title   : orientation
 Usage   : 
 Function: Provides the orientation of the contig in the clone.
 Example :
 Returns : 
 Args    :


=cut

sub orientation{
   my ($self,@args) = @_;

   #$self->throw("Object did not provide the orientation method on Contig interface!");

   my $orientation=1;
   return $orientation;

}

=head2 order

 Title   : order
 Usage   : $obj->order($newval)
 Function: 
 Returns : value of order
 Args    : newvalue (optional)


=cut

sub order{
    my ($self,@args) = @_;

   # $self->throw("Object did not provide the order method on Contig interface!");
    my $order=1;
    return $order;

}





1;
