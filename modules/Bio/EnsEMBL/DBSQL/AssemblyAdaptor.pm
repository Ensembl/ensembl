#
# Ensembl module for Bio::EnsEMBL::DBSQL::AssemblyAdaptor
#
# Cared for by James Smith <js5@sanger.ac.uk>
#
# Copyright James Smith
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::AssemblyAdaptor

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Database adaptor to provide access to Assembly objects

=head1 AUTHOR

James Smith

This modules is part of the Ensembl project http://www.ensembl.org

=head1 CONTACT

Email js5@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DBSQL::AssemblyAdaptor;
use vars qw(@ISA $AUTOLOAD);
use Bio::EnsEMBL::Assembly;
use strict;
use Bio::EnsEMBL::Utils::Eprof qw(eprof_start eprof_end eprof_dump);

# Object preamble - inherits from Bio::Root::RootI

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Eprof;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

# inherit new from BaseAdaptor


=head2 fetch_by_chromosome_start_end

 Title   : fetch_by_chromosome_start_end
 Usage   : $band_obj = $kary_adp->fetch_by_chromosome_position('chr1',10000, 2000);
 Function: Retrieves Assembly objects by chromosome ('chrN' notation)
           and its absolute "golden-path" start/end on that chromosome.
 Example :
 Returns : A Assembly object list
 Args    : Chromosome id (chrN notation) and start, end in absolute basepairs

=cut

sub new {
    my ($class,$dbobj,$ensembldb_obj) = @_;

    &eprof_start('new_ass');
    my $self = {};
    bless $self,$class;

    if( !defined $dbobj || !ref $dbobj ) {
        $self->throw("Don't have a db [$dbobj] for new adaptor");
    }

    $self->db($dbobj);

    $self->{'ens_db'} = $ensembldb_obj;
    $self->clear_caches();
    &eprof_end('new_ass');
    return $self;
}

sub clear_caches {
    my $self = shift;
    $self->{'known_clone_cache'} = [];
    $self->{'clone_cache'} = [];
    $self->{'contig_cache'} = [];
    $self->{'known_clone_flag'} = 0;
    $self->{'clone_flag'} = 0;
    $self->{'contig_flag'} = 0;
    $self->{'vg_sten_cache'} = [];
    $self->{'vg_exon_cache'} = [];
    $self->{'vg_sten_flag'} = 0;
    $self->{'vg_exon_flag'} = 0;
    $self->{'marker_cache'} = [];
    $self->{'marker_flag'} = 0;
}
sub ens_db {
    my $self = shift;
    return $self->{'ens_db'};
}

sub id { return 'xx'; }

sub chr_length( $chr ) {
    my $self = shift;
    my $chr  = shift;
    $chr =~s/chr//;
    my $length = $self->db->_db_handle->selectrow_array(
	"select max(start+length) from Fpc_Contig where chromosome = ?", {},
	$chr
    );
    return $length-1;
}

sub set_chr_start_end_from_contig {
    my ($self,$contig,$context) = @_;

    $context ||= 500;
    my($chr,$start,$end) = 
	$self->db->_db_handle->selectrow_array(
	"select chromosome,start, start+length from Fpc_Contig where contig_name = ?", {},
	$contig
    );
	
   $self->{'__error__'} = 'This contig is not in the current FPC map' unless($chr);
   $self->set_chr_start_end( $chr, $start - $context, $end + $context - 1 );
}

sub set_chr_start_end_from_contig_start_end {
    my ($self,$contig,$_start,$_end,$context) = @_;

    $context ||= 500;
	my $mid_point = int( ($_start + $_end) /2 );
    my($chr,$start) = 
	$self->db->_db_handle->selectrow_array(
	"select chromosome,start from Fpc_Contig where contig_name = ?", {},
	$contig
    );
   $self->{'__error__'} = 'This contig is not in the current FPC map' unless($chr);
   $self->set_chr_start_end( $chr, $start+$mid_point-$context, $start + $mid_point + $context-1  );
}

sub set_chr_start_end_from_clone {
    my ($self,$clone,$context) = @_;

    &eprof_start('from_clone');
    $context ||= 500;
    my($chr,$start,$start2,$fpc_size) = $self->db->_db_handle->selectrow_array(
	"select co.chromosome,
		cl.start_guess,co.start,
		cl.fpc_size
           from Fpc_Clone as cl,Fpc_Contig as co
          where cl.clone_name = ? and cl.contig_id = co.contig_id", {},
	$clone
    );
    unless( $chr) {
    ($chr,$start,$start2,$fpc_size) = $self->db->_db_handle->selectrow_array(
	"select co.chromosome,
		cl.start_guess,co.start,
		cl.fpc_size
           from Fpc_Clone as cl,Fpc_Contig as co
          where cl.embl_id = ? and cl.contig_id = co.contig_id", {},
	$clone
    );
   }
   $self->{'__error__'} = 'This clone is not in the current FPC map' unless($chr);
   $self->{'__error__'} = 'This clone is not in the current FPC map' unless($chr);
		
    $self->set_chr_start_end( $chr, $start+$start2, $start+$start2+$fpc_size-1 );
    &eprof_end('from_clone');
}

sub set_chr_start_end_by_clone {
    my ($self,$clone,$context) = @_;

    $context ||= 500000;
    my($chr,$start1,$start2,$size) = $self->db->_db_handle->selectrow_array(
	"select co.chromosome, cl.start_guess, co.start, cl.fpc_size
           from Fpc_Clone as cl,Fpc_Contig as co
          where cl.clone_name = ? and cl.contig_id = co.contig_id", {},
	$clone
    );
    unless($chr) {
      ($chr,$start1,$start2,$size) = $self->db->_db_handle->selectrow_array(
	"select co.chromosome, cl.start_guess, co.start, cl.fpc_size
           from Fpc_Clone as cl,Fpc_Contig as co
          where cl.embl_id = ? and cl.contig_id = co.contig_id", {},
	$clone
    );
    }
    my $start1 = int($start1+$start2+($size-1)/2);
	
   $self->{'__error__'} = 'This clone is not in the current FPC map' unless($chr);
    $self->set_chr_start_end( $chr, $start1-$context, $start1+$context-1 );
}

sub set_chr_start_end {
    my ($self,$chr,$start,$end) = @_;
	return unless $chr;
    $start ||= 1;				# Default start to 1
    unless(defined $end) {
       $end   = $self->chr_length( $chr );	# Default end to end of chromomsome
       $end   = 1000000 if $end > 1000000;
    }
    $end = $self->chr_length( $chr) if($end eq 'end' );
    $self->{'_global_start'} = $start;
    $self->{'_global_end'}   = $end;
    $self->{'length'}        = $end - $start +1;
    $chr=~s/chr//;
    $self->{'_chr_name'}     = "chr$chr";
    $self->{'_chr_no'}       = $chr;
    $self->clear_caches();
}

sub set_chr {
    my ($self,$chr) = @_;
    $self->set_chr_start_end( $chr, 1, $self->chr_length($chr) );
}

sub _chr_name {
    my ($self) = @_;
    return $self->{'_chr_name'};
}
sub fetch_all_contigs {
    my ($self) =@_;

    return @{$self->{'contig_cache'}} if($self->{'contig_flag'}==1);

    my $sth = $self->prepare(
	"select co.contig_id, co.contig_name,
                co.start as g_start,
                co.length
           from Fpc_Contig as co
          where co.chromosome = ? and
                co.start < ? and
                co.start+co.length > ?");
    $sth->execute( $self->{'_chr_no'}, $self->{'_global_end'}+1, $self->{'_global_start'}-1 );
    my @contigs = ();
		
    while (my $data = $sth->fetchrow_hashref()){
    	my $contig = Bio::EnsEMBL::Assembly->new(
	        'type'        => 'contig',
		'contig_id'   => $data->{'contig_id'},
		'contig_name' => $data->{'contig_name'},
		'golden_start'=> $data->{'g_start'},
		'length'      => $data->{'length'},
		'vc_start'    => $self->{'_global_start'},
		'chromomsome' => $self->{'_chr_name'},
		'chr_no' => $self->{'_chr_no'}
	);
	push @contigs, $contig;
    }
    $self->{'contig_flag'}=1;
    $self->{'contig_cache'} = \@contigs;
    return @contigs;
}

sub fetch_all_known_clones {
    my ($self) =@_;

    return @{$self->{'known_clone_cache'}} if($self->{'known_clone_flag'}==1);

    my $sth = $self->prepare(
	    "select co.contig_id, co.contig_name, cl.embl_id, cl.clone_name,
		co.start + cl.start_guess as g_start, cl.fpc_size,
		cl.seq_size, cl.state, cl.organisation, cl.sequenced
	   from Fpc_Contig as co, Fpc_Clone as cl
	  where cl.state in ('Phase0','Finish','PreDraft','Draft','Accessioned') and cl.contig_id = co.contig_id and co.chromosome = ? and
		co.start+cl.start_guess < ? and
		co.start+cl.start_guess+cl.fpc_size > ?
	  order by state desc");

    $sth->execute( $self->{'_chr_no'}, $self->{'_global_end'}+1, $self->{'_global_start'}-1 );
    my @clones = ();
       
    while (my $data = $sth->fetchrow_hashref()){
	    my $clone = Bio::EnsEMBL::Assembly->new(
		'type'		=> 'clone',
		'centre'	=> $data->{'organisation'},
		'state'		=> $data->{'state'},
		'contig_id'   => $data->{'contig_id'},
		'contig_name' => $data->{'contig_name'},
			'clone_id'    => $data->{'embl_id'},
			'clone_name'  => $data->{'clone_name'},
			'sequenced'   => $data->{'sequenced'},
			'golden_start'=> $data->{'g_start'},
			'length'	=> $data->{'fpc_size'},
			'seq_length'	=> $data->{'seq_size'},
			'vc_start'    => $self->{'_global_start'},
			'chromomsome' => $self->{'_chr_name'},
			'chr_no' => $self->{'_chr_no'}
        );
	 push @clones, $clone;
    }
    $self->{'known_clone_flag'}=1;
    $self->{'known_clone_cache'} = \@clones;
    return @clones;
}

sub fetch_all_clones {
    my ($self) =@_;

    return @{$self->{'clone_cache'}} if($self->{'clone_flag'}==1);

    my $sth = $self->prepare(
	"select co.contig_id, co.contig_name, cl.embl_id, cl.clone_name, 
                co.start + cl.start_guess as g_start, cl.fpc_size,
                cl.seq_size, cl.state, cl.organisation, cl.sequenced,
		f.embl_id as bac_f, r.embl_id as bac_r
           from Fpc_Contig as co, Fpc_Clone as cl
		left join bacend as r on r.clone_name=cl.clone_name and r.end='r'
		left join bacend as f on f.clone_name=cl.clone_name and f.end='f'
          where cl.contig_id = co.contig_id and co.chromosome = ? and
                co.start+cl.start_guess < ? and
                co.start+cl.start_guess+cl.fpc_size > ?
          order by state desc");

    $sth->execute( $self->{'_chr_no'}, $self->{'_global_end'}+1, $self->{'_global_start'}-1 );
    my @clones = ();
		
    while (my $data = $sth->fetchrow_hashref()){
    	my $clone = Bio::EnsEMBL::Assembly->new(
	        'type'        => 'clone',
		'bac_f'       => $data->{'bac_f'},
		'bac_r'       => $data->{'bac_r'},
		'centre'      => $data->{'organisation'},
		'state'       => $data->{'state'},
		'contig_id'   => $data->{'contig_id'},
		'contig_name' => $data->{'contig_name'},
		'clone_id'    => $data->{'embl_id'},
		'clone_name'  => $data->{'clone_name'},
		'sequenced'   => $data->{'sequenced'},
		'golden_start'=> $data->{'g_start'},
			'length'	=> $data->{'fpc_size'},
			'seq_length'	=> $data->{'seq_size'},
		'vc_start'    => $self->{'_global_start'},
		'chromomsome' => $self->{'_chr_name'},
		'chr_no' => $self->{'_chr_no'}
	);
	push @clones, $clone;
    }
    $self->{'clone_flag'}=1;
    $self->{'clone_cache'} = \@clones;
    return @clones;
}

sub get_all_VirtualGenes_startend
{
    my ($self)=shift;

    my $gene;
    return @{$self->{'vg_sten_cache'}} if($self->{'vg_sten_flag'}==1);
    my @genes;
   &eprof_start('get_all_Vgenes_startend');

    my $gene_arrayref = $self->db->_db_handle->selectall_arrayref(
	"select name, start-?, end-?
           from GENES
          where start <= ? and end >= ? and chr = ? and known='yes'", {},
	$self->{'_global_start'}, $self->{'_global_start'},
	$self->{'_global_end'}, $self->{'_global_start'},
	$self->{'_chr_no'} );
	
    my $entryAdaptor = $self->ens_db->get_DBEntryAdaptor();
    foreach my $row ( @$gene_arrayref ) {
    	my $start = $row->[1];
    	my $end   = $row->[2];
    	my $name  = $row->[0];
        my $gene=Bio::EnsEMBL::Gene->new();
        my @gene_xrefs = $entryAdaptor->fetch_by_gene($name);
    
        foreach my $genelink ( @gene_xrefs ) {
            $gene->add_DBLink($genelink);
        }
    	my $query1 = "select t.translation from transcript t where t.gene = ?";
    	my $sth1 = $self->{'ens_db'}->prepare($query1);
    	$sth1->execute( $name );
    	
    	while (my $transid = $sth1->fetchrow) {
    	    my @transcript_xrefs = $entryAdaptor->fetch_by_translation($transid);
	        foreach my $translink(@transcript_xrefs) {
	        	$gene->add_DBLink($translink);
    	    }
    	}


        my $vg = Bio::EnsEMBL::VirtualGene->new(
            -id     => $name,
            -gene   => $gene,
            -contigid => 'xx',
            -start  => $start,
            -end    => $end,
            -strand => 1
        );

        push @genes,$vg;
    }

    $self->{'vg_sten_flag'}=1;
    $self->{'vg_sten_cache'} = \@genes;
   &eprof_end('get_all_Vgenes_startend');
    return @genes;

}

sub get_all_Genes_exononly{
   my ($self) = @_;

   return @{$self->{'vg_exon_cache'}} if($self->{'vg_exon_flag'}==1);

   &eprof_start('get_all_Genes_exononly');
   my @clones      = $self->_get_sequenced_clones();

   my $query = "SELECT e.id, e.sticky_rank, et.rank, et.transcript, t.gene,
                       e.seq_start + co.offset + ? as start,
                       e.seq_end   + co.offset + ? as end,
		       e.strand
                  FROM exon e, exon_transcript et, transcript t, contig co, clone cl
                 WHERE cl.id = ? and cl.internal_id = co.clone and
                       co.internal_id = e.contig and
                       et.exon = e.id and et.transcript = t.id
                 ORDER BY t.gene,t.id,et.rank,e.sticky_rank";

   my $sth = $self->ens_db->_db_handle->prepare($query);
   my @out;
   my @trans;
   my $length = $self->length;

 foreach my $clone (@clones) {
   my $cl_start = $clone->start - int( ($clone->seq_length - $clone->length)/2 );
   $sth->execute( $cl_start, $cl_start, $clone->id);

   my ($exonid,$stickyrank,$rank,$transcriptid,$geneid,$start,$end,$strand);
   $sth->bind_columns(undef,\$exonid,\$stickyrank,\$rank,\$transcriptid,\$geneid,\$start,\$end,\$strand);

   my $current_transcript;
   my $current_gene;
   my $current_transcript_id;
   my $current_gene_id;
   my $previous_exon;


   while( $sth->fetch ) {

       if (($end > $length) || ($start < 1)) {
           next;
       }
       if( $geneid ne $current_gene_id ) {
           # make a new gene
           $current_gene = Bio::EnsEMBL::Gene->new;
           $current_gene->id($geneid);
           push(@out,$current_gene);
           $current_gene_id = $geneid;
       }

       if( $transcriptid ne $current_transcript_id ) {
           # make a new transcript
           $current_transcript = Bio::EnsEMBL::WebTranscript->new();
           $current_gene->add_Transcript($current_transcript);
           push(@trans,$current_transcript);
           if( $rank == 1 ) {
               $current_transcript->is_start_exon_in_context('dummy',1);
           } else {
               $current_transcript->is_start_exon_in_context('dummy',0);
           }

           $current_transcript_id = $transcriptid;
           $current_transcript->id($transcriptid);
       }

       if( $stickyrank > 1 ) {
           if( !defined $previous_exon ) {
               $self->warn("Really bad news - half-on-half off Sticky Exon. Faking it");
           }
           if( $previous_exon->end < $end ) {
               $previous_exon->end($end);
               next;
           }

       }


       my $exon = Bio::EnsEMBL::Exon->new();
       $exon->start($start);
       $exon->end($end);
       $exon->strand($strand);
       $exon->id($exonid);
       $exon->seqname('xx'); # $self->id);
       $previous_exon = $exon;
       $current_transcript->add_Exon($exon);
       $current_transcript->end_exon_rank($rank);

   }
 }
   #
   # We need to make another quick trip to the database for each
   # transcript to discover whether we have all of it or not
   #

   foreach my $trans ( @trans ) {
       my $sth2 = $self->ens_db->prepare("select max(rank) from exon_transcript where transcript = '".$trans->id."'");
       $sth2->execute;
       my ($rank) = $sth2->fetchrow_array();
       if( $rank == $trans->end_exon_rank) {
           $trans->is_end_exon_in_context('xx',1);
       } else {
           $trans->is_end_exon_in_context('xx',0);
       }

   }


   #
   # This can obviously be optimised to a better single trip
   # to the database
   #

   my $gene_obj = $self->ens_db->gene_Obj;

   foreach my $g ( @out ) {
       $gene_obj->_get_dblinks($g);
       $gene_obj->_get_description($g);
   }

   $self->{'vg_exon_flag'}=1;
   $self->{'vg_exon_cache'} = \@out;
  
   &eprof_end('get_all_Genes_exononly');
   return @out;

}

sub get_all_ExternalGenes() {
	return ();	
}

sub _get_sequenced_clones {
    my $self = shift;
    my @clones = $self->fetch_all_clones();
    
    &eprof_start('_gsc');
    my @Q = map { $_->{'sequenced'} eq 'yes' ? $_ : () } @clones;
    &eprof_end('_gsc');
    return @Q;
}

sub DESTROY { return 1; }

sub AUTOLOAD {
    my $self = shift;
#    no strict 'refs';
    my $var = $AUTOLOAD;
    $var =~ s/.*:://;		# remove class name if included...
    return $self->{$var} if (defined $self->{$var});

    print STDERR "Assembly adaptor warning - \$self->{'$var'} is undefined\n";
    $self->throw();
    return undef;
}

sub fetch_karyotype_adaptor {
   my ($self) = @_;
   return ($self->ens_db->get_KaryotypeBandAdaptor());
}

sub fetch_karyotype_band_by_name {
   my ($self,$chr, $band) = @_;

   my $kadp = $self->ens_db->get_KaryotypeBandAdaptor();
   my $kband = $kadp->fetch_by_chromosome_name($chr, $band);

   return $kband;
}

sub fetch_chromosome_length {
   my $self = shift;
   return $self->chr_length( $self->{'_chr_name'} );
}

sub get_landmark_MarkerFeatures{
    my ($self,$glob) = @_;

    return @{$self->{'marker_cache'}} if($self->{'marker_flag'}==1);
    my $chr_name   = $self->{'_chr_name'};
    $chr_name =~ s/chr//g;

    my $length     = $self->length;
   
    my $sth = $self->db->prepare(
            " SELECT  start, end, strand, name 
		        FROM  contig_landmarkMarker 
		       WHERE  chr_name = ? AND
                      start >= ? AND end <= ?
		       ORDER BY start"
   );
   $sth->execute( $chr_name, $self->{'_global_start'}, $self->{'_global_end'} );
   
   my ($start, $end, $strand, $name);
   
   $sth->bind_columns( undef, \$start, \$end,  \$strand, \$name);
   
   my @markers;
   my $prev;
   while( $sth->fetch ) {
       #############################
       # change to local coordinates
       #############################
       $start = $start - $self->{'_global_start'};
       $end   = $end   - $self->{'_global_start'};
       
       if( defined $prev && $prev->end + $glob > $start  && $prev->id eq $name ) {           
           next;
       }

       my $sf 	  = Bio::EnsEMBL::SeqFeature->new();
       $sf->start(	$start  );
       $sf->end(	$end    );
       $sf->strand(	$strand );
       $sf->id(		$name   );
       push @markers, $sf;
       $prev      = $sf;
   } 

    $self->{'marker_flag'}  = 1;
    $self->{'marker_cache'} = \@markers;
    return @markers;
}

sub get_summary {
	my($self) = @_;
	my @X = $self->db->_db_handle->selectrow_array(
		"select sum(known = 'yes'), count(name) from GENES
			where chr = ?", {}, $self->{'_chr_no'}
	);
	return ($X[0],$X[1]-$X[0]);
}
1;
