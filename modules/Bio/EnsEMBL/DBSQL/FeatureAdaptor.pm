#
# EnsEMBL module for Bio::EnsEMBL::DBSQL::FeatureAdaptor
#
# Cared for by Imre Vastrik <vastrik@ebi.ac.uk>
#
# Copyright Imre Vastrik
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::FeatureAdaptor - MySQL database adapter class for EnsEMBL Feature Objects

=head1 SYNOPSIS



=head1 DESCRIPTION



=head1 CONTACT



=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are
usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::DBSQL::FeatureAdaptor;

use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object


use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::FeatureFactory;
use Bio::EnsEMBL::TranscriptFactory;



@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

# new() is inherited from Bio::EnsEMBL::DBSQL::BaseAdaptor


=head2 delete_by_RawContig

 Title   : delete_by_RawContig_id
 Usage   : $fa->delete_by_RawContig($contig)
 Function: deletes features and repeatfeatures by Bio::EnsEMBL::DB::RawContigI.
           Gets out the internal_id from RawContig_obj and passes it to
           delete_by_RawContig_internal_id, which does the business
 Example : $fa->delete_by_RawContig($RawContig_obj)
 Returns : nothing
 Args    : Bio::EnsEMBL::DB::RawContigI

=cut


sub delete_by_RawContig {
    my ($self,$contig) = @_;
    my $contig_internal_id;
    $self->throw("$contig is not a Bio::EnsEMBL::DB::RawContigI")
        unless (defined($contig) && $contig->isa("Bio::EnsEMBL::DB::RawContigI"));
    $contig_internal_id = $contig->internal_id;
    $self->delete_by_RawContig_internal_id($contig_internal_id);
}


=head2 delete_by_RawContig_internal_id

 Title   : delete_by_RawContig_id
 Usage   : $fa->delete_by_RawContig($contig)
 Function: deletes features and repeatfeatures by Bio::EnsEMBL::DB::RawContigI
 Example : $fa->delete_by_RawContig($internal_id)
 Returns : nothing
 Args    : Bio::EnsEMBL::DB::RawContigI

=cut


sub delete_by_RawContig_internal_id {
    my ($self,$contig_internal_id) = @_;
    $contig_internal_id || $self->throw("I need contig internal id");
    $contig_internal_id =~ /^\d+$/ or $self->warn("[$contig_internal_id] does not look like internal id.");

    #print(STDERR "Deleting features for contig $contig\n");
    my $sth = $self->db->prepare("delete from feature where contig = '$contig_internal_id'");
    my $res = $sth->execute;

    #print(STDERR "Deleting repeat features for contig $contig\n");
    $sth = $self->db->prepare("delete from repeat_feature where contig = '$contig_internal_id'");
    $res = $sth->execute;
}


=head2 delete_by_RawContig_id

 Title   : delete_by_RawContig_id
 Usage   : $fa->delete_by_RawContig_id($contig_id)
 Function: deletes features and repeatfeatures by contig_id.
           Gets contig_internal_id via RawContigAdaptor and passes that
           to delete_by_RawContig_internal_id, which does the business
 Example : $fa->delete_by_RawContig_id('AC004042.00001')
 Returns : nothing
 Args    : contig_id e.g. AC004042.00001

=cut


sub delete_by_RawContig_id {
    my ($self,$contig_id) = @_;
    my $contig_internal_id = $self->db->get_RawContigAdaptor->get_internal_id_by_id($contig_id);
    unless($contig_internal_id) {
	$self->warn("Could not get internal_id for contig '$contig_id'");
	return;
    }
    $self->delete_by_RawContig_internal_id($contig_internal_id);
}


=head2 write

 Title   : write
 Usage   :
 Function: deprecates
 Example :
 Returns :
 Args    :

=cut


sub write {
    my ($self,$contig,@features) = @_;
    $self->warn("Deprecated method 'FeatureAdaptor->write', passing to 'FeatureAdaptor->store' instead");
    $self->store($contig,@features);
}


=head2 store

 Title   : store
 Usage   : $fa->store($contig,@features)
 Function: Writes a feature on the genomic sequence of a contig into the database.
           Checks that we have contig_obj, gets its internal_id. Checks that each
           feature_obj has analysis_obj attached, gets analysis_id from it.
           Checks what kind of feature(s) is/are passed in and passes it/them
           further to appropriate _store_blabla function together with contig_internal_id
           and analysis_id.
 Example :
 Returns : nothing
 Args    : Bio::EnsEMBL::DB::ContigI, Bio::EnsEMBL::SeqFeatureI

=cut


sub store {
    my ($self,$contig,@features) = @_;


    my ($p,$f,$l) = caller;
    $self->warn("$f:$l FeatureAdaptor store being phased out. It is better to use the new FeatureAdaptors directly (more type safe)");

    my $repeat_adaptor = $self->db->get_RepeatFeatureAdaptor();
    my $dna_align_adaptor = $self->db->get_DnaAlignFeatureAdaptor();
    my $protein_align_adaptor = $self->db->get_ProteinAlignFeatureAdaptor();
    my $simple_adaptor = $self->db->get_SimpleFeatureAdaptor();
#    my $prediction_adaptor = $self->db->get_PredictionFeatureAdaptor();



    # Check for contig
    $self->throw("$contig is not a Bio::EnsEMBL::DB::ContigI")
        unless (defined($contig) && $contig->isa("Bio::EnsEMBL::DB::ContigI"));
    my $contig_internal_id = $contig->dbID();


    #
    # Don't particularly like this loop, but I guess we should stick
    # with it. EB
    #

    FEATURE :
    foreach my $feature ( @features ) {	

	# Check that the thingy passed in is of "right kind"
	if( ! $feature->isa('Bio::EnsEMBL::SeqFeatureI') ) {
	    $self->throw("Feature $feature is not a feature!");
	}
# 	eval {
# 	    $feature->validate();
# 	};
# 	if ($@) {
# 	    next FEATURE;
# 	}

# 	# Check that we have Analysis
# 	my $analysisid;
# 	if (!defined($feature->analysis)) {
# 	    $self->throw("Feature " . $feature->seqname . " " .
# 			              $feature->source_tag .
# 			 " doesn't have analysis. Can't write to database");
# 	} else {
# 	    # Get AnalysisAdaptor if we haven't got one
# 	    unless($feature->analysis->adaptor) {
# 		$feature->analysis->adaptor($self->db->get_AnalysisAdaptor);
# 	    }
# 	    $analysisid = $feature->analysis->adaptor->store($feature->analysis);
# 	}

	#
	# Retarget to new adaptor scheme
	#

	if( $feature->isa('Bio::EnsEMBL::DnaPepAlignFeature') ) {
	    $protein_align_adaptor->store($contig_internal_id,$feature);
	} elsif ( $feature->isa('Bio::EnsEMBL::DnaDnaAlignFeature') ) {
	    $dna_align_adaptor->store($contig_internal_id,$feature);
	} elsif ( $feature->isa('Bio::EnsEMBL::RepeatFeature') ) {
	    $repeat_adaptor->store($contig_internal_id,$feature);
	} elsif ( $feature->isa('Bio::EnsEMBL::SimpleFeature') ) {
	    $simple_adaptor->store($contig_internal_id,$feature);
#	} elsif ( $feature->isa('Bio::EnsEMBL::PredictionFeature') ) {
#	    $prediction_adaptor->store($contig_internal_id,$feature);
	} else {
	    $self->throw("cannot store $feature - no feature adaptor that fits it!");
	}


    }
}


=head2 _store_FeaturePair

 Title   : _store_FeaturePair
 Usage   : $fa->_store_FeaturePairs($contig_internal_id,$analysisid,$feature)
 Function: internal method for storing a Bio::EnsEMBL::FeaturePair
 Example :
 Returns : nothing
 Args    : contig_internal_id, analysis_d, Bio::EnsEMBL::FeaturePair

=cut


sub _store_FeaturePair {
    my ($self,$contig_internal_id,$analysisid,$feature) = @_;
    my $homol = $feature->feature2;
    $self->_store
    (
	 'NULL',
	 $contig_internal_id,
	 $feature->start,
	 $feature->end,
	 $feature->score,
	 $feature->strand,
	 $analysisid,
	 $feature->source_tag,
	 $homol->start,
	 $homol->end,
	 $homol->seqname,
	 ((defined $feature->p_value)      ? &exponent($feature->p_value)     : 'NULL'),
	 ((defined $feature->percent_id)   ? $feature->percent_id  : 'NULL'),
	 ((defined $feature->phase)        ? $feature->phase       : 'NULL'),
	 ((defined $feature->end_phase)    ? $feature->end_phase   : 'NULL')
    );

}


=head2 _store_single_feature

 Title   : _store_single_feature
 Usage   : $fa->_store_single_feature($contig_internal_id,$analysisid,$feature)
 Function: internal method for storing a Bio::EnsEMBL::SeqFeature
 Example :
 Returns : nothing
 Args    : contig_internal_id, analysis_d, Bio::EnsEMBL::SeqFeature

=cut


sub _store_single_feature
{
    my ($self,$contig_internal_id,$analysisid,$feature) = @_;
    $self->_store
    (
	 'NULL',
	 $contig_internal_id,
	 $feature->start,
	 $feature->end,
	 $feature->score,
	 $feature->strand,
	 $analysisid,
	 $feature->source_tag,
	 -1,
	 -1,
	 "__NONE__",
	 'NULL',
	 'NULL',
	 'NULL',
	 'NULL'
   )
}


=head2 _store

 Title   : _store_single_feature
 Usage   : $fa->_store($id,$contig,$seq_start,$seq_end,$score,$strand,$analysis,
           $name,$hstart,$hend,$hid,$evalue,$perc_id,$phase,$end_phase)
 Function: internal method for which store a feature into the feature table
 Example :
 Returns : last insert id
 Args    : $id,$contig,$seq_start,$seq_end,$score,$strand,$analysis,
           $name,$hstart,$hend,$hid,$evalue,$perc_id,$phase,$end_phase

=cut


sub _store {
    my $self = shift;
    my $query = "insert into feature (id,contig,seq_start,seq_end,score,strand,analysis,name,hstart,hend,hid,evalue,perc_id,phase,end_phase) values(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)";
    my $sth = $self->db->prepare($query);
#    print STDERR join(", ", @_), "\n";
    $sth->execute(@_);
    return $sth->{mysql_insertid};
}


=head2 _store_PredictionFeature

 Title   : _store_PredictionFeature
 Usage   : $fa->_store_PredictionFeature($contig_internal_id,$analysisid,@features)
 Function: internal method for storing Bio::EnsEMBL::SeqFeature(s) with sub-SeqFeatures
 Example :
 Returns : nothing
 Args    : contig_internal_id, analysis_d, Bio::EnsEMBL::SeqFeature

=cut


sub _store_PredictionFeature {
    my ($self,$contig_internal_id,$analysisid,@features) = @_;
    foreach my $feature ( @features ) {

	# now write each sub feature

	foreach my $sub ( $feature->sub_SeqFeature ) {
	    my $last_insert_id = $self->_store
	    (
		 'NULL',
		 $contig_internal_id,
		 $sub->start,
		 $sub->end,
		 (defined($sub->score) ? $sub->score : -1000),
		 $sub->strand,
		 $analysisid,
		 $sub->source_tag,
		 -1,
		 -1,
		 ($sub->primary_tag || "__NONE__"),
		 ((defined $sub->p_value)     ?   &exponent($sub->p_value) : 'NULL'),
                 ((defined $sub->percent_id)  ?   $sub->percent_id    : 'NULL'),
                 ((defined $sub->phase)       ?   $sub->phase         : 'NULL'),
                 ((defined $sub->end_phase)   ?   $sub->end_phase     : 'NULL')
	    );
	}
    }
}


=head2 _store_Repeat

 Title   : _store_Repeat
 Usage   : $fa->_store_Repeat($contig_internal_id,$analysisid,@features)
 Function: internal method for storing Bio::EnsEMBL::RepeatI(s).
           Writes into repeat_feature table
 Example :
 Returns : nothing
 Args    : contig_internal_id, analysis_d, Bio::EnsEMBL::RepeatI

=cut


sub _store_Repeat {
    my ($self,$contig_internal_id,$analysisid,@features) = @_;
    # Since repeats have their own table we are not using &_store here

    my $sth = $self->db->prepare("insert into repeat_feature(id,contig,seq_start,seq_end,score,strand,analysis,hstart,hend,hid) values(?,?,?,?,?,?,?,?,?,?)");

    foreach my $feature ( @features ){	
	if( ! $feature->isa('Bio::EnsEMBL::RepeatI') ) {
	    $self->throw("Feature $feature is not a Repeat!");
	}
	eval {
	    $feature->validate();
	};
	if ($@)	{
	    next;
	}
	my $homol = $feature->feature2;
	$sth->execute(
		      'NULL',
		      $contig_internal_id,
		      $feature->start,
		      $feature->end,
		      $feature->score,
		      $feature->strand,
		      $analysisid,
		      $homol->start,
		      $homol->end,
		      $homol->seqname
		      );
    }
}


=head2 find_GenomeHits

 Title   : find_GenomeHits
 Usage   :
 Function: deprecated
 Example :
 Returns :
 Args    :

=cut


sub find_GenomeHits {
    my ($self,$arg) = @_;
    $self->warn("Deprecated method 'FeatureAdaptor->find_GenomeHits', passing to 'FeatureAdaptor->fetch_by_hid' instead");
    $self->fetch_by_hid($arg);
}


=head2 fetch_by_hid

 Title   : fetch_by_hid
 Usage   : my @features = $fa->fetch_by_hid($hid)
 Function: fetches feature_objs from feature table by hid coloumn
 Example : my @features = $fa->fetch_by_hid('IL5_MACMU')
 Returns : array of Bio::Ensembl::FeaturePair
 Args    : hid (string)

=cut


sub fetch_by_hid {
    my ($self,$arg) = @_;

    $self->throw("No hit id input") unless defined($arg);
    my $query = "select c.id, " .
	                "f.seq_start, " .
			"f.seq_end, "   .
			"f.score, "     .
			"f.strand, "    .
			"f.analysis, "  .
			"f.name, "      .
			"f.hstart, "    .
			"f.hend, "      .
			"f.hid, "       .
			"f.evalue, "    .
			"f.perc_id, "   .
			"f.phase, "     .
			"f.end_phase "   .
	        "from   feature as f,contig as c " .
		"where  f.hid = '$arg' and " .
		        "c.internal_id = f.contig";

    my $sth   = $self->db->prepare($query);
    my $res   = $sth->execute;
    my ($contig,$start,$end,$score,$strand,$analysisid,$name,$hstart,$hend,$hid,
	$evalue,$perc_id,$phase,$end_phase);
    $sth->bind_columns(undef,\$contig,\$start,\$end,\$score,\$strand,\$analysisid,
		       \$name,\$hstart,\$hend,\$hid,\$evalue,\$perc_id,\$phase,
		       \$end_phase);

    my @features;
    my $analysisadaptor = $self->db->get_AnalysisAdaptor;
    while($sth->fetch) {
	my $out;
	$name = defined($name) ? $name : 'no_source';
	$out = Bio::EnsEMBL::FeatureFactory->new_feature_pair();
	$out->set_all_fields($start,$end,$strand,$score,$name,'similarity',$contig,
			     $hstart,$hend,1,$score,$name,'similarity',$hid,$evalue,
			     $perc_id,$phase,$end_phase);
	$out->analysis($analysisadaptor->fetch_by_dbID($analysisid));
	$out->validate;
	push(@features,$out);
	
    }
    return @features;
}


=head2 get_PredictionFeature_by_id

 Title   : get_PredictionFeature_by_id
 Usage   :
 Function: deprecated
 Example :
 Returns :
 Args    :

=cut


sub get_PredictionFeature_by_id {
    my ($self,$genscan_id)=@_;
    $self->warn("Deprecated method 'FeatureAdaptor->get_PredictionFeature_by_id', passing to 'FeatureAdaptor->fetch_PredictionFeature_by_id' instead");
    $self->fetch_PredictionFeature_by_id($genscan_id);
}


=head2 fetch_PredictionFeature_by_id

 Title   : get_PredictionFeature_by_id
 Usage   : $fa->get_PredictionFeature_by_id($genscan_id)
 Function: Fetches a genscan prediction by fset_feature.fset
 Example : my $feature = $fa->fetch_PredictionFeature_by_id(194643)
 Returns : Bio::EnsEMBL::SeqFeature with sub-SeqFeatures
           Throws an exception if fset with genscan_id does not exist in the db
 Args    : genscan_id

=cut


sub fetch_PredictionFeature_by_id {

    my ($self,$genscan_id) = @_;

    $genscan_id || $self->throw("I need a genscan id");

#    my $query = "select f.id,f.seq_start,f.seq_end,f.strand,f.score,f.analysis,f.name,f.hid,fset.id,c.id,f.phase " .
#       "from feature f, fset fset,fset_feature ff,contig c where ff.feature = f.id and fset.id = ff.fset ".
#        " and c.internal_id=f.contig and ff.fset ='$genscan_id' and name = 'genscan'";
#    my $sth = $self->db->prepare($query);
#    $sth->execute();
#    my ($fid,$start,$end,$strand,$score,$analysisid,$name,$hid,$fsetid,$contig,$phase);
#    $sth->bind_columns(undef,\$fid,\$start,\$end,\$strand,\$score,\$analysisid,\$name,\$hid,\$fsetid,\$contig,\$phase);
#    my $analysisadaptor = $self->db->get_AnalysisAdaptor;
#    my $current_fset;
#    while( $sth->fetch ) {
#	unless($current_fset) {
#	    $current_fset = Bio::EnsEMBL::FeatureFactory->new_feature;
#	    $current_fset = new Bio::EnsEMBL::SeqFeature;
#	    $current_fset->source_tag($name);
#	    $current_fset->primary_tag('prediction');
#	    $current_fset->analysis($analysisadaptor->fetch_by_dbID($analysisid));
#	    $current_fset->seqname($contig);
#	    $current_fset->raw_seqname($contig);
#	    $current_fset->id($fsetid);
#	    $current_fset->score(defined($score) ? $score : undef);
#	    $current_fset->strand($strand);
#        }

#	my $out = Bio::EnsEMBL::FeatureFactory->new_feature;
#	$out->seqname($contig);
#	$out->start($start);
#	$out->end($end);
#	$out->strand($strand);
#	$out->phase($phase);
#	$out->source_tag($name);
#	$out->primary_tag($hid);
#	$out->score(defined($score) ? $score : undef);
#	$out->analysis($analysisadaptor->fetch_by_dbID($analysisid));

	# Final check that everything is ok.
#	$out->validate();

#	$current_fset->add_sub_SeqFeature($out,'EXPAND');
#    }
#    $current_fset || $self->throw("Fset $genscan_id does not exist in the database");

    my $contigid = $genscan_id;

    $contigid =~ s/(.*?)\..*/$1/;

    my $contig =  $self->db->get_Contig_by_internal_id($contigid);

    my @fsets = $contig->get_all_PredictionFeatures;

    foreach my $fset (@fsets) {
      if (defined($fset->sub_SeqFeature)) {
	my @f = $fset->sub_SeqFeature;

	if ($f[0]->id eq $genscan_id) {
	  return $fset;
	}
      }
    }
    $self->throw("Fset $genscan_id does not exist in the dataabase");

}




=head2 get_PredictionFeature_as_Transcript

 Title   : get_PredictionFeature_as_Transcript
 Usage   :
 Function: depredated
 Example :
 Returns :
 Args    :


=cut

sub get_PredictionFeature_as_Transcript {
    my ($self,$genscan_id)=@_;
    $self->warn("Deprecated method 'FeatureAdaptor->get_PredictionFeature_as_Transcript', passing to 'FeatureAdaptor->fetch_PredictionFeature_as_Transcript' instead");
    $self->fetch_PredictionFeature_as_Transcript($genscan_id);
}


=head2 fetch_PredictionFeature_as_Transcript

 Title   : fetch_PredictionFeature_as_Transcript
 Usage   : $fa->fetch_PredictionFeature_as_Transcript($genscan_id)
 Function: Passes the genscan_id to fetch_PredictionFeature_by_id and returns the
           resulting SeqFeature as Bio::EnsEMBL::Transcript.
 Example : my $transcript = $fa->fetch_PredictionFeature_as_Transcript(194643)
 Returns : Bio::EnsEMBL::Transcript
 Args    : genscan_id

=cut


sub fetch_PredictionFeature_as_Transcript{
    my ($self,$genscan_id)=@_;
    $genscan_id || $self->throw("I need a genscan id");

    my $ft=$self->fetch_PredictionFeature_by_id($genscan_id);
    my $contig=$self->db->get_Contig($ft->seqname);

    return &Bio::EnsEMBL::TranscriptFactory::fset2transcript($ft,$contig);
}


sub exponent {
    my ($number) = @_;

    my ($exp) = sprintf("%.3e", $number);
    return $exp;
}


=head2 delete

 Title   : delete
 Usage   :
 Function: deletes all features from a contig;
 Example :
 Returns :
 Args    :


=cut

sub delete {
    my ($self,$contig) = @_;
    $self->warn("Deprecated method FeatureAdaptor->delete, passing to 'FeatureAdaptor->delete_by_RawContig_internal_id' instead");
    $self->delete_by_RawContig_internal_id($contig);
}


=head2 fetch_RepeatFeatures_by_RawContig

 Title   : fetch_RepeatFeatures_by_RawContig
 Usage   : foreach my $rf ($fa->fetch_all_RepeatFeatures($contig))
 Function: Gets all the repeat features on a contig.
           If the thingy passed in is Bio::EnsEMBL::DB::ContigI, gets the internal_id from
           there. Otherwise assumes that the thingy is contig id, which is used to get
           contig internal_id via RawContigAdaptor.
 Example :
 Returns : Array of Bio::EnsEMBL::Repeat
 Args    : Bio::EnsEMBL::DB::ContigI or contig id as a string


=cut


sub fetch_RepeatFeatures_by_RawContig {
   my ($self, $contig) = @_;
   my $contig_internal_id;
   if (ref($contig) && $contig->isa("Bio::EnsEMBL::DB::ContigI")) {
       # we have ContigI object
       $contig_internal_id = $contig->internal_id;
   } elsif (defined($contig)) {
       # assume that the thing passed in is Contig id
       $contig_internal_id = $self->db->get_RawContigAdaptor->get_internal_id_by_id($contig);
   } else {
       $self->throw("I need contig id");
   }

   my @array;

   # make the SQL query
   my $statement = "select id,seq_start,seq_end,strand,score,analysis,hstart,hend,hid " .
		   "from repeat_feature where contig = '$contig_internal_id'";

   my $sth = $self->db->prepare($statement);

   $sth->execute();

   my ($fid,$start,$end,$strand,$score,$analysisid,$hstart,$hend,$hid);

   # bind the columns
   $sth->bind_columns(undef,\$fid,\$start,\$end,\$strand,\$score,\$analysisid,\$hstart,\$hend,\$hid);

   my $analysisadaptor = $self->db->get_AnalysisAdaptor;

   while( $sth->fetch ) {
       my $out;
       if( $hid ne '__NONE__' ) {
	   # is a paired feature
	   # build EnsEMBL features and make the FeaturePair

	   $out = Bio::EnsEMBL::FeatureFactory->new_repeat();
	   $out->set_all_fields($start,$end,$strand,$score,'repeatmasker','repeat',$fid,
				$hstart,$hend,1,$score,'repeatmasker','repeat',$hid);

	   $out->analysis($analysisadaptor->fetch_by_dbID($analysisid));
	   $out->id($fid);
       } else {
	   $self->warn("Repeat feature does not have a hid. bad news....");
       }
       $out->validate();
       push(@array,$out);
  }
  return @array;
}                                       # get_all_RepeatFeatures


1;
