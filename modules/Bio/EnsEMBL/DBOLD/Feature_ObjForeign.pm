#
# EnsEMBL module for Bio::EnsEMBL::DBOLD::Feature_ObjForeign
#
# Cared for by Elia Stupka <elia@ebi.ac.uk>, Philip lijnzaad@ebi.ac.uk
#
# Copyright EnsEMBL
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBOLD::Feature_ObjForeign, specialization of Feature_Obj that
allows you to translate table names to other names. The intended use is
for an database that reuses (readonly, not writing!) part of another
database in order to save on speed (no loading of all the assembly related
data) and disk usage.


=head1 SYNOPSIS

  use Bio::EnsEMBL::DBOLD::Obj;
  use Bio::EnsEMBL::DBOLD::Feature_ObjForeign;

  $db = new Bio::EnsEMBL::DBOLD::Obj( -user => 'root', -db => 'pog' , -host => 'caldy' , -driver => 'mysql' );
  my $feature_obj =Bio::EnsEMBL::Feature_ObjForeign->new($obj);

  my $mapping = { feature => myfeature, 
                  static_golden_path => ens081.static_golden_path
                 };
  my $feature_obj->table_name_tranlations $mapping;

  # From here on, all the application code first translates the table
  # names, then does the queries on possibly translated names.

=head1 DESCRIPTION

Bio::EnsEMBL::DBOLD::Feature_ObjForeign, specialization of Feature_Obj that
allows you to translate table names to other names. The intended use is
for an database that reuses (readonly, not writing!) part of another
database in order to save on speed (no loading of all the assembly related
data) and disk usage.

=head1 CONTACT

Elia Stupka: elia@ebi.ac.uk, Philip lijnzaad@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal
methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::DBOLD::Feature_ObjForeign;

use vars qw(@ISA $AUTOLOAD);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::EnsEMBL::Root;
#use Bio::EnsEMBL::DBOLD::Obj;

use Bio::EnsEMBL::Ghost;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::FeatureFactory;
use DBI;
use Bio::EnsEMBL::DBOLD::Utils;
use Bio::EnsEMBL::DBOLD::DummyStatement;


@ISA = qw(Bio::EnsEMBL::DBOLD::ObjForeign
          Bio::EnsEMBL::Root 
          Bio::EnsEMBL::DBOLD::Feature_Obj
         );

sub new {  # comes from ObjForeign
    my ($class, @args) = @_;
    my $self= $class->SUPER::new(@_);
    bless $self, $class;
    $self;
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

    if (ref( $contig) && $contig->isa("Bio::EnsEMBL::DB::ContigI")) {
	$self->throw("You have to give a contig id, not a contig object!");
    }
    

    my $sth = $self->_db_obj->prepare("SELECT fs.feature,fs.fset " .
			     "from   ".$self->_lookup_table_name('fset_feature')." as fs, " .
			     "       ".$self->_lookup_table_name('feature')." as f " .
			     "where  fs.feature = f.id " .
			     "and    f.contig = '$contig'");

    my $res = $sth->execute || $self->warn("Could not find features for contig $contig");

    my %fset;

    while (my $rowhash = $sth->fetchrow_hashref) {
	$fset{$rowhash->{fset}} = 1;
    }
    
    my @fset = keys %fset;
    
    if ($#fset >= 0) {
	my $fsstr = "";

	foreach my $fs (@fset) {
	    $fsstr .= $fs . ",";
	}

	chop($fsstr);
	

	
	$sth = $self->_db_obj->prepare("DELETE from ".$self->_lookup_table_name('fset')." where id in ($fsstr)");
	$res = $sth->execute;
	
	$sth = $self->_db_obj->prepare("DELETE from ".$self->_lookup_table_name('fset_feature')." where fset in ($fsstr)");
	$res = $sth->execute;
    }
    
    #print(STDERR "Deleting features for contig $contig\n");
    
    $sth = $self->_db_obj->prepare("DELETE from ".$self->_lookup_table_name('feature')." where contig = '$contig'");
    $res = $sth->execute;
    
    #print(STDERR "Deleting repeat features for contig $contig\n");
    
    $sth = $self->_db_obj->prepare("DELETE from ".$self->_lookup_table_name('repeat_feature')." where contig = '$contig'");
    $res = $sth->execute;

}

=head2 write

 Title   : write
 Usage   : $obj->write($contig,@features)
 Function: Writes a feature on the genomic sequence of a contig into the database
 Example :
 Returns : nothing
 Args    : Bio::EnsEMBL::SeqFeatureI


=cut

sub write {
    my ($self,$contig,@features) = @_;

    #
    # Yes - we need to rethink how we are writing features into the
    # database. This is a little obtuse and clunky now
    #

    $self->throw("$contig is not a Bio::EnsEMBL::DB::ContigI") 
        unless (defined($contig) && $contig->isa("Bio::EnsEMBL::DB::ContigI"));
    
    my $contigid = $contig->id;
    my $analysis;

    
    # Put the repeats in a different table, and also things we need to write
    # as fsets.
    my @repeats;
    my @fset;

    FEATURE :
    foreach my $feature ( @features ) {	

	if( ! $feature->isa('Bio::EnsEMBL::SeqFeatureI') ) {
	    $self->throw("Feature $feature is not a feature!");
	}

	eval {
	    $feature->validate();
	};

	if ($@) {

	    next FEATURE;
	}
	
	if($feature->isa('Bio::EnsEMBL::RepeatI')) {
	    push(@repeats,$feature);
	} elsif ( $feature->sub_SeqFeature ) {
	    push(@fset,$feature);
	} else {    
	    if (!defined($feature->analysis)) {
		$self->throw("Feature " . $feature->seqname . " " . $feature->source_tag ." doesn't have analysis. Can't write to database");
	    } else {
		$analysis = $feature->analysis;
	    }


	    my $analysisid = $self->write_Analysis($analysis);

	    if ( $feature->isa('Bio::EnsEMBL::FeaturePair') ) {
		my $homol = $feature->feature2;
	    my $sth = $self->_db_obj->prepare(  
                  "INSERT into ".$self->_lookup_table_name('feature')."(id,contig,seq_start,seq_end,score,strand,name,analysis,hstart,hend,hid,perc_id,evalue,phase,end_phase) ".
                  "values ('NULL',"
                  .$contig->internal_id     .","
                  .$feature->start          .","
                  .$feature->end            .","
                  .$feature->score          .","
                  .$feature->strand         .","
                  ."'".$feature->source_tag."',"
                  .$analysisid              .","
                  .$homol->start            .","
                  .$homol->end              .","
                  ."'".$homol->seqname      ."',"
                  .((defined $feature->percent_id)   ? $feature->percent_id  : 'NULL')  .","    
                  .((defined $feature->p_value)      ? &exponent($feature->p_value)     : 'NULL')  .","
                  .((defined $feature->phase)        ? $feature->phase       : 'NULL')  .","
                  .((defined $feature->end_phase)    ? $feature->end_phase   : 'NULL')  .")");
        
        $sth->execute();
        
        
        } else {
		my $sth = $self->_db_obj->prepare(  "INSERT into ".$self->_lookup_table_name('feature')."(id,contig,seq_start,seq_end,score,strand,name,analysis,hstart,hend,hid,perc_id,evalue,phase,end_phase) ".
                                        "values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)");
    
        $sth->execute('NULL',
			          $contig->internal_id,
			          $feature->start,
			          $feature->end,
			          $feature->score,
			          $feature->strand,
			          $feature->source_tag,
			          $analysisid,
			          -1,
			          -1,
			          "__NONE__",
                      'NULL',
                      'NULL',
                      'NULL',
                      'NULL' );
	    }
	}
    }

    my $sth2 = $self->_db_obj->prepare("INSERT into ".$self->_lookup_table_name('repeat_feature')."(id,contig,seq_start,seq_end,score,strand,analysis,hstart,hend,hid) values(?,?,?,?,?,?,?,?,?,?)");

    foreach my $feature (@repeats) {
	if (!defined($feature->analysis)) {
	    $self->throw("Feature " . $feature->seqname . " " . $feature->source_tag ." doesn't have analysis. Can't write to database");
	} else {
	    $analysis = $feature->analysis;
	}


	my $analysisid = $self->write_Analysis($analysis);
	my $homol = $feature->feature2;

	$sth2->execute('NULL',
		       $contig->internal_id,
		       $feature->start,
		       $feature->end,
		       $feature->score,
		       $feature->strand,
		       $analysisid,
		       $homol->start,
		       $homol->end,
		       $homol->seqname);
    }


    # Now the predictions
    # we can't block do these as we need to get out the id wrt to the features
    foreach my $feature ( @fset ) {
#	print STDERR "Adding in a fset feature ",$feature->gff_string,"\n";

	if (!defined($feature->analysis)) {

	    $self->throw("Feature " . $feature->seqname . " " . 
			              $feature->source_tag . 
			 " doesn't have analysis. Can't write to database");
	} else {
	    $analysis = $feature->analysis;
	}

	my $analysisid = $self->write_Analysis($analysis);
	my $score      = $feature->score();

	if( !defined $score ) { $score = "-1000"; }

	my $sth3 = $self->_db_obj->prepare("INSERT into ".$self->_lookup_table_name('fset')."(id,score) values ('NULL',$score)");
	   $sth3->execute();

	# get out this id. This looks really clunk I know. Any better ideas... ?

	my $sth4 = $self->_db_obj->prepare("SELECT LAST_INSERT_ID()");
	   $sth4->execute();

	my $arr = $sth4->fetchrow_arrayref();
	my $fset_id = $arr->[0];

	# now write each sub feature
	my $rank = 1;

	foreach my $sub ( $feature->sub_SeqFeature ) {
	    my $sth5 = $self->_db_obj->prepare("INSERT into ".$self->_lookup_table_name('feature')."(id,contig,seq_start,seq_end,score,strand,analysis,name,hstart,hend,hid,evalue,perc_id,phase,end_phase) "
                      ."values('NULL','"
                      .$contig->internal_id ."',"
				      .$sub->start          .","
				      .$sub->end            . ","
				      .$sub->score          . ","
				      .$sub->strand         . ","
				      .$analysisid          . "," 
				      ."'".$sub->source_tag ."',"
                      ."-1,-1,"
                      ."'".($sub->primary_tag || "__NONE__")."',"
                      . ((defined $sub->p_value)     ?   &exponent($sub->p_value)       : 'NULL')  .","
                      . ((defined $sub->percent_id)  ?   $sub->percent_id    : 'NULL')  .","
                      . ((defined $sub->phase)       ?   $sub->phase         : 'NULL')  .","
                      . ((defined $sub->end_phase)   ?   $sub->end_phase     : 'NULL')  .")");
	    
        $sth5->execute();
	    my $sth6 = $self->_db_obj->prepare("INSERT into ".$self->_lookup_table_name('fset_feature')."(fset,feature,rank) values ($fset_id,LAST_INSERT_ID(),$rank)");
	    $sth6->execute();
	    $rank++;
	    }
    }
    
    return 1;
}

=head2 get_Protein_annseq

 Title   : get_Protein_annseq
 Usage   : get_Protein_annseq ($ENSP); 
 Function: Creates an annseq object for a particular peptide, storing the peptide
           sequence in $annseq->primary_seq, and adding all the protein features as generic
           Seqfeatures
 Example : 
 Returns : $annseq
 Args    : $ENSP


=cut

sub get_Protein_annseq{
    my ($self,$ENSP) = @_;

    my $annseq = Bio::EnsEMBL::AnnSeq->new();
    
    my $sth     = $self->_db_obj->prepare("SELECT id from ".$self->_lookup_table_name('transcript')." where translation = '$ENSP'");
    my $res     = $sth->execute();
    my $rowhash = $sth->fetchrow_hashref;
    
    my $gene_obj=Bio::EnsEMBL::DBOLD::Gene_Obj->new($self->_db_obj);
    my $transcript = $gene_obj->get_Transcript($rowhash->{'id'});
    my $translation = $gene_obj->get_Translation($ENSP);

    $transcript->translation($translation);

    my $seq = $transcript->translate();
    $annseq->primary_seq($seq);

    $sth = $self->_db_obj->prepare("SELECT * from ".$self->_lookup_table_name('proteinfeature')." where translation = '$ENSP'");
    $res = $sth->execute();

    while( my $rowhash = $sth->fetchrow_hashref) {
	my $analysis = $rowhash->{'analysis'};
	my $sth2     = $self->_db_obj->prepare("SELECT * from ".$self->_lookup_table_name('analysis')." where id = '$analysis'");
	my $res2     = $sth2->execute();
	my $rowhash2 = $sth2->fetchrow_hashref;

	my $feature  = new Bio::SeqFeature::Generic ( -start   => $rowhash->{'seq_start'}, 
						      -end     => $rowhash->{'seq_end'},
						      -score   => $rowhash->{'score'},
						      -primary => $rowhash2->{'gff_feature'},
						      -source  => $rowhash2->{'gff_source'});

	$annseq->add_SeqFeature($feature);
    }
    
    return $annseq;   
}

=head2 write_all_Protein_features

 Title   : write_all_Protein_features
 Usage   : $obj->write_all_Protein_features($ENSP)
 Function: writes all protein features of a particular peptide into the database          
 Example :
 Returns : 
 Args    :


=cut

sub write_all_Protein_features {
    my ($self,$prot_annseq,$ENSP) = @_;
    
    my $c=0;
    foreach my $feature ($prot_annseq->all_SeqFeatures()) {
	my $sth = $self->_db_obj->prepare("INSERT into ".$self->_lookup_table_name('proteinfeature')." (id,seq_start, seq_end, score, analysis, translation) values (NULL,"
				 .$feature->start().","
				 .$feature->end().","
				 .$feature->score().",'"
				 .$c."','"
				 .$ENSP."')");
	$sth->execute();
	
	my $sth2 = $self->_db_obj->prepare("INSERT into ".$self->_lookup_table_name('analysis')." (id,db,db_version,program,program_version,gff_source,gff_feature) values ('$c','testens',1,'elia_program',1,'"
				  .$feature->source_tag()."','"
				  .$feature->primary_tag()."')");
	 $sth2->execute();
	$c++;
    }
}

=head2 write_Protein_feature

 Title   : write_Protein_feature
 Usage   : $obj->write_Protein_feature($ENSP, $feature)
 Function: writes a protein feature object of a particular peptide into the database          
 Example :
 Returns : 
 Args    :


=cut

sub write_Protein_feature {
    my ($self,$ENSP,$feature) = @_;
    
    my $sth = $self->_db_obj->prepare("INSERT into ".$self->_lookup_table_name('proteinfeature')." (seq_start, seq_end, score, translation) values ("
			     .$feature->start()." ,"
			     .$feature->end()." ,'"
			     .$feature->score()." ,'"
			     .$ENSP."'
				)");
    $sth->execute();
}

=head2 get_Analysis

 Title   : get_Analysis
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Analysis {
    my ($self,$id) = @_;

    my $sth = $self->_db_obj->prepare("SELECT db,db_version,program,program_version,gff_source,gff_feature,id from ".$self->_lookup_table_name('analysis')." where id = $id");
    my $rv  = $sth->execute;
    my $rh  = $sth->fetchrow_hashref;

    if ($sth->rows) {
	my $anal = Bio::EnsEMBL::FeatureFactory->new_analysis();

	if( defined $rh->{'db'} ) {
	    $anal->db($rh->{'db'});
	}
	if( defined $rh->{'db_version'} ) {
	    $anal->db_version($rh->{'db_version'});
	}

	$anal->program        ($rh->{'program'});
	$anal->program_version($rh->{'program_version'});
	$anal->gff_source     ($rh->{gff_source});
	$anal->gff_feature    ($rh->{gff_feature});
	my $mid = $rh->{'id'};

	$anal->dbID("$mid");
	return $anal;
    }  else {
	$self->throw("Can't fetch analysis id $id\n");
    }
    
}

=head2 exists_Analysis

 Title   : exists_Analysis
 Usage   : $obj->exists_Analysis($anal)
 Function: Tests whether this feature already exists in the database
 Example :
 Returns : Analysis id if the entry exists
 Args    : Bio::EnsEMBL::Analysis


=cut

sub exists_Analysis {
    my ($self,$anal) = @_;

    
    $self->throw("Object is not a Bio::EnsEMBL::AnalysisI") unless $anal->isa("Bio::EnsEMBL::AnalysisI");
    # If all the attributes of the analysis object are not set it's existence can't be tested 
    $self->throw("program property of analysis object not defined") unless ($anal->program); 
    $self->throw("program_version property of analysis object not defined") unless ($anal->program_version);    
    $self->throw("gff_source property of analysis object not defined") unless ($anal->gff_source); 
    $self->throw("gff_feature property of analysis object not defined") unless ($anal->gff_feature);   
        
    my $query;

    if ($anal->has_database == 1) {
            $query = "SELECT id from ".$self->_lookup_table_name('analysis')." where db = \""      . $anal->db              . "\" and" .
                " db_version = \""      . $anal->db_version      . "\" and " .
                " program =    \""      . $anal->program         . "\" and " .
                " program_version = \"" . $anal->program_version . "\" and " .
                " gff_source = \""      . $anal->gff_source      . "\" and" .
                " gff_feature = \""     . $anal->gff_feature     . "\"";
    } else {
        $query = "SELECT id from ".$self->_lookup_table_name('analysis')." where " .
                " program =    \""      . $anal->program         . "\" and " .
                " program_version = \"" . $anal->program_version . "\" and " .
                " gff_source = \""      . $anal->gff_source      . "\" and" .
                " gff_feature = \""     . $anal->gff_feature     . "\"";
    }
    
    if( exists $self->_db_obj->_analysis_cache->{$query} ) {
	return $self->_db_obj->_analysis_cache->{$query};
    }

    my $sth = $self->_db_obj->prepare($query);
    my $rv = $sth->execute();

    if ($rv && $sth->rows > 0) {
	my $rowhash = $sth->fetchrow_hashref;
	my $anaid = $rowhash->{'id'}; 
	$self->_db_obj->_analysis_cache->{$query} = $anaid;
	return $anaid;
    } else {
	return 0;
    }
}

=head2 write_Analysis

 Title   : write_Analysis
 Usage   : $obj->write_Analysis($anal)
 Function: Writes analysis details to the database
           Checks first whether this analysis entry already exists
 Example :
 Returns : int
 Args    : Bio::EnsEMBL::AnalysisI

=cut

sub write_Analysis {
    my ($self,$anal) = @_;

    $self->throw("Argument is not a Bio::EnsEMBL::AnalysisI") unless $anal->isa("Bio::EnsEMBL::AnalysisI");
    # If all the attributes of the analysis object are not set it shouldn't be written
    $self->throw("program property of analysis object not defined") unless ($anal->program); 
    $self->throw("program_version property of analysis object not defined") unless ($anal->program_version);    
    $self->throw("gff_source property of analysis object not defined") unless ($anal->gff_source); 
    $self->throw("gff_feature property of analysis object not defined") unless ($anal->gff_feature); 
        
    # First check whether this entry already exists.
    my $query;
    my $analysisid = $self->exists_Analysis($anal);
    return $analysisid if $analysisid;

    
    if ($anal->has_database == 1) {
        local $^W = 0;
	$query = "INSERT into ".$self->_lookup_table_name('analysis')."(id,db,db_version,program,program_version,gff_source,gff_feature) values (NULL,\"" .
                $anal->db               . "\",\""   .
                $anal->db_version       . "\",\""   .
		$anal->program          . "\",\"" .
		$anal->program_version  . "\",\"" .
                $anal->gff_source       . "\",\"" .
                $anal->gff_feature      . "\")";
    } else {
	$query = "INSERT into ".$self->_lookup_table_name('analysis')."(id,program,program_version,gff_source,gff_feature) values (NULL,\"" .
                $anal->program          . "\",\"" .
                $anal->program_version  . "\",\"" .
                $anal->gff_source       . "\",\"" .
                $anal->gff_feature      . "\")";
    }

    my $sth  = $self->_db_obj->prepare($query);
    my $rv   = $sth->execute;
    
    
    $sth = $self->_db_obj->prepare("SELECT last_insert_id()");
    $rv  = $sth->execute;
    
    if ($sth->rows == 1) {
	my $rowhash = $sth->fetchrow_hashref;
	return $rowhash->{'last_insert_id()'};
    } else {
	$self->throw("Wrong number of rows returned : " . $sth->rows . " : should be 1");
    }

}

=head2 find_GenomeHits 

 Title   : find_GenomeHits
 Usage   : $obj->find_GenomeHits($hitid)
 Function: 
 Example : 
 Returns : 
 Args    : 


=cut


sub find_GenomeHits {
    my ($self,$arg) = @_;

    $self->throw("No hit id input") unless defined($arg);

    my $query = "SELECT c.id, " .
	                "f.seq_start, " . 
			"f.seq_end, "   . 
			"f.score, "     .
			"f.strand, "    .
			"f.analysis, "  .
			"f.name, "      .
			"f.hstart, "    .
			"f.hend, "      .
			"f.hid "       .
	        "from   ".$self->_lookup_table_name('feature')." as f,".$self->_lookup_table_name('contig')." as c " .
		"where  f.hid = '$arg' and " . 
		        "c.internal_id = f.contig";

    my $sth   = $self->_db_obj->prepare($query);
    my $res   = $sth->execute;
    
    my ($contig,$start,$end,$score,$strand,$analysisid,$name,$hstart,$hend,$hid);
    

    $sth->bind_columns(undef,\$contig,\$start,\$end,\$score,\$strand,\$analysisid,
		       \$name,\$hstart,\$hend,\$hid);
    

    my %analhash;         # Stores all the analysis objects
    my @features;

    while($sth->fetch) {
	my $out;
	my $analysis;
	
	if (!$analhash{$analysisid}) {
	   
	    $analysis = $self->get_Analysis($analysisid);
	    $analhash{$analysisid} = $analysis;
	   
	} else {
	    $analysis = $analhash{$analysisid};
	}
       
	if( !defined $name ) {
	    $name = 'no_source';
	}
	 
	$out = Bio::EnsEMBL::FeatureFactory->new_feature_pair();
	$out->set_all_fields($start,$end,$strand,$score,$name,'similarity',$contig,
			     $hstart,$hend,1,$score,$name,'similarity',$hid);

	$out->analysis($analysis);
	$out->validate;
       
      push(@features,$out);
	
    }

    return @features;
}



=head2 get_PredictionFeature_by_id 

 Title   : get_PredictionFeature_by_id
 Usage   : $obj->get_PredictionFeature_by_id($id)
 Function: 
 Example : 
 Returns : 
 Args    : 


=cut






sub get_PredictionFeature_by_id {
   my ($self,$genscan_id) = @_;
 
   unless ($genscan_id){$self->throw("I need a genscan id");}
   
   my $fsetid;
   my %analhash;
  
   my $query = "SELECT f.id,f.seq_start,f.seq_end,f.strand,f.score,f.analysis,fset.id,c.id " .
       "from ".$self->_lookup_table_name('feature')." f, ".$self->_lookup_table_name('fset')." fset,".$self->_lookup_table_name('fset_feature')." ff,".$self->_lookup_table_name('contig')." c where ff.feature = f.id and fset.id = ff.fset ".
        " and c.internal_id=f.contig and ff.fset ='$genscan_id' and name = 'genscan'";
 
       
   my $sth = $self->_db_obj->prepare($query);   
        
   $sth->execute();
   
   my ($fid,$start,$end,$strand,$score,$analysisid,$contig);
           
    $sth->bind_columns(undef,\$fid,\$start,\$end,\$strand,\$score,\$analysisid,\$fsetid,\$contig);

   my $current_fset;
   if( $sth->fetch ) {
       my $out;

       my $analysis;

       if (!$analhash{$analysisid}) {
       
           
           $analysis = $self->get_Analysis($analysisid);

           $analhash{$analysisid} = $analysis;
       
       } else {
           $analysis = $analhash{$analysisid};
       }
          
           $current_fset = new Bio::EnsEMBL::SeqFeature;
           $current_fset->source_tag('genscan');
           $current_fset->primary_tag('prediction');
           $current_fset->analysis($analysis);
           $current_fset->seqname($contig);
           $current_fset->id($fsetid);
        

   $out = new Bio::EnsEMBL::SeqFeature;
 
       $out->seqname   ($contig);
       $out->start     ($start);
       $out->end       ($end);
       $out->strand    ($strand);  
 
       $out->source_tag('genscan');
       $out->primary_tag('prediction');
 
       if( defined $score ) {
           $out->score($score);
       }

       $out->analysis($analysis);

       # Final check that everything is ok.
           
       $out->validate();
       $current_fset->add_sub_SeqFeature($out,'EXPAND');
       $current_fset->strand($strand);
   }
   
     else { $self->throw("Fset $genscan_id does not exist in the database");}


return $current_fset;

}




=head2 get_PredictionFeature_as_Transcript 

 Title   : get_PredictionFeature_as_Transcript
 Usage   : $obj->get_PredictionFeature_as_Transcript($id)
 Function: 
 Example : 
 Returns : 
 Args    : 


=cut





sub get_PredictionFeature_as_Transcript{
    my ($self,$genscan_id)=@_;
    unless ($genscan_id){$self->throw("I need a genscan id");}

    my $ft=$self->get_PredictionFeature_by_id($genscan_id);
    my $contig=$self->_db_obj->get_Contig($ft->seqname);
 
    &Bio::EnsEMBL::DBOLD::Utils::fset2transcript($ft,$contig);


}



=head2 _db_obj

 Title   : _db_obj
 Usage   : $obj->_db_obj($newval)
 Function: 
 Example : 
 Returns : value of _db_obj
 Args    : newvalue (optional)


=cut

sub _db_obj{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_db_obj'} = $value;
    }
    return $self->{'_db_obj'};

}

sub exponent {
    my ($number) = @_;

    my ($exp) = sprintf("%.3e", $number);
    return $exp;
}
