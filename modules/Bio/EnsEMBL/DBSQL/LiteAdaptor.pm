# EnsEMBL Gene reading writing adaptor for mySQL
#
# Copyright EMBL-EBI 2001
#
# Author: James Smith
#

=head1 NAME

Bio::EnsEMBL::DBSQL::LiteAdaptor - MySQL Database queries to generate and store gens.

=head1 SYNOPSIS

=head1 CONTACT

  Arne Stabenau: stabenau@ebi.ac.uk
  James Smith  : js5@sanger.ac.uk

=head1 APPENDIX

=cut


package Bio::EnsEMBL::DBSQL::LiteAdaptor;
use strict;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::Annotation::DBLink;
use vars '@ISA';

@ISA = ( 'Bio::EnsEMBL::DBSQL::BaseAdaptor' );

sub new {
    my ($class,$dbobj) = @_;

    my $self = {};
    bless $self,$class;

    if( !defined $dbobj || !ref $dbobj ) {
        $self->throw("Don't have a db [$dbobj] for new adaptor");
    }

    $self->db($dbobj);

    $self->{'_lite_db_name'} = $dbobj->{'_lite_db_name'};
    return $self;
}

#sub db {
#    my ($self, $arg) = @_;
#    $self->{'_db'} = $arg if ($arg);
#    return($self->{'_db'});
#}

sub fetch_virtualtranscripts_start_end {
    my ( $self, $chr, $vc_start, $vc_end, $database ) =@_;
    my $_db_name = $self->{'_lite_db_name'};
    $database      ||= '';
    my $cache_name = "_$database"."_vtrans_cache_$chr"."_$vc_start"."_$vc_end";
    return $self->{$cache_name} if( $self->{$cache_name} );
    my $sth = $self->prepare(
        "select *
           from $_db_name.transcript
          where chr_name = ? and chr_start <= ? and chr_start >= ? and
                chr_end >= ?".
            ( $database ne '' ? " and db = '$database'" : "" )
    );
    
    eval {
        $sth->execute( "$chr", $vc_end, $vc_start-3000000, $vc_start );
    };
    return [] if($@);
    my @transcripts;
    while( my $row = $sth->fetchrow_hashref() ) {
        push @transcripts, {
            'db'            => $row->{'db'},
            'type'          => $row->{'type'},
            'transcript'    => $row->{'transcript_id'},
            'stable_id'     => $row->{'transcript_name'},
            'translation'   => $row->{'translation_name'},
            'gene'          => $row->{'gene_name'},
            'chr_name'      => $row->{'chr_name'},
            'chr_start'     => $row->{'chr_start'},
            'chr_end'       => $row->{'chr_end'},
            'coding_start'  => ( $row->{'coding_start'}||$row->{'chr_start'} ) - $vc_start + 1 ,
            'coding_end'    => ( $row->{'coding_end'}  ||$row->{'chr_end'}   ) - $vc_start + 1 ,
            'strand'        => $row->{'chr_strand'},
            'start'         => $row->{'chr_start'} - $vc_start + 1,
            'end'           => $row->{'chr_end'}   - $vc_start + 1,
            'external_db'   => $row->{'external_db'},
            'synonym'       => $row->{'external_name'},
            'exon_structure' => [ split ':', $row->{'exon_structure'} ],
        };
    }
    return $self->{$cache_name} = \@transcripts;
    return \@transcripts
}

sub fetch_virtualgenscans_start_end {
    my ( $self, $chr, $vc_start, $vc_end ) =@_;
    my $_db_name = $self->{'_lite_db_name'};
    my $cache_name = "_virtualgenscans_cache_$chr"."_$vc_start"."_$vc_end";
    
    return $self->{$cache_name} if( $self->{$cache_name} );
    my $sth = $self->prepare(
        "select *
           from $_db_name.genscan
          where chr_name = ? and chr_start <= ? and chr_start >= ? and
                chr_end >= ?"
    );
    eval {
        $sth->execute( "$chr", $vc_end, $vc_start-2000000, $vc_start );
    };
    return [] if($@);
    my @transcripts;
    while( my $row = $sth->fetchrow_hashref() ) {
        push @transcripts, {
            'genscan'   => $row->{'name'},
            'chr_name'  => $row->{'chr_name'},
            'chr_start' => $row->{'chr_start'},
            'chr_end'   => $row->{'chr_end'},
            'start'     => $row->{'chr_start'} - $vc_start + 1,
            'end'       => $row->{'chr_end'}   - $vc_start + 1,
            'strand'    => $row->{'chr_strand'},
            'exon_structure' => [ split ':', $row->{'exon_structure'} ]
        };
    }
    return $self->{$cache_name} = \@transcripts;
    return \@transcripts
}

sub fetch_virtualgenes_start_end {
    my ( $self, $chr, $vc_start, $vc_end, $database ) =@_;
    my $_db_name = $self->{'_lite_db_name'};
    $database      ||= '';
    my $cache_name = "_virtualgenes_$database"."_cache_$chr"."_$vc_start"."_$vc_end";
    return $self->{$cache_name} if( $self->{$cache_name} );
    my $sth = $self->prepare(
        "select *
           from $_db_name.gene as g 
          where chr_name = ? and g.chr_start <= ? and g.chr_start >= ? and
                g.chr_end >= ?".
            ( $database ne '' ? " and db = '$database'" : "" )
    );
    eval {
        $sth->execute( "$chr", $vc_end, $vc_start-3000000, $vc_start );
    };
    return [] if($@);
    my @genes;
    while( my $row = $sth->fetchrow_hashref() ) {
        push @genes, {
            'db'        => $row->{'db'},
            'type'      => $row->{'type'},
            'gene'      => $row->{'gene_id'},
            'stable_id' => $row->{'gene_name'},
            'chr_name'  => $row->{'chr_name'},
            'chr_start' => $row->{'chr_start'},
            'chr_end'   => $row->{'chr_end'},
            'start'     => $row->{'chr_start'} - $vc_start + 1,
            'end'       => $row->{'chr_end'}   - $vc_start + 1,
            'strand'    => $row->{'chr_strand'},
            'external_db' => $row->{'external_db'},
            'synonym'     => $row->{'external_name'},
        };
    }
    return $self->{$cache_name} = \@genes;
    return \@genes
}
                
# if count is -ve, get genes before start, else after.
sub fetch_virtualgenes_start_count {
    my ( $self, $chr, $vc_start, $count ) =@_;
    my $_db_name = $self->{'_lite_db_name'};
    my $cache_name = "_virtualgenes_count_cache_$chr"."_$vc_start"."_$count";
    return $self->{$cache_name} if( $self->{$cache_name} );

    my $sql;
    if ($count < 0){
	$sql =  "select g.gene_id, g.gene_name, 
                g.chr_name, g.chr_start, g.chr_end,
                g.chr_strand, g.external_name, g.external_db
           from $_db_name.gene as g 
          where g.chr_name = ? and g.chr_start < ?
		and g.db = 'core'
	  order by g.chr_start desc
          limit ?"
    }
    else {
	$sql = "select g.gene_id, g.gene_name, 
                g.chr_name, g.chr_start, g.chr_end,
                g.chr_strand, g.external_name, g.external_db
           from $_db_name.gene as g 
          where g.chr_name = ? and g.chr_start >= ?
		and g.db = 'core'
	  order by g.chr_start
          limit ?"
    }
    
    my $sth = $self->prepare( $sql );
    eval {
        $sth->execute( "$chr", $vc_start, abs($count) );
    };
    return [] if($@);
    my @genes;
    while( my $row = $sth->fetchrow_arrayref() ) {
        push @genes, {
            'gene'      => $row->[0],
            'stable_id' => $row->[1],
            'chr_name'  => $row->[2],
            'chr_start' => $row->[3],
            'chr_end'   => $row->[4],
            'start'     => $row->[3]-$vc_start+1,
            'end'       => $row->[4]-$vc_start+1,
            'strand'    => $row->[5],
            'synonym'   => $row->[6],
            'db'        => $row->[7]
        };
    }
    @genes = reverse @genes if $count < 0;
    return $self->{$cache_name} = \@genes;
    return \@genes
}

sub fetch_virtualRepeatFeatures_start_end {
    my ( $self, $chr, $vc_start, $vc_end, $type, $glob_bp ) =@_;
    $type ||= '';
    my $cache_name = "_repeats_$type"."_cache_$chr"."_$vc_start"."_$vc_end";
    return $self->{$cache_name} if( $self->{$cache_name} );
	$glob_bp ||= 0;
    my $_db_name = $self->{'_lite_db_name'};

    my $sth = $self->prepare(
        "select *
           from $_db_name.repeat
          where chr_name = ? and chr_start <= ? and chr_start >= ? and chr_end >= ?".
		  	( (defined $type && $type ne '') ? " and type = '$type'" : '' ).
          " order by chr_start"            
    );

    eval {
        $sth->execute( "$chr", $vc_end, $vc_start-1000000, $vc_start);
    };
    return [] if($@);

	my @repeats;
	my $old_start = -99999999999999999;
	my $old_end   = -99999999999999999;
	while( my $row = $sth->fetchrow_hashref() ) {
      	my $end = $row->{'chr_end'};
## Glob results! 
        next if($end < $old_end );
    	$old_end   = $end;
    	if( $end-$old_start < $glob_bp/2 ) {
			$repeats[-1]->{'chr_end'} = $end; 
			$repeats[-1]->{'end'}     = $end - $vc_start + 1; 
	  	}	else {
    	  	$old_start = $row->{'chr_start'};
			push @repeats, {
				'id'        => $row->{'id'},
				'hid'       => $row->{'hid'},
    	        'chr_name'  => $row->{'chr_name'},
        	    'chr_start' => $old_start,
            	'chr_end'   => $end,
        	    'start'     => $old_start-$vc_start+1,
            	'end'       => $end      -$vc_start+1,
	            'strand'    => $row->{'chr_strand'},
			};
		}
    }
    return $self->{$cache_name} = \@repeats;
    return \@repeats;
}


sub fetch_snp_features {

 my ( $self, $chr, $vc_start, $vc_end,$glob ) =@_;


 #lists of variations to be returned
    my @variations;
    my %hash;
    my $string; 

    my $_db_name = $self->{'_lite_db_name'};
   
    my $query = qq{

        SELECT   chr_start, chr_strand,
                 refsnpid, tscid, hgbaseid
        FROM   	 $_db_name.snp
        WHERE  	 chr_name='$chr' 
        AND      chr_start>$vc_start
    	AND      chr_start<$vc_end
              };

    #&eprof_start('snp-sql-query');

    my $sth = $self->prepare($query);

    eval {
        $sth->execute( );
    };
    return () if($@);
    #&eprof_end('snp-sql-query');

    my $snp;
    my $cl;

    #&eprof_start('snp-sql-object');

  SNP:
    while( (my $arr = $sth->fetchrow_arrayref()) ) {
        
        my ($snp_start, $chrom_strand,$snpuid,$tscid, $hgbaseid) = @{$arr};
            
  # globbing
        
        my $key=$snpuid; #.$acc;           # for purpose of filtering duplicates
        my %seen;                       # likewise
        
        
        if ( ! $seen{$key} )  {
            ## we're grabbing all the necessary stuff from the db in one
            ## SQL statement for speed purposes, so we have to do some
            ## duplicate filtering here.

            $seen{$key}++;
            
            #Variation
            $snp = new Bio::EnsEMBL::ExternalData::Variation
              (-start => $snp_start-$vc_start +1 ,
               -end => $snp_start-$vc_start +1,
               -strand => $chrom_strand,
#               -original_strand => $strand,
               -score => 1,
               -source_tag => 'dbSNP',
              );
            
            my $link = new Bio::Annotation::DBLink;
            $link->database('dbSNP');
            $link->primary_id($snpuid);
            #add dbXref to Variation
            $snp->add_DBLink($link);
	    if ($hgbaseid) {
	      my $link2 = new Bio::Annotation::DBLink;
	      $link2->database('HGBASE');
	      $link2->primary_id($hgbaseid);
	      $snp->add_DBLink($link2);
	    }
	    if ($tscid) {
	      my $link3 = new Bio::Annotation::DBLink;
	      $link3->database('TSC-CSHL');
	      $link3->primary_id($tscid);
	      #add dbXref to Variation
	      $snp->add_DBLink($link3);
	    }
#            $cl=$acc;
            # set for compatibility to Virtual Contigs
            #add SNP to the list
            push(@variations, $snp);
        }                               # if ! $seen{$key}
      }                                    # while a row from select statement

    #&eprof_end('snp-sql-object');
    
    return @variations;

}

sub fetch_virtualfeatures {
    my ( $self, $chr, $vc_start, $vc_end, $type, $score, $glob ) =@_;
    my $_db_name = $self->{'_lite_db_name'};
    my $cache_name = "_$type"."_cache_$chr"."_$vc_start"."_$vc_end"."_$score";
    return $self->{$cache_name} if( $self->{$cache_name} );
    my $sth = $self->prepare(
        "select *
           from $_db_name.$type
          where chr_name=? and chr_start<=? and chr_start >= ? and chr_end >= ? and
                score >= ?"
    );
    eval {
        $sth->execute( "$chr", $vc_end, $vc_start-1000000, $vc_start, $score );
    };
    return [] if($@);
    my @features;
    while( my $row = $sth->fetchrow_hashref() ) {
        push @features, {
            'chr_name'  => $row->{'chr_name'},
            'chr_start' => $row->{'chr_start'},
            'chr_end'   => $row->{'chr_end'},
            'start'     => $row->{'chr_start'} - $vc_start + 1,
            'end'       => $row->{'chr_end'} - $vc_start + 1,
            'strand'    => $row->{'chr_strand'},
            'id'        => $row->{'id'},
            'score'     => $row->{'score'}
        };
    }
    return $self->{$cache_name} = \@features;
}
    
sub fetch_virtualsnps {
    my ( $self, $chr, $vc_start, $vc_end, $glob_bp ) =@_;
    $glob_bp||=0;
    my $_db_name = $self->{'_lite_db_name'};
    my $cache_name = "_snp_cache_$chr"."_$vc_start"."_$vc_end";
    return $self->{$cache_name} if( $self->{$cache_name} );
    my $sth = $self->prepare(
        "select *
           FROM $_db_name.snp
          WHERE chr_name = ? AND      chr_start>=? AND      chr_start<=?
        order by chr_start"
    );
    eval {
        $sth->execute( "$chr", $vc_start, $vc_end );
    };
    return [] if($@);
    my @variations;
	my $old_start = -99999999999999999;
	while( my $row = $sth->fetchrow_hashref() ) {
      	my $start = $row->{'chr_start'};
## Glob results! 
        next if($start < $old_start );
    	if($start < $old_start + $glob_bp/2) {
			$variations[-1]->{'end'}     = $start - $vc_start + 1;
			$variations[-1]->{'chr_end'} = $start;
	  	}	else {
            push @variations, {
                'chr_name'  => $row->{'chr_name'},
                'chr_start' => $start,
                'chr_end'   => $start,
                'start'     => $start - $vc_start + 1,
                'end'       => $start - $vc_start + 1,
                'strand'    => $row->{'chr_strand'},
                'id'        => $row->{'refsnpid'},
                'tscid'     => $row->{'tscid'},
                'hgbaseid'  => $row->{'hgbaseid'},
                'anosnpid'  => $row->{'anosnpid'},
                'type'      => $row->{'type'},
            };
    	  	$old_start = $start;
		}
    }

    return $self->{$cache_name} = \@variations;
}

1;
__END__

