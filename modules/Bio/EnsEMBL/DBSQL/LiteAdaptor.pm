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
    $database      ||= 'ensembl';
    my $cache_name = "_$database"."_vtrans_cache_$chr"."_$vc_start"."_$vc_end";
    return $self->{$cache_name} if( $self->{$cache_name} );
    my $sth = $self->prepare(
        "select transcript_id, transcript_name, translation_name, gene_name,
                chr_start, chr_end, chr_strand, external_name, external_db,
                exon_structure, type
           from $_db_name.js5_transcript
          where chr_name = ? and chr_start <= ? and chr_start >= ? and
                chr_end >= ? and db = ?"
    );
    
    eval {
        $sth->execute( $chr, $vc_end, $vc_start-3000000, $vc_start, $database );
    };
    return [] if($@);
    my @transcripts;
    while( my $row = $sth->fetchrow_arrayref() ) {
        push @transcripts, {
            'transcript'=> $row->[0],
            'stable_id' => $row->[1],
            'translation'=> $row->[2],
            'gene'      => $row->[3],
            'chr_start' => $row->[4],
            'chr_end'   => $row->[5],
            'start'     => $row->[4]-$vc_start+1,
            'end'       => $row->[5]-$vc_start+1,
            'strand'    => $row->[6],
            'synonym'   => $row->[7],
            'db'        => $row->[8],
            'exon_structure' => [ split ':', $row->[9] ],
            'type'      => $row->[10]
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
        "select id, chr_name, chr_start, chr_end, chr_strand, exon_structure
           from $_db_name.js5_genscan_3
          where chr_name = ? and chr_start <= ? and chr_start >= ? and
                chr_end >= ?"
    );
    eval {
        $sth->execute( $chr, $vc_end, $vc_start-1000000, $vc_start );
    };
    return [] if($@);
    my @transcripts;
    while( my $row = $sth->fetchrow_arrayref() ) {
        push @transcripts, {
            'genscan'   => $row->[0],
            'chr_start' => $row->[2],
            'chr_end'   => $row->[3],
            'start'     => $row->[2]-$vc_start+1,
            'end'       => $row->[3]-$vc_start+1,
            'strand'    => $row->[4],
            'exon_structure' => [ split ':', $row->[5] ]
        };
    }
    return $self->{$cache_name} = \@transcripts;
    return \@transcripts
}

sub fetch_virtualgenes_start_end {
    my ( $self, $chr, $vc_start, $vc_end ) =@_;
    my $_db_name = $self->{'_lite_db_name'};
    my $cache_name = "_virtualgenes_cache_$chr"."_$vc_start"."_$vc_end";
    return $self->{$cache_name} if( $self->{$cache_name} );
    my $sth = $self->prepare(
        "select g.gene, g.name, 
                g.chr_name, g.gene_chrom_start, g.gene_chrom_end,
                g.chrom_strand, gx.display_id, gx.db_name
           from $_db_name.gene as g, $_db_name.gene_xref as gx
          where g.gene = gx.gene and
                g.chr_name = ? and g.gene_chrom_start <= ? and g.gene_chrom_start >= ? and
                g.gene_chrom_end >= ?"
    );
    eval {
        $sth->execute( $chr, $vc_end, $vc_start-1000000, $vc_start );
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
    return $self->{$cache_name} = \@genes;
    return \@genes
}
                
sub fetch_EMBLgenes_start_end {
    my ( $self, $chr, $vc_start, $vc_end ) =@_;
    my $_db_name = $self->{'_lite_db_name'};
    my $cache_name = "_emblgenes_cache_$chr"."_$vc_start"."_$vc_end";
    return $self->{$cache_name} if( $self->{$cache_name} );
    my $sth = $self->prepare(
        "select g.gene, g.name, 
                g.chr_name, g.gene_chrom_start, g.gene_chrom_end,
                g.chrom_strand, gx.display_id, gx.db_name, g.type
           from $_db_name.embl_gene as g, $_db_name.embl_gene_xref as gx
          where g.gene = gx.gene and
                g.chr_name = ? and g.gene_chrom_start <= ? and g.gene_chrom_start >= ? and
                g.gene_chrom_end >= ?"
    );
    eval {
        $sth->execute( $chr, $vc_end, $vc_start-1000000, $vc_start );
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
            'db'        => $row->[7],
            'type'      => $row->[8]
        };
    }
    return $self->{$cache_name} = \@genes;
    return \@genes
}

sub fetch_SangerGenes_start_end {
    my ( $self, $chr, $vc_start, $vc_end ) =@_;
    my $_db_name = $self->{'_lite_db_name'};
    my $cache_name = "_sangergenes_cache_$chr"."_$vc_start"."_$vc_end";
    return $self->{$cache_name} if( $self->{$cache_name} );
    my $sth = $self->prepare(
        "select g.gene, g.name,
                g.chr_name, g.gene_chrom_start, g.gene_chrom_end,
                g.chrom_strand, gx.display_id, gx.db_name, g.type
           from $_db_name.sanger_gene as g, $_db_name.sanger_gene_xref as gx
          where g.gene = gx.gene and
                g.chr_name = ? and g.gene_chrom_start <= ? and
                g.gene_chrom_end >= ?"
    );
    eval {
        $sth->execute( $chr, $vc_end, $vc_start );
    };
    return [] if($@);
    my @genes;
    while( my $row = $sth->fetchrow_arrayref() ) {
        next unless $row->[8]=~/HUMACE/;
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
            'db'        => $row->[7],
            'type'      => $row->[8]
        };
    }
    return $self->{$cache_name} = \@genes;
    return \@genes
}

sub fetch_virtualRepeatFeatures_start_end {
    my ( $self, $chr, $vc_start, $vc_end, $type, $glob_bp ) =@_;
    my $cache_name = "_repeats_$type"."_cache_$chr"."_$vc_start"."_$vc_end";
    return $self->{$cache_name} if( $self->{$cache_name} );
	my $glob_bp ||= 0;
    my $_db_name = $self->{'_lite_db_name'};

    my $sth = $self->prepare(
        "select r.id, r.hid,  r.chr_name, r.repeat_chrom_start, r.repeat_chrom_end, r.repeat_chrom_strand
           from $_db_name.repeat as r
          where r.chr_name = ? and r.repeat_chrom_start <= ? and r.repeat_chrom_start >= ? and r.repeat_chrom_end >= ?".
		  	( (defined $type && $type ne '') ? " and r.type = '$type'" : '' )
    );

    eval {
        $sth->execute( $chr, $vc_end, $vc_start-1000000, $vc_start);
    };
    return [] if($@);

	my @repeats;
	my $old_end = -99999999999999999;
	while( my $row = $sth->fetchrow_arrayref() ) {
      	my $end = $row->[4];
## Glob results! 
        next if($end < $old_end );
    	if($end < $old_end + $glob_bp) {
			$repeats[-1]->{'chr_end'} = $end;
	  	}	else {
			push @repeats, {
				'id'        => $row->[0],
				'hid'       => $row->[1],
    	        'chr_name'  => $row->[2],
        	    'chr_start' => $row->[3],
            	'chr_end'   => $end,
        	    'start'     => $row->[3]-$vc_start+1,
            	'end'       => $end     -$vc_start+1,
	            'strand'    => $row->[5],
			};
		}
	  	$old_end = $end;
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

        SELECT   snp_chrom_start,strand,chrom_strand,
                 refsnpid,
                 tscid, hgbaseid,clone 
        FROM   	 $_db_name.snp
        WHERE  	 chr_name='$chr' 
        AND      snp_chrom_start>$vc_start
	AND      snp_chrom_start<$vc_end
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
        
        my ($snp_start, $strand,$chrom_strand,$snpuid,$tscid, $hgbaseid,$acc) = @{$arr};
            
  # globbing
        
        my $key=$snpuid.$acc;           # for purpose of filtering duplicates
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
               -original_strand => $strand,
               -score => 1,
               -source_tag => 'dbSNP',
              );
            
            my $link = new Bio::Annotation::DBLink;
            $link->database('dbSNP');
            $link->primary_id($snpuid);
           $link->optional_id($acc);
            #add dbXref to Variation
            $snp->add_DBLink($link);
	    if ($hgbaseid) {
	      my $link2 = new Bio::Annotation::DBLink;
	      $link2->database('HGBASE');
	      $link2->primary_id($hgbaseid);
	      $link2->optional_id($acc);
	      $snp->add_DBLink($link2);
	    }
	    if ($tscid) {
	      my $link3 = new Bio::Annotation::DBLink;
	      $link3->database('TSC-CSHL');
	      $link3->primary_id($tscid);
	      $link3->optional_id($acc);
	      #add dbXref to Variation
	      $snp->add_DBLink($link3);
	    }
            $cl=$acc;
            # set for compatibility to Virtual Contigs
            $snp->seqname($acc);
            #add SNP to the list
            push(@variations, $snp);
        }                               # if ! $seen{$key}
      }                                    # while a row from select statement

    #&eprof_end('snp-sql-object');
    
    return @variations;

}







1;
__END__

