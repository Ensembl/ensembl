# EnsEMBL Gene reading writing adaptor for mySQL
#
# Copyright EMBL-EBI 2001
#
# Author: Arne Stabenau
# based on 
# Elia Stupkas Gene_Obj
# 
# Date : 20.02.2001
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
                g.chr_name = ? and g.gene_chrom_start <= ? and
                g.gene_chrom_end >= ?"
    );
    eval {
        $sth->execute( $chr, $vc_end, $vc_start );
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
                g.chr_name = ? and g.gene_chrom_start <= ? and
                g.gene_chrom_end >= ?"
    );
    eval {
        $sth->execute( $chr, $vc_end, $vc_start );
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

sub fetch_virtualRepeatFeatures_start_end {
    my ( $self, $chr, $vc_start, $vc_end, $type, $glob_bp ) =@_;
    my $cache_name = "_repeats_$type"."_cache_$chr"."_$vc_start"."_$vc_end";
    return $self->{$cache_name} if( $self->{$cache_name} );
	my $glob_bp ||= 0;
    my $_db_name = $self->{'_lite_db_name'};
    my $sth = $self->prepare(
        "select r.id, r.hid,  r.chr_name, r.repeat_chrom_start, r.repeat_chrom_end, r.repeat_chrom_strand
           from $_db_name.repeat as r
          where r.chr_name = ? and r.repeat_chrom_start <= ? and r.repeat_chrom_end >= ?".
		  	( (defined $type && $type ne '') ? " and r.type = '$type'" : '' )
    );

    eval {
        $sth->execute( $chr, $vc_end, $vc_start);
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
            	'end'       => $end-$vc_start+1,
	            'strand'    => $row->[5],
			};
		}
	  	$old_end = $end;
    }
    return $self->{$cache_name} = \@repeats;
    return \@repeats;
}

1;
__END__

