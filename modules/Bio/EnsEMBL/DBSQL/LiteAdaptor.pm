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
    $_db_name = $self->{'_lite_db_name'};
    print STDERR "DB NAME: $_db_name\n\n";
    my $sth = $self->prepare(
        "select g.gene, g.name, 
                g.chr_name, g.gene_chrom_start, g.gene_chrom_end,
                g.chrom_strand, gx.display_id, gx.db_name
           from $_db_name.gene as g, $_db_name.gene_xref as gx
          where g.gene = gx.gene and
                g.chr_name = ? and g.gene_chrom_start <= ? and
                g.gene_chrom_end >= ?"
    );
    $sth->execute( $chr, $vc_end, $vc_start );
    my @genes;
    while( my $row = $sth->fetchrow_arrayref() ) {
        push @genes, {
            'gene'      => $row->[0],
            'stable_id' => $row->[1],
            'chr_name'  => $row->[2],
            'chr_start' => $row->[3]-$vc_start,
            'chr_end'   => $row->[4]-$vc_start,
            'strand'    => $row->[5],
            'synonym'   => $row->[6],
            'db'        => $row->[7]
        };
    }
    return \@genes
}
                

1;
__END__

