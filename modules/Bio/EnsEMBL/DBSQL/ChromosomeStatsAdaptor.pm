package Bio::EnsEMBL::DBSQL::ChromosomeStatsAdaptor;

use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::ChromosomeStats;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);


sub get_stats_per_chromosome
{

    my ($self,$chromosome)=@_;

    $self->throw("I need a chromosome") unless defined $chromosome;

    my $query="select chromosome_id,id,name,species_id,known_genes,unknown_genes,snps,length 
               from chromosome where chromosome_id='$chromosome'";

    my $sth = $self->db->prepare($query);
    my $res = $sth->execute();

    my ($chr_id,$id,$name,$species,$known_genes,$unknown_genes,$snps,$length);

    $sth->bind_columns(undef,\$chr_id,\$id,\$name,\$species,\$known_genes,\$unknown_genes,\$snps,\$length);

    my @density_array;

    while( $sth->fetch() ) {
	my $cs=Bio::EnsEMBL::ChromosomeStats->new();
	$cs->chromosome_id($chr_id);
	$cs->id($id);
	$cs->name($name);
	$cs->species($species);
	$cs->known_genes($known_genes);
	$cs->unknown_genes($unknown_genes);
	$cs->snps($snps);
	$cs->length($length);

	push @density_array,$cs;
    }

    return @density_array;
}


