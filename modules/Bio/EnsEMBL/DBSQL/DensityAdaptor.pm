package Bio::EnsEMBL::DBSQL::DensityAdaptor;

use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::DensityWindow;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);


sub get_density_per_chromosome_type
{

    my ($self,$chromosome,$type)=@_;

    $self->throw("I need a chromosome") unless defined $chromosome;
    $self->throw("I need a type") unless defined $type;

    my $query="select chr_start,chr_end,value from map_density 
               where chromosome_id='$chromosome' and type='$type'";

    my $sth = $self->db->prepare($query);
    my $res = $sth->execute();

    my ($chr_start,$chr_end,$value);

    $sth->bind_columns(undef,\$chr_start,\$chr_end,\$value);

    my @density_array;

    while( $sth->fetch() ) {
	my $dw=Bio::EnsEMBL::DensityWindow->new();
	$dw->start($chr_start);
	$dw->end($chr_end);
	$dw->value($value);
	push @density_array,$dw;
    }

    return @density_array;
}


