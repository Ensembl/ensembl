#
# DensityAdaptor module
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

DensityAdaptor

=head1 SYNOPSIS

my $obj= Bio::EnsEMBL::DBSQL::Obj->new(-dbname=>'ens08',-user=>'ensadmin',-host=>'ensrv5');

my $da= Bio::EnsEMBL::DBSQL::DensityAdaptor->new($obj);

my $binvalueset = $da->get_density_per_chromosome_type('1','gene');

=head1 DESCRIPTION

Database adaptor for the density plot objects.

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Lite::DensityPlot::DensityAdaptor;

use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::DensityPlot::BinValue;
use Bio::EnsEMBL::DensityPlot::BinValueSet;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);


=head2 get_density_per_chromosome_type

 Title   : get_density_per_chromosome_type
 Usage   : my $BinValueSet = $DensAdapt->get_density_per_chromosome_type(1,"gene");
 Function:
 Example :
 Returns : A DensityPlot::BinValueSet containing BinValue objects loaded with the data from the db for the specified chromosome and type. 
 Args    : A chromosome id (1-24) and data type (snp, gene, kngene, etc)

=cut

sub get_density_per_chromosome_type
{

    my ($self,$chromosome,$type)=@_;

    $self->throw("I need a chromosome") unless defined $chromosome;
    $self->throw("I need a type") unless defined $type;

    my $chr_adaptor = $self->db->get_ChromosomeAdaptor();
    my $chr_id = $chr_adaptor->fetch_by_chr_name($chromosome)->dbID();

    my $query="	SELECT	chr_start,
			chr_end,
			value 
		FROM	map_density 
		WHERE	chromosome_id = $chr_id 
		AND	type = '$type' 
		ORDER BY chr_start";

    my $sth = $self->db->prepare($query);
    my $res = $sth->execute();

    my ($chr_start,$chr_end,$value);

    $sth->bind_columns(undef,\$chr_start,\$chr_end,\$value);

    my $valueset = new Bio::EnsEMBL::DensityPlot::BinValueSet;
    while( $sth->fetch() ) {
    	my $binvalues = new Bio::EnsEMBL::DensityPlot::BinValue;
    	$binvalues->chromosomestart($chr_start);
    	$binvalues->chromosomeend($chr_end);
    	$binvalues->value($value);
    	$valueset->add_binvalue($binvalues);
    }

    return $valueset;
}

1;
