#
# Ensembl module for Bio::EnsEMBL::DBSQL::Assembly
#
# Cared for by James Smith <js5@sanger.ac.uk>
#
# Copyright James Smith
#
# You may distribute this module under the same terms as perl itself

=head1 NAME

Bio::EnsEMBL::Assembly

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Assembly objects encapsulate data pertaining to a single karyotype band.
Access these objects through a Bio::EnsEMBL::DBSQL::AssemblyAdaptor

=head1 AUTHOR

James Smith

This modules is part of the Ensembl project http://www.ensembl.org

=head1 CONTACT

Email js5@sanger.ac.uk

=cut

package Bio::EnsEMBL::Assembly;

use strict;
use vars qw(@ISA $AUTOLOAD);
use Bio::Root::RootI;
@ISA = qw(Bio::Root::RootI);


sub new {
    my ($class) = shift;
    my $self = {};
    my %data = @_;

    if($data{'type'} eq 'clone') {
        $self = {
            'contig_id'	    => $data{'contig_id'},
	    'contig_name'   => $data{'contig_name'},
	    'id'	    => $data{'clone_id'},
	    'name'	    => $data{'clone_name'},
	    'golden_start'  => $data{'golden_start'},
	    'golden_end'    => $data{'golden_start'} + $data{'length'} - 1,
	    'start'         => $data{'golden_start'} - $data{'vc_start'},
	    'end'           => $data{'golden_start'} - $data{'vc_start'} + $data{'length'} - 1,
	    'length'        => $data{'length'},
	    'seq_length'    => $data{'seq_length'},
	    'seq_start'     => $data{'golden_start'} - $data{'vc_start'} + int( ($data{'length'} - $data{'seq_length'})/2 ), 
	    'seq_end'       => $data{'golden_start'} - $data{'vc_start'} + int( ($data{'length'} - $data{'seq_length'})/2 ) + $data{'seq_length'} -1,
	    'state'         => $data{'state'},
	    'centre'        => $data{'centre'},
	    'chr'           => $data{'chromomsome'},
	    'sequenced'     => $data{'sequenced'},
	    'bac_f'         => $data{'bac_f'},
	    'bac_r'         => $data{'bac_r'},
        };
    } else {
        $self = {
	    'id'	    => $data{'contig_id'},
	    'name'	    => $data{'contig_name'},
	    'golden_start'  => $data{'golden_start'},
	    'golden_end'    => $data{'golden_start'} + $data{'length'} - 1,
	    'start'         => $data{'golden_start'} - $data{'vc_start'},
	    'end'           => $data{'golden_start'} - $data{'vc_start'} + $data{'length'} - 1,
	    'length'        => $data{'length'},
	    'chr'           => $data{'chromomsome'}
        };
    }
    bless $self,$class;

    return $self;
}

sub DESTROY { return 1; }

sub AUTOLOAD {
    my $self = shift;
#    no strict 'refs';
    my $var = $AUTOLOAD;
    $var =~ s/.*:://;		# remove class name if included...
    return $self->{$var} if (defined $self->{$var});

    print STDERR "Assembly warning ($self)- \$self->{'$var'} is undefined\n" unless ( $var eq  'bac_f' || $var eq 'bac_r');
    return undef;
}

1;
