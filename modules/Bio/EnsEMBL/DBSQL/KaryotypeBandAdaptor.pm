#
# Ensembl module for Bio::EnsEMBL::DBSQL::KaryotypeBandAdaptor
#
# Cared for by James Stalker <jws@sanger.ac.uk>
#
# Copyright James Stalker
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::KaryotypeBandAdaptor

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Database adaptor to provide access to KaryotypeBand objects

=head1 AUTHOR

James Stalker

This modules is part of the Ensembl project http://www.ensembl.org

=head1 CONTACT

Email jws@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DBSQL::KaryotypeBandAdaptor;
use vars qw(@ISA);
use Bio::EnsEMBL::KaryotypeBand;
use strict;

# Object preamble - inherits from Bio::Root::RootI

use Bio::EnsEMBL::DBSQL::BaseAdaptor;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

# inherit new from BaseAdaptor


=head2 fetch_by_chromosome_position

 Title   : fetch_by_chromosome_position
 Usage   : $band_obj = $kary_adp->fetch_by_chromosome_position('chr1',10000);
 Function: Retrieves a KaryotypeBand obj by its chromosome ('chrN' notation)
           and its absolute "golden-path" position on that chromosome.
 Example :
 Returns : A KaryotypeBand object
 Args    : Chromosome id (chrN notation) and position in absolute basepairs

=cut

sub fetch_by_chromosome_position{
    my ($self,$chr,$position) = @_;

    $self->throw("Need both chromosome and position") unless defined $position;

    my $sth = $self->prepare("	SELECT	chr_start,
					chr_end,
					band,
					stain
				FROM	karyotype 
				WHERE	chr_name = '$chr'
				AND	$position <= chr_end 
				AND	$position > chr_start 
			     ");

    $sth->execute;
    my ($chr_start,$chr_end,$band,$stain) = $sth->fetchrow_array;

    return undef unless defined $band;

    my $band_obj = Bio::EnsEMBL::KaryotypeBand->new();
    $band_obj->name($band);
    $band_obj->chromosome($chr);
    $band_obj->start($chr_start);
    $band_obj->end($chr_end);
    $band_obj->stain($stain);

    return $band_obj;
}




=head2 fetch_by_chromosome_name

 Title   : fetch_by_chromosome_name
 Usage   : $band_obj = $kary_adp->fetch_by_chromosome_name('chr1','q23.1');
 Function: Retrieves a KaryotypeBand obj by its chromosome ('chrN' notation)
           and its band name (e.g. 'q23.1').
 Example :
 Returns : A KaryotypeBand object
 Args    : Chromosome id (chrN notation) and band name 

=cut

sub fetch_by_chromosome_name{
    my ($self,$chr,$name) = @_;

    $self->throw("Need both chromosome and name") unless defined $name;

    my $sth = $self->prepare("	SELECT	chr_start,
					chr_end,
					stain
				FROM	karyotype 
				WHERE	chr_name = '$chr' 
				AND	band = '$name'
			     ");

    $sth->execute;
    my ($chr_start,$chr_end,$stain) = $sth->fetchrow_array;

    return undef unless defined $chr_start;

    my $band_obj = Bio::EnsEMBL::KaryotypeBand->new();
    $band_obj->name($name);
    $band_obj->chromosome($chr);
    $band_obj->start($chr_start);
    $band_obj->end($chr_end);
    $band_obj->stain($stain);

    return $band_obj;
}




=head2 fetch_all_by_chromosome

 Title   : fetch_all_by_chromosome
 Usage   : @band_obj = $kary_adp->fetch_all_by_chromosome('chr1');
 Function: Retrieves all KaryotypeBand objects on a chromosome ('chrN' notation)
 Example :
 Returns : An array of KaryotypeBand objects
 Args    : Chromosome id (chrN notation) 

=cut

sub fetch_all_by_chromosome{
    my ($self,$chr) = @_;

    $self->throw("Need a chromosome") unless defined $chr;

    my $sth = $self->prepare("	SELECT	chr_start,
					chr_end,
					band,
					stain
				FROM	karyotype 
				WHERE	chr_name = '$chr' 
				ORDER BY chr_start
			     ");

    $sth->execute;
    my ($chr_start,$chr_end,$band,$stain);
    $sth->bind_columns(undef,\$chr_start,\$chr_end,\$band,\$stain);
    
    my @bands;

    while ($sth->fetch()){
	my $band_obj = Bio::EnsEMBL::KaryotypeBand->new();
	$band_obj->name($band);
	$band_obj->chromosome($chr);
	$band_obj->start($chr_start);
	$band_obj->end($chr_end);
	$band_obj->stain($stain);

	push @bands,$band_obj;
    }

    return @bands;
}

1;
