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


=head2 fetch_by_chromosome_start_end

 Title   : fetch_by_chromosome_start_end
 Usage   : $band_obj = $kary_adp->fetch_by_chromosome_position('chr1',10000, 2000);
 Function: Retrieves KaryotypeBand objects by chromosome ('chrN' notation)
           and its absolute "golden-path" start/end on that chromosome.
 Example :
 Returns : A KaryotypeBand object list
 Args    : Chromosome id (chrN notation) and start, end in absolute basepairs

=cut

sub fetch_by_chromosome_start_end{
    my ($self,$chr,$start,$end) = @_;

    $self->throw("Need both chromosome and start/end") unless (defined $start && defined $end);

    my $sth = $self->prepare("	SELECT	chr_start,
					chr_end,
					band,
					stain
				FROM	karyotype 
				WHERE	chr_name = '$chr'
				AND	$start <= chr_end 
				AND	$end > chr_start 
			     ");

    $sth->execute;
	my @bands = ();
	my ($chr_start,$chr_end,$band,$stain) = ();
		
    while (($chr_start,$chr_end,$band,$stain) = $sth->fetchrow_array()){
    	last unless defined $band;
    	my $band_obj = Bio::EnsEMBL::KaryotypeBand->new();
    	$band_obj->name($band);
    	$band_obj->chromosome($chr);
    	$band_obj->start($chr_start);
    	$band_obj->end($chr_end);
    	$band_obj->stain($stain);
		#print STDERR "Kary Get: $chr_start,$chr_end,$band,$stain\n";
		push (@bands, $band_obj);
	}

    return @bands;
}


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
	
	my ($chr_start,$chr_end,$band,$stain)  = $sth->fetchrow_array();
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

    my $sth = $self->prepare(
        "select	chr_start, chr_end, stain
           from karyotype 
          where chr_name = ? and band = ?"
    );
    $sth->execute( $chr, $name );

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

=head2 fetch_by_chromosome_name_virtual

 Title   : fetch_by_chromosome_name_virtual
 Usage   : $band_obj = $kary_adp->fetch_by_chromosome_name('chr1','q23.1');
 Function: Retrieves a KaryotypeBand obj by its chromosome ('chrN' notation)
           and its band name (e.g. 'q23.1'). { can also use a non-fully declared band:
	   e.g. q23, q23.3
 Example :
 Returns : A KaryotypeBand object BUT DOES NOT RETURN STAIN INFORMATION!!!!!
         : (as we can not get this from the information we retrieve)
 Args    : Chromosome id (chrN notation) and (partial?) band name 

=cut

sub fetch_by_chromosome_name_virtual {
    my ($self,$chr,$name) = @_;

    $self->throw("Need both chromosome and name") unless defined $name;

    my $sth;
    $sth = $self->prepare(
        "select	chr_start, chr_end
           from karyotype 
          where chr_name = ? and band = ?"
    );
    $sth->execute( $chr, $name );
    my ( $chr_start, $chr_end ) = $sth->fetchrow_array;

	unless($chr_start) {
	    $sth = $self->prepare(
            "select min(chr_start), max(chr_end)
               from karyotype 
              where chr_name = '$chr' and band like '$name\%'"
        );
	    $sth->execute;
	    ( $chr_start, $chr_end ) = $sth->fetchrow_array;
	}

    return undef unless defined $chr_start;

    my $band_obj = Bio::EnsEMBL::KaryotypeBand->new();
    $band_obj->name($name);
    $band_obj->chromosome($chr);
    $band_obj->start($chr_start);
    $band_obj->end($chr_end);
    $band_obj->stain('white');

    return $band_obj;
}
=head2 fetch_by_name

 Title   : fetch_by_name
 Usage   : $band_obj = $kary_adp->fetch_by_chromosome_name('q23.1');
 Function: Retrieves a KaryotypeBand obj by its band name (e.g. 'q23.1').
 Example :
 Returns : A KaryotypeBand object
 Args    : Band name 

=cut


sub fetch_by_name{
    my ($self,$chr_name,$name) = @_;

    $self->throw("Need both name") unless defined $name;

    my $sth = $self->prepare(
        "select	chr_name, chr_start, chr_end, stain
           from karyotype 
          where band = ?"
    );
    $sth->execute( $name );

    my ($chr, $chr_start,$chr_end,$stain) = $sth->fetchrow_array;

    return undef unless defined $chr_start;

    my $band_obj = Bio::EnsEMBL::KaryotypeBand->new();
    $band_obj->name($name);
    $band_obj->chromosome($chr);
    $band_obj->start($chr_start);
    $band_obj->end($chr_end);
    $band_obj->stain($stain);

    return $band_obj;
}

=head2 fetch_by_name_virtual

 Title   : fetch_by_name_virtual
 Usage   : $band_obj = $kary_adp->fetch_by_name('q23.1');
 Function: Retrieves a KaryotypeBand obj by its 
           and its band name (e.g. 'q23.1'). { can also use a non-fully declared band:
	   e.g. q23, q23.3 }
 Example :
 Returns : A KaryotypeBand object BUT DOES NOT RETURN STAIN INFORMATION!!!!!
         : (as we can not get this from the information we retrieve)
 Args    : (partial?) band name 

=cut

sub fetch_by_name_virtual {
    my ($self,$chr_name,$name) = @_;

    $self->throw("Need both chromosome and name") unless defined $name;

    my $sth;
    $sth = $self->prepare(
        "select	chr_name, chr_name, chr_start, chr_end
           from karyotype 
          where band = ?"
    );
    $sth->execute( $name );
    my ( $chr, $chr_start, $chr_end ) = $sth->fetchrow_array;

	unless($chr_start) {
	    $sth = $self->prepare(
            "select chr_name, min(chr_start), max(chr_end)
               from karyotype 
              where band like '$name\%'"
        );
	    $sth->execute;
	    ( $chr, $chr_start, $chr_end ) = $sth->fetchrow_array;
	}

    return undef unless defined $chr_start;

    my $band_obj = Bio::EnsEMBL::KaryotypeBand->new();
    $band_obj->name($name);
    $band_obj->chromosome($chr);
    $band_obj->start($chr_start);
    $band_obj->end($chr_end);
    $band_obj->stain('white');

    return $band_obj;
}

=head2 fetch_chromosome_length

 Title   : fetch_chromosome_length
 Usage   : 
 Function: Returns length of a chromosome ('chrN' notation)
 Example :
 Returns : SV
 Args    : Chromosome id (chrN notation) 

=cut

sub fetch_chromosome_length {
    my ($self,$chr) = @_;

    $self->throw("Need a chromosome") unless defined $chr;

	# return a cached copy of the chromosome bands
	if (exists $self->{"_karyotype_band_cache_$chr"}){
		my @tmp = @{$self->{"_karyotype_band_cache_$chr"}};
		return $tmp[-1]->end();
	}

    my $sth = $self->prepare("	SELECT
											max(chr_end)
								FROM		karyotype 
								WHERE		chr_name = '$chr' 
			     			");

    $sth->execute;
    my ($chr_end) = $sth->fetchrow_array();
	return($chr_end);

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

    # return a cached copy of the chromosome bands
    if (exists $self->{"_karyotype_band_cache_$chr"}){
 	return(@{$self->{"_karyotype_band_cache_$chr"}});
    }

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

	# save a copy in the local cache
	$self->{"_karyotype_band_cache_$chr"} = \@bands;

    return @bands;
}

1;
