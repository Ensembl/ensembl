#
# bioperl module for output from MSPcrunch
#
# Cared for by Michele Clamp  <michele@sanger.ac.uk>
#
# Copyright Michele Clamp
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Analysis::MSPcrunch - reads and stores output from MSPcrunch

=head1 SYNOPSIS

    my $msp    = new Bio::EnsEMBL::Analysis::MSPcrunch(-file => $myfile,
                                                       -type => $homol_type,
						       -source_tag => $source,
						       )

Extracting data

    my @homols = $msp->each_Homol;

Returns an array of Bio::EnsEMBL::FeaturePair;

=head1 DESCRIPTION

Object to parse and store output from MSPcrunch

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Analysis::MSPcrunch;

use vars qw(@ISA);
use strict;


# Object preamble - inherits from Bio::Root::Object;

use Bio::Root::Object;

use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Analysis::MSPType;
use Bio::EnsEMBL::FeaturePair;
#used to be Pep_SeqFeature, (removed by mcvicker)
use Bio::EnsEMBL::SeqFeature;

use FileHandle;
# Inherits from the base bioperl object
@ISA = qw(Bio::Root::Object);

# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;
  
  #my $make = $self->SUPER::_initialize;

  # Input variables
  # ---------------
  # Filename.
  # Type of msp file

  my ($contig_id,$mspfile,$type,$source_tag) = $self->_rearrange([qw(CONTIG_ID
								     FILE
								     TYPE
								     SOURCE_TAG
								     )],@args);
  # Stored data
  # -----------
  # These are the homols parsed from the mspfile

  $self->{_homols} = [];

  $self->source_tag($source_tag);
  $self->type      ($type);
  $self->mspfile   ($mspfile);
  $self->contig_id ($contig_id);
  $self->_parse;

  return $self; # success - we hope!
}

=head2 _parse

  Title   : _parse
  Usage   : $self->_parse
  Function: Internal method to parse a MSPcrunch file
  Returns : Nothing
  Args    : None

=cut

sub _parse {
    my ($self) = @_;

    my $file = $self->mspfile;
    
    $self->_make_analysis($file);

    my $fh = new FileHandle("<$file") || $self->throw("Can't open $file");

    while (defined(my $line = <$fh>)) {
	eval {
	    my $homol = $self->_read_Homol($line);
	    $self->add_Homol($homol);
	};
	if ($@) {
	    print STDERR "Couldn't parse line [$line] [$@]\n";
	}
    }

    $fh->close;

}

=head2 _read_Homol

  Title   : _read_Homol
  Usage   : $self->_read_Homol($line);
  Function: Converts a line from a MSPcrunch output file to a Bio::EnsEMBL::FeaturePair
  Returns : Nothing
  Args    : Bio::SeqFeature::String

=cut

sub _read_Homol {
    my ($self,$line) = @_;

    my ($score, $pid,
        $qstart, $qend, $id1,
        $hstart, $hend, $id2,
        $desc) = split(' ', $line, 9);

    my $strand1 = 1;
    my $strand2 = 1;

    # Make sure the start coords are less than the end coords.

    if ($qstart > $qend ) {
	my $tmp    = $qend;
  	   $qend   = $qstart;
	   $qstart = $tmp;

	$strand1 = -1;
    }
    if ($hstart > $hend ) {

	my $tmp    = $hend;
	$hend   = $hstart;
	$hstart = $tmp;

	$strand2 = -1;
    }


    print STDERR "Strands are $strand1 $strand2\n";

    my $strand = 1;

    if ($strand1 == -1 || $strand2 == -1) {
	$strand = -1;
    }

    my ($type1,$type2) = $self->get_types;


    my $sf1;
    my $sf2;

#(mcvicker) removed Pep_SeqFeature garbage
#    if ($type1 eq "PEP") {

#	$sf1 = new Bio::EnsEMBL::Pep_SeqFeature(-start  => $qstart,
#						-end    => $qend,
#						-strand => 1
#						);
#	#$sf1->start_frac(1);
#	#$sf1->end_frac  (3);

#    } else {
	$sf1 = new Bio::EnsEMBL::SeqFeature(-start  => $qstart,
					    -end    => $qend,
					    -strand => 1);
#    }
    
#    if ($type2 eq "PEP") {
#	$sf2 = new Bio::EnsEMBL::Pep_SeqFeature(-start  => $hstart,
#						-end    => $hend,
#						-strand => $strand);
#	#$sf2->start_frac(1);
#	#$sf2->end_frac  (3);

#    } else {
	$sf2 = new Bio::EnsEMBL::SeqFeature(-start  => $hstart,
					    -end    => $hend,
					    -strand => $strand);
#    }

    $sf1->score($score);
    $sf2->score($score);

    $sf1->seqname($self->contig_id);
    $sf2->seqname($id2);

    $sf1->primary_tag('similarity');
    $sf2->primary_tag('similarity');

    if ($self->source_tag eq "") {
	print(STDERR "ERROR: No source tag in mspfile " . $self->mspfile . "\n");
    }
    $sf1->source_tag($self->source_tag);
    $sf2->source_tag($self->source_tag);


    my $fp = new Bio::EnsEMBL::FeaturePair(-feature1 => $sf1,
					   -feature2 => $sf2,
					   );
	
    my $anal = $self->analysis;

    $self->throw("Can't make analysis\n") unless defined($self->analysis);
    $sf1->analysis($self->analysis);
    $sf2->analysis($self->analysis);
    $fp->analysis ($self->analysis);
    return ($fp);
}
	
=head2 each_Homol

  Title   : each_Homol
  Usage   : my @homols = $self->each_Homol
  Function: Returns all the MSPcrunch output lines as homol features
  Returns : array of Bio::SeqFeature::Homol
  Args    : none

=cut

sub each_Homol {
    my ($self) = @_;

    $self->throw("No homol array defined in Bio::EnsEMBL::Analysis::MSPcrunch") unless defined($self->{_homols});
    
    return @{$self->{_homols}};
}

=head2 add_Homol

  Title   : add_Homol
  Usage   : $self->add_Homol($homol);
  Function: Adds a homol feature to the object.
  Returns : Nothing
  Args    : Bio::SeqFeature::Homol

=cut

sub add_Homol {
    my ($self,$homol) = @_;
    
    
    $self->throw("Argument to Bio::EnsEMBL::Analysis::MSPcrunch->add_Homol is not a Bio::EnsEMBL::FeaturePair") unless $homol->isa("Bio::EnsEMBL::FeaturePair");
    
    $homol->validate;

    push(@{$self->{_homols}},$homol);

}


=head2 mspfile

  Title   : mspfile
  Usage   : $self->mspfile
  Function: Get/set method for the MSPcrunch filename
  Returns : String
  Args    : String

=cut

sub mspfile {
    my ($self,$file) = @_;
    
    if (defined($file)) {
	$self->throw("MSPcrunch output file $file doesn't exist") unless -e $file;
	
	$self->{_mspfile} = $file;
    }
    
    return $self->{_mspfile};
}

=head2 contig_id

  Title   : contig_id
  Usage   : $self->contig_id
  Function: Get/set method for the contig id
  Returns : String
  Args    : String

=cut

sub contig_id {
    my ($self,$id) = @_;
    
    if (defined($id)) {
	$self->{_contig_id} = $id;
    }
    
    return $self->{_contig_id};
}


=head2 source_tag

  Title   : source_tag
  Usage   : $self->source_tag
  Function: Get/set method for the source_tag
  Returns : String
  Args    : String

=cut

sub source_tag {
    my ($self,$arg) = @_;


    if (defined($arg)) {
	$self->{_source_tag} = $arg;
    }

    return $self->{_source_tag};
}



=head2 type

  Title   : type
  Usage   : my $type = $msp->type;
  Function: Get/set method for the type of MSP output we have
  Returns : String
  Args    : none

=cut


BEGIN {

    my %allowed_type = map {$_, 1} qw(DNA-DNA DNA-PEP PEP-DNA PEP-PEP);

    sub type {
        my ($self, $type) = @_;

        if (defined($type)) {

	    unless ($allowed_type{$type}) {
	        $self->throw("Wrong type '$type' for MSPcrunch entered : allowed values are "
                             . join(' ', sort keys %allowed_type));
	    }
	    $self->{_type} = $type;
        }
        
        #use Carp;
        #use Data::Dumper;
        #confess("No type", Dumper($self)) unless defined $self->{_type};
        
        return $self->{_type};
    }
}

=head2 get_types

  Title   : get_types
  Usage   : my ($type1,$type2) = $msp->get_types
  Function: Fetches the sequence type of the query and target sequence in the MSP hit
  Returns : String
  Args    : none

=cut

sub get_types {
    my ($self) = @_;

    my $typestring = $self->type;

    my ($type1,$type2) = split(/-/,$typestring);

    return ($type1,$type2);
}


sub analysis {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->throw("Argument is not Bio::EnsEMBL::Analysis::Analysis object") 
	    unless $arg->isa("Bio::EnsEMBL::Analysis");

	$self->{_analysis} = $arg;
    }

    return $self->{_analysis};
}

sub _make_analysis {
    my ($self,$mspfile) = @_;

    (my $ext  = $mspfile) =~ s/.*(\..*.msptmp)/$1/;

    my $test    = Bio::EnsEMBL::Analysis::MSPType->each_MSPType;
    my $MSPType = Bio::EnsEMBL::Analysis::MSPType->extension2MSPType($ext);
    my $anal    = new Bio::EnsEMBL::Analysis;
    
    $anal->db             ($MSPType->[2]);
    $anal->db_version     ($MSPType->[6]);
    $anal->program        ($MSPType->[1]);
    $anal->program_version($MSPType->[7]);
    $anal->gff_source     ($MSPType->[1]);
    $anal->gff_feature    ($MSPType->[8]);
    
    $self->analysis($anal);
}

1;
