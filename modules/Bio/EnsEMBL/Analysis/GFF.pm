#
# Object for parsing/storing GFF files
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

Bio::EnsEMBL::Analysis::GFF - Parses and stores data from GFF files

=head1 SYNOPSIS

    my $obj    = new Bio::EnsEMBL::Analysis::GFF(-file => $file)


=head1 DESCRIPTION

Object to store details of an analysis run

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Analysis::GFF;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object;

use Bio::Root::Object;
use Bio::SeqFeature::Generic;
use Bio::EnsEMBL::Repeat;
use FileHandle;

use Bio::EnsEMBL::FeaturePair;

# Inherits from the base bioperl object
@ISA = qw(Bio::Root::Object);

# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;
  
  #my $make = $self->SUPER::_initialize;

  my ($file,$type) = $self->_rearrange([qw(FILE
					   TYPE)],@args);

  # The array to store all the feature objects
  $self->{_features} = [];
  $self->type($type);
  $self->GFFFile($file);

  # Parse the GFF File
  $self->_parse();

  return $self; # success - we hope!
}


=head2 each_Feature

  Title   : each_Feature
  Usage   : @homols = $self->each_Feature
  Function: Returns each GFF feature stored in the object
  Returns : @Bio::SeqFeature::Generic
  Args    : none

=cut

sub each_Feature {
    my ($self) = @_;

    my @f = @{$self->{_features}};

    return @{$self->{_features}};
}

=head2 add_Feature {

  Title   : add_Feature
  Usage   : $self->add_Feature
  Function: Adds feature to the object
  Returns : Nothing
  Args    : Bio::SeqFeature::Generic

=cut

sub add_Feature {
    my ($self,$arg) = @_;
    
    my @f = @{$self->{_features}};

    if (defined($arg) && $arg->isa("Bio::EnsEMBL::SeqFeatureI")) {
	$arg->validate;
	push(@{$self->{_features}},$arg);
    } else {
	$self->throw("Feature : $arg : is not a Bio::EnsEMBL::SeqFeatureI");
    }

}

=head2 type

  Title   : type
  Usage   : $self->type($type);
  Function: Get/set method for the type of feature object to create
  Returns : Bio::SeqFeature::Generic
  Args    : string

=cut

sub type {
    my ($self,$arg) = @_;


    if (defined($arg)) {
	$self->{_type} = $arg;
    }

    return $self->{_type};
}

=head2 _parse

  Title   : _parse
  Usage   : $self->_parse
  Function: Parses the input GFF file and stores the features
            Is called automatically when a new filename is set in
            the object
  Returns : nothing
  Args    : none

=cut

sub _parse {
    my ($self) = @_;

    $self->throw("No GFF file input") unless $self->GFFFile;

    my $file = $self->GFFFile;
    
    $self->_make_analysis;
    
    open(IN,"<$file") || $self->throw("Can't open $file");
    while (defined( my $line = <IN>)) {
	if ($line !~ /^\#/) {
	    my $feature = $self->_parse_line($line);
	    $self->add_Feature($feature);
	}
    }
    
    close(IN);
}



=head2 _parse_line

  Title   : _parse_line
  Usage   : $self->_parse_line
  Function: Parses one line from a GFF file and returns a feature
  Returns : Bio::SeqFeature::Generic
  Args    : String

=cut

sub _parse_line {
    my ($self,$line) = @_;

    my ($seqname,$source,$feature,$start,$end,$score,$strand,$frame,$attrib) = split(/\t/,$line);

    my @attrib = split(' ',$attrib);

    if ($start > $end) {
	my $tmp = $start;
	$start = $end;
	$end   = $tmp;
    }

    if ($strand eq "+") { $strand = 1;}
    if ($strand eq "-") { $strand = -1;}
    
    if ($frame eq ".") {  $frame = "";}

    my $f = new Bio::EnsEMBL::SeqFeature(-seqname     => $seqname,
					 -source_tag  => $source,
					 -primary_tag => $feature,
					 -start       => $start,
					 -end         => $end,
					 -strand      => $strand,
					 -frame       => $frame,
					 -score       => $score,
					 );

    $f->seqname ($seqname);
    $f->analysis($self->analysis);

    if ($attrib[0] eq "Target" && $feature eq "similarity") {
	$f = $self->_parse_attrib($f,$feature,@attrib);
    }
    
    return $f;
}

    
=head2 _parse_attrib

  Title   : _parse_attrib
  Usage   : my $f = $self->parse_attrib($feature,$feature_tag,@attrib)
  Function: Get/set method for the GFFfile file name
  Returns : String
  Args    : String

=cut

sub _parse_attrib {
    my ($self,$feature,$feature_tag,@attrib) = @_;
    
    # I am not overly sure whether this should be done here
    if ($feature_tag eq "similarity" && $attrib[0] eq "Target") {
	my $hname = $attrib[1];
  	   $hname =~ s/\"//g;

	my $hstart = $attrib[2];
	my $hend   = $attrib[3];

	if ($hstart < 1 || $hend < $hstart) { 
	    $self->warn("Invalid homol coordinates $hstart - $hend. Skipping feature\n");
	    return $feature;
	}

	my $homol = new Bio::EnsEMBL::SeqFeature(-start       => $feature->start,
						 -end         => $feature->end,
						 -strand      => $feature->strand,
						 -frame       => $feature->frame,
						 -primary_tag => $feature->primary_tag,
						 -source_tag  => $feature->source_tag,
						 -score       => $feature->score,
						 -seqname     => $feature->seqname,
						 );
	
	$homol->seqname($feature->seqname);
	
	my $newf = new Bio::EnsEMBL::SeqFeature  (-start       => $hstart,
						  -end         => $hend,
						  -strand      => $feature->strand,
						  -frame       => $feature->frame,
						  -seqname     => $hname,
						  -source_tag  => $feature->source_tag,
						  -score       => $feature->score,
						  -primary_tag => $feature->primary_tag,
						  );
	$newf ->seqname($hname);
	$newf ->analysis($self->analysis);
	$homol->analysis($self->analysis);

	my $fp;

	if ($self->type eq "Repeat") {
	    $fp = new Bio::EnsEMBL::Repeat(-feature1 => $homol,
					   -feature2 => $newf,
					   );

	} else {
	    $fp = new Bio::EnsEMBL::FeaturePair(-feature1 => $homol,
						-feature2 => $newf,
						);
	}
	$fp->analysis($self->analysis);
	return $fp;
    }  else {
	$self->throw("Can't make homol object from attrib tags");
    }
}

=head2 GFFFile

  Title   : GFFFile
  Usage   : $self->GFFfile($file)
  Function: Get/set method for the GFFfile file name
  Returns : String
  Args    : String

=cut

sub GFFFile {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_GFFFile} = $arg;
    }

    return $self->{_GFFFile};
}

sub _make_analysis {
    my ($self) = @_;
    
    (my $ext  = $self->GFFFile) =~ s/.*(\..*.out.gff)/$1/;
    
    my $test    = Bio::EnsEMBL::Analysis::MSPType->each_MSPType;
    my $MSPType = Bio::EnsEMBL::Analysis::MSPType->extension2MSPType($ext);
    my $anal    = new Bio::EnsEMBL::Analysis::Analysis;
    
    $anal->db             ($MSPType->[2]);
    $anal->db_version     ($MSPType->[6]);
    $anal->program        ($MSPType->[1]);
    $anal->program_version($MSPType->[7]);
    $anal->gff_source     ($MSPType->[1]);
    $anal->gff_feature    ($MSPType->[8]);
    
    $self->analysis($anal);
}

sub analysis {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_analysis} = $arg;
    } 

    return $self->{_analysis};
}
1;
