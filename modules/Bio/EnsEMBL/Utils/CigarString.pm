=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=head1 AUTHOR

Juguang Xiao <juguang@tll.org.sg>

=cut

=head1 NAME

Bio::EnsEMBL::Utils::CigarString, a utilites module to generate cigar
strings

=head1 DESCRIPTION

Sequence alignment hits were previously stored within the core database
as ungapped alignments. This imposed 2 major constraints on alignments:

a) alignments for a single hit record would require multiple rows in the
   database, and
b) it was not possible to accurately retrieve the exact original
   alignment.

Therefore, in the new branch sequence alignments are now stored as
ungapped alignments in the cigar line format (where CIGAR stands for
Concise Idiosyncratic Gapped Alignment Report).

In the cigar line format alignments are stored as follows:

  M: Match
  D: Deletion
  I: Insertion

An example of an alignment for a hypthetical protein match is shown
below:


  Query:   42 PGPAGLP----GSVGLQGPRGLRGPLP-GPLGPPL...
              PG    P    G     GP   R      PLGP
  Sbjct: 1672 PGTP*TPLVPLGPWVPLGPSSPR--LPSGPLGPTD...

This would be stored in the protein_align_feature table as the following
cigar line:

  7M4D12M2I2MD7M

=cut

package Bio::EnsEMBL::Utils::CigarString;

use strict;
use vars qw(@ISA);

=head2 split_hsp

    Name  : split_hsp (this name is derived from the original sub in BlastWorn)
    Usage : my $hsp; # a ready Bio::Search::HSP::GenericHSP object.
my $factory = new Bio::EnsEMBL::Utils::CigarString;
my $cigar_string = $factory->split_hsp($hsp);
  
    Function: generate cigar string.
    Argument: a HSP object.
    Returns : a text string.
  
=cut

sub split_hsp {
    my ($self, $hsp) = @_;

    $self->throw("a defined object needed") unless($hsp && defined($hsp));
    unless(ref($hsp) && $hsp->isa('Bio::Search::HSP::GenericHSP')){
        $self->throw("a HSP object needed");
    }

    my ($qtype, $htype) = $self->_findTypes($hsp);
    my ($qstrand, $hstrand) = $self->_findStrands($hsp);
    my ($qinc, $hinc) = $self->_findIncrements($qstrand,$hstrand,$qtype,$htype);

    my @gaps = ();
    my @qchars = split(//, $hsp->query_string);
    my @hchars = split(//, $hsp->hit_string);
    my $qstart;
    if($qstrand == 1){
        $qstart = $hsp->query->start;
    }elsif($qstart == -1){
        $qstart = $hsp->query->end;
    }else{
        $self->warn("[$qstart], invalid strand value on query");
        $qstart = $hsp->query->start; 
        # Is this a SearchIO's bug???
    }
    
    my $hstart; 
    if($hstrand == 1){
        $hstart = $hsp->subject->start;
    }elsif($hstrand != -1){
        $hstart = $hsp->subject->end;
    }else{
        $self->throw("[$hstart], invalid strand value on subject");
    }

    my $qend = $qstart;
    my $hend = $hstart;
    my $count = 0;
    my $found = 0;

    my @align_coordinates = ();
    while($count <= $#qchars){
        if($qchars[$count] ne '-' && $hchars[$count] ne '-') {
            $qend += $qinc;
            $hend += $hinc;
            $found = 1;
        }else{ # gapped region
            push(@align_coordinates, [$qstart, $hstart]) if($found == 1);

            $qstart = $qend;
            $qstart += $qinc if($qchars[$count] ne '-');

            $hstart = $hend;
            $hstart += $hinc if($hchars[$count] ne '-');

            $qend = $qstart;
            $hend = $hstart;
            $found = 0;
        }
        $count++;
    }

    if($found){
        push(@align_coordinates, [$qstart, $hstart]);
    }

    my $cigar_string = "";
    my $last = $#align_coordinates;
    if($last >= 0){
        for(my $i=0; $i<$last; $i++){
            my $q_this_start = $align_coordinates[$i]->[0];
            my $q_next_start = $align_coordinates[$i+1]->[0];
            my $q_length = ($q_next_start-$q_this_start-1)*$qinc;
            $q_length = abs($q_length);
            my $h_this_start = $align_coordinates[$i]->[1];
            my $h_next_start = $align_coordinates[$i+1]->[1];
            my $h_length = ($h_next_start-$h_this_start-1)*$hinc;
            $h_length = abs($h_length);

            my $diff = $q_length - $h_length;
            if($diff > 0){ # Insertion
                $cigar_string .= $diff unless($diff == 1);
                $cigar_string .= 'I';
            }elsif($diff < 0){ # Deletion
                $cigar_string .= -$diff unless($diff == -1);
                $cigar_string .= 'D';
            }else{ # e.g $diff == 0, Match
                $cigar_string .= $q_length unless($q_length == 1);
                $cigar_string .= 'M';
            }
                
        } # for
    } # if
    
    return $cigar_string;
}


sub _findStrands {
    my ($self,$hsp) = @_;
    
    my $qstrand = $hsp->query->strand;
    unless($qstrand == 1 || $qstrand == -1){
        $self->warn("query's strand value is neither 1 or -1");
        $qstrand = 1;
    }
    
    my $hstrand = $hsp->subject->strand;
    unless($hstrand == 1 || $hstrand == -1){
        $self->warn("subject's strand value is neither 1 or -1");
        $hstrand = 1;
    }
    
    return ( $qstrand, $hstrand);
}

sub _findTypes {
    my ($self,$hsp) = @_;

    my $type1;
    my $type2;
    my $len1 = $hsp->query->length();
    my $len2 = $hsp->subject->length();

    if ($len1/$len2 > 2) {
        $type1 = 'dna';
        $type2 = 'pep';
    } elsif ($len2/$len1 > 2) {
        $type1 = 'pep';
        $type2 = 'dna';
    } else {
        $type1 = 'dna';
        $type2 = 'dna';
    }

    return ($type1,$type2);
}

sub _findIncrements {
    my ($self,$qstrand,$hstrand,$qtype,$htype) = @_;

    my $qinc   = 1 * $qstrand;
    my $hinc   = 1 * $hstrand;

    if ($qtype eq 'dna' && $htype eq 'pep') {
    $qinc = 3 * $qinc;
    }
    if ($qtype eq 'pep' && $htype eq 'dna') {
    $hinc = 3 * $hinc;
    }

    return ($qinc,$hinc);
}

# This is a core logic of cigar string. The finite state machine theory is 
# apply. See the below table, x-axis represents the input, with 3 options: 
# (+/+) -- Both current query and subject bases are non-gap. Match
# (-/+) -- The current query base is gap, but subject not. Deletion
# (+/-) -- The current subject base is gap, but query not. Insertion
# While the y-axis means the current state with letter 'M', 'D', 'I'
#
# The content of this table is the action taken in response of the input and 
# the current state.
# R     remain the state, counter increment.
# G;X   generate the cigar line based on the current state and counter;
#       clear the counter to zero and change to the state X
#       
#       ||  +/+  |  -/+ |  +/-  |
# -------+----------------------+
#    M  ||   R   |  G;D |   G;I |
# ------------------------------+
#    D  ||  G;M  |   R  |   G;I |
# ------------------------------+   
#    I  ||  G;M  |  G;D |    R  |
#

=head2 generate_cigar_string

  Name : generate_cigar_string
  Usage: $cigar_string = $self->generate_cigar_string(\@qchars, \@hchars);
  Function: generate the cigar string for a piece of alignment.
  Args:     2 array references. The lengths of 2 arrays are the same
  Return:   a cigar string

=cut

# Developer's Note: The method is originally abstracted from the concept of 
# cigar string. It only asks the essential information of 2 sequence characters
# of the alignment, while the BlastWorn::split_HSP asks more unused information
# for cigar string, which is useful to form align_coordinates. - Juguang

my ($count, $state); # strictly only used in the following 2 subs

sub generate_cigar_string {

#   my ($self, $qstart, $hstart, $qinc, $hinc, $qchars_ref, $hchars_ref) = @_;

    my ($self, $qchars_ref, $hchars_ref) = @_;
    
    my @qchars = @{$qchars_ref};
    my @hchars = @{$hchars_ref};

    unless(scalar(@qchars) == scalar(@hchars)){
        $self->throw("two sequences are not equal in lengths");
    }
    
    $count = 0;
    $state = 'M';    # the current state of gaps, (M, D, I) 
   
    my $cigar_string = '';
    for(my $i=0; $i <= $#qchars; $i++){
        my $qchar = $qchars[$i];
        my $hchar = $hchars[$i];
        if($qchar ne '-' && $hchar ne '-'){ # Match
            $cigar_string .= $self->_sub_cigar_string('M');
        }elsif($qchar eq '-'){ # Deletion
            $cigar_string .= $self->_sub_cigar_string('D');
        }elsif($hchar eq '-'){ # Insertion
            $cigar_string .= $self->_sub_cigar_string('I');
        }else{
            $self->throw("Impossible state that 2 gaps on each seq aligned");
        }
    }
    $cigar_string .= $self->_sub_cigar_string('X'); # not forget the tail.
    return $cigar_string;
}

sub _sub_cigar_string {
    my ($self, $new_state) = @_;
    my $sub_cigar_string = '';
    if($state eq $new_state){
        $count++; # Remain the state and increase the counter
    }else{
        $sub_cigar_string .= $count unless $count == 1;
        $sub_cigar_string .= $state;
        $count = 1;
        $state = $new_state;
    }
    return $sub_cigar_string;
}

=head2 generate_cigar_string_by_hsp
  
  Name :    generate_cigar_string_by_hsp
  Usage :   my $hsp; # a ready GenericHSP object
my $cigar_string = $self->generate_cigar_string_by_hsp($hsp);
  Function: generate a cigar string by given HSP object.
  Args :    a GenericHSP object
  Returns:  a text string of cigar string

=cut

sub generate_cigar_string_by_hsp {
    my ($self, $hsp) = @_;

    unless(ref($hsp) && $hsp->isa('Bio::Search::HSP::GenericHSP')){
        $self->throw("A GenericHSP object needed");
    }

    my @qchars = split(//, $hsp->query_string);
    my @hchars = split(//, $hsp->hit_string);

    return $self->generate_cigar_string(\@qchars, \@hchars);
}

1;
