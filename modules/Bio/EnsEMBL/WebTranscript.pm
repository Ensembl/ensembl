
#
# BioPerl module for Bio::EnsEMBL::WebTranscript
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::WebTranscript - Specialisation of Transcript for specific web calls

=head1 SYNOPSIS

   # not appropiate

=head1 DESCRIPTION

This object is a sneaky specialisation specifically for web displays,
via get_all_Genes_exononly.  We need to remember whether things are on
or off this context...

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioperl-l@bio.perl.org

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Ewan Birney

Email birney@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::WebTranscript;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::EnsEMBL::Root

use Bio::EnsEMBL::Transcript;


@ISA = qw(Bio::EnsEMBL::Transcript);


# new is inherieted from Transcript


=head2 is_start_exon_in_context

 Title   : is_start_exon_in_context
 Usage   : $obj->is_start_exon_in_context($newval)
 Function: 
 Example : 
 Returns : value of is_start_exon_in_context
 Args    : newvalue (optional)


=cut

sub is_start_exon_in_context{
   my ($obj,$dummy,$value) = @_;
   if( defined $value) {
      $obj->{'is_start_exon_in_context'} = $value;
    }
    return $obj->{'is_start_exon_in_context'};

}

       
=head2 is_end_exon_in_context

 Title   : is_end_exon_in_context
 Usage   : $obj->is_end_exon_in_context($newval)
 Function: 
 Example : 
 Returns : value of is_end_exon_in_context
 Args    : newvalue (optional)


=cut

sub is_end_exon_in_context{
   my ($obj,$dummy,$value) = @_;
   if( defined $value) {
      $obj->{'is_end_exon_in_context'} = $value;
    }
    return $obj->{'is_end_exon_in_context'};

}


=head2 end_exon_rank

 Title   : end_exon_rank
 Usage   : $obj->end_exon_rank($newval)
 Function: 
 Example : 
 Returns : value of end_exon_rank
 Args    : newvalue (optional)


=cut

sub end_exon_rank{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'end_exon_rank'} = $value;
    }
    return $obj->{'end_exon_rank'};

}




