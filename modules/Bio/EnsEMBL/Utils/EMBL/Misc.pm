use strict;

package Bio::EnsEMBL::Utils::EMBL::Misc;

use Exporter;

use vars qw(@ISA @EXPORT_OK);

@ISA = qw(Exporter);

@EXPORT_OK = qw(features2join_string);



=head2 features2join_string

  Arg [1]    : listref of Bio::EnsEMBL::SeqFeatures
  Example    : $location = features2join_string(\@features);
  Description: Constructs an EMBL location string from a list of features
  Returntype : string
  Exceptions : none
  Caller     : EMBL dumping scripts

=cut

sub features2join_string {
  my $features = shift;
  
  my @join = ();

  foreach my $f (@$features) {
    my $ctg = $f->contig;

    if($ctg->isa('Bio::EnsEMBL::Slice') && 
       $f->start >= 1 && $f->end <= $ctg->length) {
      #this feature in on a slice and doesn't lie outside the boundary
      
      if($f->strand() == 1) {
	push @join, $f->start()."..".$f->end();
      } else {
	push @join, "complement(".$f->start()."..".$f->end().")";
      }
    } else {
      my @fs = ();
      #this feature is outside the boundary of the dump, 
      #convert to clone coords
      if($ctg->isa('Bio::EnsEMBL::Slice')) {
	#copy the feature, and transform to contig coords
	my $new_f;
	%$new_f = %$f;
	bless $new_f, ref $f;
	push @fs, $new_f->transform;
      }

      if($fs[0]->isa('Bio::EnsEMBL::StickyExon')) {
	@fs = @{$fs[0]->get_all_component_Exons};
      }

      #use the accession:x-y format
      foreach my $feat (@fs) {
	my $contig = $feat->contig;
	my $clone_acc = $contig->clone->id;
	my $clone_ver = $contig->clone->embl_version;
	my $start = $feat->start + $contig->embl_offset - 1;
	my $end   = $feat->end + $contig->embl_offset - 1;
	my $strand = $feat->strand;
	if($strand == 1) {
	  push @join, "$clone_acc.$clone_ver:$start..$end";
	} else {
	  push @join, "complement($clone_acc.$clone_ver:$start..$end)";
	}
      }
    }
  }

  my $out = join ',', @join;

  if(scalar @join > 1) {
    $out = "join($out)";
  }

  return $out;
}


1;
