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

# Parse UniProt (SwissProt & SPTrEMBL) files to create xrefs.
#
# Files actually contain both types of xref, distinguished by ID line;
#
# ID   CYC_PIG                 Reviewed;         104 AA.  Swissprot
# ID   Q3ASY8_CHLCH            Unreviewed;     36805 AA.  SPTrEMBL



package XrefParser::UniProtParser_descriptions_only;

use strict;
use warnings;
use Carp;
use POSIX qw(strftime);
use File::Basename;

use base qw( XrefParser::BaseParser );


# --------------------------------------------------------------------------------

sub run {


  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $files        = $ref_arg->{files};
  my $release_file = $ref_arg->{rel_file};
  my $verbose      = $ref_arg->{verbose};
  my $dbi          = $ref_arg->{dbi};
  $dbi = $self->dbi unless defined $dbi;

  if((!defined $source_id) or (!defined $species_id) or (!defined $files) or (!defined $release_file)){
    croak "Need to pass source_id, species_id, files and rel_file as pairs";
  }
  $verbose |=0;

  my $file = @{$files}[0];

  my ( $sp_source_id, $sptr_source_id, $sp_release, $sptr_release );

  $sp_source_id =
    $self->get_source_id_for_source_name('Uniprot/SWISSPROT',"sequence_mapped", $dbi);
  $sptr_source_id =
    $self->get_source_id_for_source_name('Uniprot/SPTREMBL', undef, $dbi);

  print "SwissProt source id for $file: $sp_source_id\n" if($verbose);
  print "SpTREMBL source id for $file: $sptr_source_id\n" if($verbose);
 

  my @xrefs =
    $self->create_xrefs( $sp_source_id, $sptr_source_id, $species_id, $file, $verbose, $dbi );

  if ( !@xrefs ) {
      return 1;    # 1 error
  }

#  # delete previous if running directly rather than via BaseParser
#  if (!defined(caller(1))) {
#    print "Deleting previous xrefs for these sources\n" if($verbose);
#    $self->delete_by_source(\@xrefs);
#  }

  # upload
  if(!defined($self->upload_xref_object_graphs(@xrefs, $dbi))){
    return 1; 
  }

    if ( defined $release_file ) {
        # These two lines are duplicated from the create_xrefs() method
        # below...
        my $sp_pred_source_id =
          $self->get_source_id_for_source_name(
            'Uniprot/SWISSPROT_predicted', undef, $dbi);
        my $sptr_pred_source_id =
          $self->get_source_id_for_source_name(
            'Uniprot/SPTREMBL_predicted', undef, $dbi);

        # Parse Swiss-Prot and SpTrEMBL release info from
        # $release_file.
        my $release_io = $self->get_filehandle($release_file);
        while ( defined( my $line = $release_io->getline() ) ) {
            if ( $line =~ m#(UniProtKB/Swiss-Prot Release .*)# ) {
                $sp_release = $1;
                print "Swiss-Prot release is '$sp_release'\n" if($verbose);
            } elsif ( $line =~ m#(UniProtKB/TrEMBL Release .*)# ) {
                $sptr_release = $1;
                print "SpTrEMBL release is '$sptr_release'\n" if($verbose);
            }
        }
        $release_io->close();

        # Set releases
        $self->set_release( $sp_source_id,        $sp_release, $dbi );
        $self->set_release( $sptr_source_id,      $sptr_release, $dbi );
        $self->set_release( $sp_pred_source_id,   $sp_release, $dbi );
        $self->set_release( $sptr_pred_source_id, $sptr_release, $dbi );
    }


  return 0; # successfull
}


# --------------------------------------------------------------------------------
# Parse file into array of xref objects

sub create_xrefs {
  my ($self, $sp_source_id, $sptr_source_id, $species_id, $file, $verbose, $dbi ) = @_;

  my $num_sp = 0;
  my $num_sptr = 0;
  my $num_sp_pred = 0;
  my $num_sptr_pred = 0;


  # Get predicted equivalents of various sources used here
  my $sp_pred_source_id =
    $self->get_source_id_for_source_name('Uniprot/SWISSPROT_predicted', undef, $dbi);
  my $sptr_pred_source_id =
    $self->get_source_id_for_source_name('Uniprot/SPTREMBL_predicted', undef, $dbi);

  print "Predicted SwissProt source id for $file: $sp_pred_source_id\n" if($verbose);
  print "Predicted SpTREMBL source id for $file: $sptr_pred_source_id\n" if($verbose);

  my $uniprot_io = $self->get_filehandle($file);
  if ( !defined $uniprot_io ) { return }

  my @xrefs;

  local $/ = "//\n";

  # Create a hash of all valid taxon_ids for this species
  my %species2tax = $self->species_id2taxonomy($dbi);
  my @tax_ids = @{$species2tax{$species_id}};
  my %taxonomy2species_id = map{ $_=>$species_id } @tax_ids;

  while ( $_ = $uniprot_io->getline() ) {

    # if an OX line exists, only store the xref if the taxonomy ID that the OX
    # line refers to is in the species table
    # due to some records having more than one tax_id, we need to check them 
    # all and only proceed if one of them matches.
    #OX   NCBI_TaxID=158878, 158879;
    #OX   NCBI_TaxID=103690;

    my ($ox) = $_ =~ /OX\s+[a-zA-Z_]+=([0-9 ,]+)/;
    my @ox = ();
    my $found = 0;

    if ( defined $ox ) {
      @ox = split /\, /, $ox;

      # my %taxonomy2species_id = $self->taxonomy2species_id();

      foreach my $taxon_id_from_file (@ox) {
        $taxon_id_from_file =~ s/\s//;
	if ( exists $taxonomy2species_id{$taxon_id_from_file} ){
	  $found = 1;
	}
      }
    }

    next if (!$found); # no taxon_id's match, so skip to next record
    my $xref;

    # set accession (and synonyms if more than one)
    # AC line may have primary accession and possibly several ; separated synonyms
    # May also be more than one AC line
    my ($acc) = $_ =~ /(AC\s+.+)/s; # will match first AC line and everything else
    my @all_lines = split /\n/, $acc;

    # extract ^AC lines only & build list of accessions
    my @accessions;
    foreach my $line (@all_lines) {
      my ($accessions_only) = $line =~ /^AC\s+(.+)/;
      push(@accessions, (split /;\s*/, $accessions_only)) if ($accessions_only);
    }

    $xref->{ACCESSION} = $accessions[0];
    for (my $a=1; $a <= $#accessions; $a++) {
      push(@{$xref->{"SYNONYMS"} }, $accessions[$a]);
    }

    # Check for CC (caution) lines containing certain text
    # if this appears then set the source of this and and dependent xrefs to the predicted equivalents
    my $is_predicted = /CC.*EMBL\/GenBank\/DDBJ whole genome shotgun \(WGS\) entry/;

    my ($label, $sp_type) = $_ =~ /ID\s+(\w+)\s+(\w+)/;

    # SwissProt/SPTrEMBL are differentiated by having STANDARD/PRELIMINARY here
    if ($sp_type =~ /^Reviewed/i) {

      $xref->{SOURCE_ID} = $sp_source_id;
      if ($is_predicted) {
	$xref->{SOURCE_ID} = $sp_pred_source_id;
	$num_sp_pred++;
      } else {
	$xref->{SOURCE_ID} = $sp_source_id;
	$num_sp++;
      }
    } elsif ($sp_type =~ /Unreviewed/i) {

      if ($is_predicted) {
	$xref->{SOURCE_ID} = $sptr_pred_source_id;
	$num_sptr_pred++;
      } else {
	$xref->{SOURCE_ID} = $sptr_source_id;
	$num_sptr++;
      }

    } else {

      next; # ignore if it's neither one nor t'other

    }



    # some straightforward fields
    $xref->{LABEL} = $label;
    $xref->{SPECIES_ID} = $species_id;
    $xref->{SEQUENCE_TYPE} = 'peptide';
    $xref->{STATUS} = 'experimental';

    # May have multi-line descriptions
    my ($description_and_rest) = $_ =~ /(DE\s+.*)/s;
    @all_lines = split /\n/, $description_and_rest;

    # extract ^DE lines only & build cumulative description string
    my $description = " ";
    my $name        = "";
    my $sub_description = "";
    foreach my $line (@all_lines) {

      next if(!($line =~ /^DE/));


      # get the data
      if($line =~ /^DE   RecName: Full=(.*);/){
        $name .= $1;
      }
      elsif($line =~ /RecName: Full=(.*);/){
        $description .= $1;
      }
      elsif($line =~ /SubName: Full=(.*);/){
        $sub_description .= $1;
      }


      $description =~ s/^\s*//g;
      $description =~ s/\s*$//g;

      my $desc = $name.$description;
      if(!length($desc)){
	$desc = $sub_description;
      }
      $desc =~ s/\(\s*EC\s*\S*\)//g;
      $xref->{DESCRIPTION} = $desc;


      push @xrefs, $xref;
    
    }
  }

  $uniprot_io->close();

  print "Read $num_sp SwissProt xrefs and $num_sptr SPTrEMBL xrefs from $file\n" if($verbose);
  print "Found $num_sp_pred predicted SwissProt xrefs and $num_sptr_pred predicted SPTrEMBL xrefs\n" if (($num_sp_pred > 0 || $num_sptr_pred > 0) and $verbose);
  
  return \@xrefs;

  #TODO - currently include records from other species - filter on OX line??
}

1;
