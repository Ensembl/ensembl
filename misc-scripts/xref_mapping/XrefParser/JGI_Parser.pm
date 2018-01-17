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

package XrefParser::JGI_Parser;

use strict;
use warnings;
use Carp;
use File::Basename;

use base qw( XrefParser::BaseParser );

# JGI protein file with gene predictons  - FASTA FORMAT  
#
#
# This is the parser that provides most functionality, subclasses 
# (JGI_ProteinParser) just set sequence type)

sub run {

  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $files        = $ref_arg->{files};
  my $release_file = $ref_arg->{rel_file};
  my $verbose      = $ref_arg->{verbose};
  my $dbi          = $ref_arg->{dbi};
  $dbi = $self->dbi unless defined $dbi;

  if((!defined $source_id) or (!defined $species_id) or (!defined $files) ){
    croak "Need to pass source_id, species_id and files as pairs";
  }
  $verbose |=0;

  my $file = @{$files}[0];

  
  my $source_name = $self->get_source_name_for_source_id ($source_id, $dbi) ;
  # the source name defines how to parse the header 

  # different formats for different sources (all have entries in external_db.txt and populate_metadata.sql )
  # 
  #
  # SOURCES  :  cint_jgi_v1  AND cint_aniseed_jgi_v1 
  # -------------------------------------------------
  #
  #  filename :  ciona.prot.fasta 
  #  >ci0100130001
  #  MLPIVDFKQCRPSVEASDKEINETAKLLVDALSTVGFAYLKNCGIKKNCRRSQKHRG*MGGVRYLYYPPI
  #  RVNIPDDEVKRNSIRRSIGYFVFPDDDVVINQPLQFKGDADVPDPVKDPITALKYIQQKLSHTCQNT*
  # 
  #
  #
  # SOURCES :  cint_jgi_v2  && cint_aniseed_jgi_v2  
  # ------------------------------------------------
  #
  # filename : FM1.aa.fasta
  #   >jgi|Cioin2|201001|fgenesh3_pm.C_chr_01p000019
  #  MQQQQQDDLVVKLVLVGDGGVGKTTFVKRHLTGEFEKKYVATLGVEVHPIVFQTQRGRIRFNVWDTAGQE
  #  DEDDDL*
  #



  my @xrefs;

  local $/ = "\n>";

  my $file_io = $self->get_filehandle($file);

  if ( !defined $file_io ) {
    print STDERR "ERROR: Could not open $file\n";
    return 1;    # 1 is an error
  }

  while ( $_ = $file_io->getline() ) {

    next if (/^File:/);   # skip header

    my $xref;

    my ($header, $sequence) = $_ =~ /^>?(.+?)\n([^>]*)/s or warn("Can't parse FASTA entry: $_\n"); 

      my @attr = split/\|/, $header ; 
      my $acession = $header ;
    
    # split header in different ways according to source name : 
    my ( $version,$label  ) ;      
    if ($source_name=~m/cint_jgi_v1/) { 
      # header format is  >ci0100146277
      # we want  146277
      ($acession = $header)  =~s/\w{6}//;  
      $version = "JGI 1.0" ;  
    
    } elsif ($source_name=~m/cint_aniseed_.*v1/) { 
      # header format is  >ci0100146277, we want this 
      $version = "JGI 1.0" ;  


    } elsif ($source_name=~m/cint_jgi_v2/) { 
      $acession = $attr[2]  ; 
      $version = "JGI 2.0" ;  
      $label = $attr[3] ; 
    } elsif ($source_name=~m/cint_aniseed_.*v2/) { 
      my $aniseed_prefix = "ci0200" ; 
      $acession = $aniseed_prefix . $attr[2]  ; 
      $version = "JGI 2.0" ;  

    }else { 
      print STDERR "WARNING : The source-name specified ($source_name) in the populate_metatable.sql file is\n" .
        "WARNING : not matching the differnt cases specified in JGI_Parser.pm - plese\n".
          "WARNING : edit the parser \n" ; 
      return 1;
    } 
    #print "ACCESSION $acession\n" ;  

    # make sequence into one long string
    $sequence =~ s/\n//g;

    # build the xref object and store it
    $xref->{ACCESSION}     = $acession;
    $xref->{LABEL}         = $acession;
    $xref->{SEQUENCE}      = $sequence;
    $xref->{SOURCE_ID}     = $source_id;
    $xref->{SPECIES_ID}    = $species_id;
    $xref->{SEQUENCE_TYPE} = $self->get_sequence_type();
    $xref->{STATUS}        = 'experimental';

    # pull cg_name from peptide files as well and create dependent xrefs
#    if ($self->get_sequence_type() =~ /peptide/) {
#      my ($cg_name) = $cg =~ /cg_name=(.*)/;
#      my %dep;
#      $dep{SOURCE_NAME} = 'JGI__Gene';
#      $dep{LINKAGE_SOURCE_ID} = $xref->{SOURCE_ID};
#      $dep{SOURCE_ID} = $celera_gene_source_id;
#      $dep{ACCESSION} = $cg_name;
#      push @{$xref->{DEPENDENT_XREFS}}, \%dep; # array of hashrefs
#    }
#
    push @xrefs, $xref;

  }

  $file_io->close();

  print scalar(@xrefs) . " JGI_ xrefs succesfully parsed\n" if($verbose);

  $self->upload_xref_object_graphs(\@xrefs, $dbi);

  return 0; # successful
}

1;
