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

package XrefParser::GOParser;

use strict;
use warnings;
use POSIX qw(strftime);
use File::Basename;
use Carp;
use base qw( XrefParser::BaseParser );

sub run {

  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $files        = $ref_arg->{files};
  my $release_file = $ref_arg->{rel_file};
  my $verbose      = $ref_arg->{verbose};

  my $gene_source_id = $self->get_source_id_for_source_name('GO_to_gene');


  if((!defined $source_id) or (!defined $species_id) or (!defined $files) ){
    croak "Needs to pass source_id, species_id and  files as pairs";
  }
  $verbose |=0;


  my $file = @{$files}[0];
  my $file_desc = @{$files}[1];

  #
  # Get the descriptions from the desc file.
  #

  my %go_to_desc;
  if($file_desc){
    my $go_desc_io = $self->get_filehandle($file_desc);

    if ( !defined $go_desc_io ) {
      print STDERR "ERROR: Could not open description file, $file_desc\n";
      return 1;    # 1 error
    }


    print "description file for GO\n" if($verbose);
    my $term = undef;
    my $desc = undef;
    while ( $_ = $go_desc_io->getline() ) {
      if(/\<id\>   # start of id tag
          (GO:\d+)  # GO: followed by the id
          \<\/id\>  # end of id tag
	  /x){
        $term = $1;
      }
      elsif(/\<name\>   # start of name tag
              (.*)       # the name we want
              \<\/name\> # end of name tag
	     /x){
        if(defined($term)){
          $go_to_desc{$term} = $1;
        }
        $term = undef;
      }
    }
    $go_desc_io->close();
  }

  my %wrongtype;

  #get the "main" GO source id.
  $source_id = $self->get_source_id_for_source_name("GO","main");


  #get the mapping that are already there so that we don't get lots of duplicates.
  # stored in the global hash xref_dependent_mapped.
  $self->get_dependent_mappings($source_id);

  my $swiss_miss=0;
  my (%swiss) = %{$self->get_valid_codes("uniprot/",$species_id)};

  my $refseq_miss=0;
  my (%refseq) = %{$self->get_valid_codes("refseq",$species_id)};

  my %mouse;
  my $mouse_set;

  my %worm;
  my $wormset;
  my $worm_separator;

  my %fish;
  my $fishset;

  my %cerevisiae;
  my $cerevisiae_set;

  my %fly;
  my $fly_set;

  my $count = 0;
  my $refseq_count = 0;
  my $uniprot_count = 0;
  my $sgd_count = 0;
  my $yeast_count = 0;
  my $worm_count = 0;
  my $mgi_count = 0;
  my $zfish_count = 0;
  my $flyb_count = 0;



  my %sp2tax     =  $self->species_id2taxonomy();   #some species have multiple
                                                    #tax_id i.e. strains

  my @tax_ids = @{$sp2tax{$species_id}};
  foreach my $tax_id ( @tax_ids){

    my $go_io = $self->get_filehandle($file);

    if ( !defined $go_io ) {
      print STDERR "ERROR: Could not open $file\n";
      return 1;    # 1 error
    }

    print "processing for taxon: $tax_id\n" if($verbose);
    my $taxon_line = "taxon:".$tax_id;
    my $miss =0;
    while ( $_ = $go_io->getline() ) {
      next if(!(/$taxon_line/x));

      chomp;
      my @array = split (/\t/x,$_);

      # Skip "NOT" terms entirely
      next if ($array[3] eq "NOT");
      my $info_text = $array[14];

      $array[9] =~ s/\'/\\\'/gx; # replace ' with \'
      my $master=0;
      if($array[0] =~ /ENSEMBL/){
        #these might be good for a check
        # match GO to Uniprot
        # match Uniprot to ENSEMBL
        # check ENSEMBL's are the same.
      }
      elsif($array[0] =~ /RefSeq/){
        if($refseq{$array[1]}) {
          foreach my $xref_id (@{$refseq{$array[1]}}){
            $self->add_dependent_xref({ master_xref_id => $xref_id,
              acc            => $array[4],
              label          => $array[4],
              desc           => $go_to_desc{$array[4]} || '',
              linkage        => $array[6],
              info_text      => $info_text,
              source_id      => $source_id,
              species_id     => $species_id} );
            $refseq_count++;
            $count++;
          }
        }
        else{
          $refseq_miss++;
        }
      }
      elsif($array[0] =~ /UniProt/){
        if($swiss{$array[1]}){
          foreach my $xref_id (@{$swiss{$array[1]}}){
            if ($array[7] =~ /Ensembl/) { next; }
            $self->add_dependent_xref({ master_xref_id => $xref_id,
              acc            => $array[4],
              label          => $array[4],
              desc           => $go_to_desc{$array[4]} || '',
              linkage        => $array[6],
              info_text      => $info_text,
              source_id      => $source_id,
              species_id     => $species_id} );
            $uniprot_count++;
            $count++;
          }
        }
        else{
          $swiss_miss++;
        }
      }
      elsif($array[0] =~ /^WB/x){
        #WB  CE20707 ZYG-9  GO:0008017  WB:WBPaper00003099|PMID:9606208 ISS  F  protein  taxon:6239  20030829  WB
        if(!defined($wormset)){
          $wormset = 1;
          %worm = %{$self->get_valid_codes("wormbase_gene",$species_id)};
        }

        my $worm_acc = $array[1];

        if(defined($worm{$worm_acc})){
          foreach my $xref_id (@{$worm{$worm_acc}}) {
            $self->add_dependent_xref({ master_xref_id => $xref_id,
                                        acc            => $array[4],
                                        label          => $array[4],
                                        desc           => $go_to_desc{$array[4]} || '',
                                        linkage        => $array[6],
                                        info_text      => $info_text,
                                        source_id      => $gene_source_id,
                                        species_id     => $species_id} );
            $worm_count++;
            $count++;
          }

        }
        else{
          $miss++;
        }
      }
      elsif($array[0] =~ /^ZFIN/x){
        #ZFIN  ZDB-GENE-030131-5418  rfng  GO:0030902  ZFIN:ZDB-PUB-050125-4|PMID:15659486  IMP  ZFIN:ZDB-MRPHLNO-050308-5  radical fringe homolog (Drosophila)  gene  taxon:7955  20050310  ZFIN
        if(!defined($fishset)){
          $fishset = 1;
          %fish = %{$self->get_valid_codes("ZFIN_ID",$species_id)};
        }
        if($fish{$array[1]}){
          foreach my $xref_id (@{$fish{$array[1]}}) {
            $self->add_dependent_xref({ master_xref_id => $xref_id,
                acc            => $array[4],
                label          => $array[4],
                desc           => $go_to_desc{$array[4]} || '',
                linkage        => $array[6],
                info_text      => $info_text,
                source_id      => $source_id,
                species_id     => $species_id} );
            $zfish_count++;
            $count++;
          }
        }
      }

      elsif($array[0] =~ /MGI/x){
        # MGI	MGI:1923501	0610007P08Rik 		GO:0004386	MGI:MGI:1354194	IEA
        # 0   1           2             3   4           5               6
        if(!defined($mouse_set)){
          $mouse_set = 1;
          # Todo: Make sure we get this hash populated
          %mouse = %{$self->get_valid_codes("MGI",$species_id)};

          print "Got " . keys (%mouse) . " MGI ids\n";

        }
        if ( $mouse{$array[1]} ){
          foreach my $xref_id ( @{$mouse{$array[1]}} ) {
            $self->add_dependent_xref({ master_xref_id => $xref_id,
                acc            => $array[4],
                label          => $array[4],
                desc           => $go_to_desc{$array[4]} || '',
                linkage        => $array[6],
                info_text      => $info_text,
                source_id      => $source_id,
                species_id     => $species_id} );
            $mgi_count++;
            $count++;
          }
        }
      }
      # SGD GO code
      elsif ($array[0] =~ /SGD/x) {

        if(!defined($cerevisiae_set)){
          $cerevisiae_set = 1;
          # Todo: Make sure we get this hash populated
          %cerevisiae = %{$self->get_valid_codes("sgd_translation",$species_id)};

          print STDERR "Got " . keys (%cerevisiae) . " cerevisiae ids\n";

        }

        if($cerevisiae{$array[1]}){
          foreach my $xref_id (@{$cerevisiae{$array[1]}}){

            my $label = $array[2];
            #print STDERR "GO SGD label: $label\n";

            # Only keep GO annotations for protein_coding genes
            # as the other annotations would get attached to transcript objects, instead of translations,
            # GO attached to Transcripts used to break the webcode display and Biomart, although not a problem aymore !?

            if (($label !~ /^t\w\(/) && ($label !~ /^\d+/) && ($label !~ /^RDN/)
            && ($label !~ /^snR/) && ($label !~ /^LSR/) && ($label !~ /^R|^T|^S|^P|^I|^H/)
            && ($label !~ /^EMT\d/) && ($label !~ /^FDH\d/) && ($label !~ /^NME\d/) && ($label !~ /^CDC\d+/)) {
              $self->add_dependent_xref({ master_xref_id => $xref_id,
                      acc            => $array[4],
                      label          => $array[4],
                      desc           => $go_to_desc{$array[4]} || '',
                      linkage        => $array[6],
                      info_text      => $info_text,
                      source_id      => $source_id,
                      species_id     => $species_id} );
              $sgd_count++;
              $count++;
            }
          }
        }
      }
      elsif ($array[0] =~ /^FB/x) {
        if(!defined($fly_set)){
          $fly_set = 1;
          %fly = %{$self->get_valid_codes("flybase_gene_id", $species_id)};
          print STDERR "Got " . keys (%fly) . " fruitfly ids\n";
        }

        my $fly_id = $array[1];
        if ( $fly{$fly_id} ) {
          foreach my $xref_id (@{$fly{$fly_id}}) {
            $self->add_dependent_xref({ master_xref_id => $xref_id,
                                        acc            => $array[4],
                                        label          => $array[4],
                                        desc           => $go_to_desc{$array[4]} || '',
                                        linkage        => $array[6],
                                        info_text      => $info_text,
                                        source_id      => $gene_source_id,
                                        species_id     => $species_id} );
            $flyb_count++;
            $count++;
          }
        }
        else {
          $miss++;
        }
      }
      elsif(!defined($wrongtype{$array[0]})){
        print STDERR "WARNING: unknown type ".$array[0]."\n" if($verbose);
        $wrongtype{$array[0]} = 1;
      }
    }

    $go_io->close();

    print "\t$count GO dependent xrefs added $refseq_miss refseq not found and $swiss_miss Swissprot not found \n" if($verbose);
    print "Added $refseq_count refseq, $uniprot_count uniprot and $zfish_count zfin\n";
  }
  if ( defined $release_file ) {
    # Parse and set release information from $release_file.
    my $release_io = $self->get_filehandle($release_file);
    my %release_hash;
    my %date_hash;
    my %species2name = $self->species_id2name();
    my @names   = @{$species2name{$species_id}};
    my %name2species_id     = map{ $_=>$species_id } @names;

    while (my $line = $release_io->getline() ) {
      my ($species, $version, $date) = split('\t', $line);
      $release_hash{$species} = $version;
      $date_hash{$species} = $date;
    }

    my $found = 0;
    foreach my $species (keys %release_hash) {
      if ($name2species_id{$species}) {
        print "GO release: " . $release_hash{$species} . "\n" if ($verbose);
        $self->set_release( $source_id, $release_hash{$species});
        $found = 1;
        last;
      }
    }
    if (!$found) {
      print "Uniprot GO release: " . $release_hash{'uniprot'} . "\n" if ($verbose);
      $self->set_release( $source_id, $release_hash{'uniprot'});
    }
  }

  return 0;
}

1;
