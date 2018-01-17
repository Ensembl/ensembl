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

package XrefParser::MGIParser;

use strict;
use warnings;
use Carp;
use DBI;

use base qw(XrefParser::BaseParser);

sub run {

  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $files        = $ref_arg->{files};
  my $verbose      = $ref_arg->{verbose};
  my $dbi          = $ref_arg->{dbi};
  $dbi = $self->dbi unless defined $dbi;

  if((!defined $source_id) or (!defined $species_id) or (!defined $files) ){
    croak "Need to pass source_id, species_id and files as pairs";
  }
  $verbose |=0;

  my $file = @{$files}[0];

  my $file_io = $self->get_filehandle($file);
  if ( !defined $file_io ) {
    print STDERR "ERROR: Could not open $file\n";
    return 1;    # 1 is an error
  }
  my %label;
  my %version;
  my %description;
  my %accession;

  my $sql = 'select source_id, priority_description from source where name like "MGI"';
  my $sth = $dbi->prepare($sql);
  
  $sth->execute();
  my ($mgi_source_id, $desc);
  $sth->bind_columns(\$mgi_source_id, \$desc);
  my @arr;
  while($sth->fetch()){
    push @arr, $mgi_source_id;
  }
  $sth->finish;
  
  $sql = "select accession, label, version,  description from xref where source_id in (".join(", ",@arr).")";

  $sth = $dbi->prepare($sql);
  $sth->execute();
  my ($acc, $lab, $ver);
  $sth->bind_columns(\$acc, \$lab, \$ver, \$desc);
  while (my @row = $sth->fetchrow_array()) {
    if(defined($desc)){
      $accession{$lab} = $acc;
      $label{$acc} = $lab;
      $version{$acc} = $ver;
      $description{$acc} = $desc;
    }
  }
  $sth->finish;

  #synonyms
  $sql = "insert ignore into synonym (xref_id, synonym) values (?, ?)";
  my $add_syn_sth = $dbi->prepare($sql);    

  my $syn_hash = $self->get_ext_synonyms("MGI", $dbi);

  my $count = 0;
  my $syn_count = 0;
  while ( my $line = $file_io->getline() ) {
#MGI:1915941	1110028C15Rik	RIKEN cDNA 1110028C15 gene	33.61	1	ENSMUSG00000026004	ENSMUST00000042389 ENSMUST00000068168 ENSMUST00000113987 ENSMUST00000129190 ENSMUST00000132960	ENSMUSP00000036975 ENSMUSP00000063843 ENSMUSP00000109620 ENSMUSP00000118603

    if($line =~ /(MGI:\d+).*(ENSMUSG\d+)/){
      my $acc = $1;
      my $ensid = $2;
      my $xref_id = $self->add_xref({ acc        => $acc,
				      version    => $version{$acc},
				      label      => $label{$acc},
				      desc       => $description{$acc},
				      source_id  => $source_id,
                                      dbi        => $dbi,
				      species_id => $species_id,
				      info_type  => "DIRECT"} );

      $self->add_direct_xref( $xref_id, $ensid, "Gene", '', $dbi);
      if(defined($syn_hash->{$acc})){
	foreach my $syn (@{$syn_hash->{$acc}}){
	  $add_syn_sth->execute($xref_id, $syn);
	}
      }
      $count++;
    }
    else{
      print STDERR "PROBLEM: $line";
    }

  }
  print "$count direct MGI xrefs added\n";
  return 0;

}

1;
