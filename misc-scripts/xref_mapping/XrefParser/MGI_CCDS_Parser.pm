=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

package XrefParser::MGI_CCDS_Parser;

use strict;
use warnings;
use Carp;
use DBI;

use base qw( XrefParser::BaseParser );


sub run_script {

  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $file         = $ref_arg->{file};
  my $verbose      = $ref_arg->{verbose};

  if((!defined $source_id) or (!defined $species_id) or (!defined $file) ){
    croak "Need to pass source_id, species_id and file as pairs";
  }
  $verbose |=0;

  my $wget = "";

  if($file =~ /wget[=][>](\S+?)[,]/){
    $wget = $1;
  }


  my %label;
  my %version;
  my %description;
  my %accession;

  my $dbi = $self->dbi();

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



  #
  # Get master xref ids via the ccds label.
  #

  $sql = 'select x.label, x.xref_id from xref x, source s where x.source_id = s.source_id and s.name ="CCDS"';
  
  my %ccds_label_to_xref_id;
  $sth = $dbi->prepare($sql);
  $sth->execute();
  my ($xref_id);
  $sth->bind_columns(\$lab, \$xref_id);
  while (my @row = $sth->fetchrow_array()) {
    $ccds_label_to_xref_id{$row[0]} = $row[1];
  }
  $sth->finish;



  my $ua = LWP::UserAgent->new();
  $ua->timeout(10);
  $ua->env_proxy();
  

  my $count = 0;
  my $ccds_missing = 0;
  my $entrezgene_missing = 0;

  my $response = $ua->get($wget);
  
  if ( !$response->is_success() ) {
    die $response->status_line;
  }
  else{
    #
    #
    ##chromosome	g_accession	gene	gene_id	ccds_id	ccds_status	cds_strand	cds_from	cds_to	cds_locations	match_type
    #1	NC_000067.5	Xkr4	497097	CCDS14803.1	Public	-	3206102	3661428	[3206102-3207048, 3411782-3411981, 3660632-3661428]	Identical
    #1	NC_000067.5	Rp1h	19888	CCDS14804.1	Public	-	4334680	4342905	[4334680-4340171, 4341990-4342161, 4342282-4342905]	Identical
    my @lines = split(/\n/,$response->content);
    foreach my $line (@lines){
      my($chrom, $g_acc, $gene_name, $entrez_id, $ccds, @junk) = split(/\t/,$line);
      if(defined($ccds_label_to_xref_id{$ccds})){ 
	if(defined($accession{$gene_name}) and
	   defined($label{$accession{$gene_name}})){
	  my $acc = $accession{$gene_name};
	  $self->add_dependent_xref({ master_xref_id => $ccds_label_to_xref_id{$ccds},
				      acc            => $acc,
				      version        => $version{$acc},
				      label          => $label{$acc},
				      desc           => $description{$acc},
				      source_id      => $source_id,
				      species_id     => $species_id });

	  $count++;
	}
	else{
	  $entrezgene_missing++;
	}
      }
      else{
	$ccds_missing++;
      }
    }
  }
  print "$ccds_missing ccds not resolved, $entrezgene_missing mgi not found. Added $count MGI xrefs via CCDS\n" if($verbose);
  return 0;
}

1;

