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

package XrefParser::EG_DBParser;

use strict;
use warnings;
use Carp;
use File::Basename;

use base qw( XrefParser::BaseParser );
use XrefParser::Database;

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

  print STDERR "parsing EG_DB_Xrefs...\n";

  my ($type, $my_args) = split(/:/,$file);

  my $user  ="ensro";
  my $host;
  my $port;
  my $dbname;
  my $pass;

  if($my_args =~ /host[=][>](\S+?)[,]/){
    $host = $1;
  }
  if($my_args =~ /port[=][>](\S+?)[,]/){
    $port =  $1;
  }
  if($my_args =~ /dbname[=][>](\S+?)[,]/){
    $dbname = $1;
  }
  if($my_args =~ /pass[=][>](\S+?)[,]/){
    $pass = $1;
  }
  if($my_args =~ /user[=][>](\S+?)[,]/){
    $user = $1;
  }

  print STDERR "species_id, $species_id\n";

  my %source;

  # Todo: get the whole list of sources for All EG set
  
  my $sources_aref = ['Mycgr3_jgi_v2.0_gene', 'BROAD_U_maydis','CADRE','CADRE_Afum_A1163','AspGD','ENA_GENE','BROAD_F_oxysporum','BROAD_G_moniliformis','BROAD_G_zeae','GeneDB','BROAD_P_infestans','phyra_jgi_v1.1','physo1_jgi_v1.1','PGD_GENE','phatr_jgi_v2_bd','phatr_jgi_v2','thaps_jgi_v2','thaps_jgi_v2_bd','SCHISTODB','triad_jgi_v1.0'];
  foreach my $source_name (@$sources_aref) {
      $source{$source_name}     = $self->get_source_id_for_source_name($source_name)    || die "Could not get source_id for $source_name\n";
  }
 
  my $db =  XrefParser::Database->new({ host   => $host,
					port   => $port,
					user   => $user,
					dbname => $dbname,
					pass   => $pass});

  my $dbi2 = $db->dbi();
  if(!defined($dbi2)){
      print STDERR "failed to connect to EG_Xrefs database!\n";
    return 1;
  }

  my $sql=(<<"SQL");
  SELECT gene_stable_id, transcript_stable_id, gene_dbid, transcript_dbid, 
         source, xref_name, xref_primary_id, xref_description
    FROM EG_Xref
      WHERE taxonomy_id = $species_id
SQL

  my ($gene_stable_id, $transcript_stable_id, $gene_dbid, $transcript_dbid, $source_name, $label, $acc, $desc);
  my $sth = $dbi2->prepare($sql);
  $sth->execute() or croak( $dbi2->errstr() );
  $sth->bind_columns(\$gene_stable_id, \$transcript_stable_id, \$gene_dbid, \$transcript_dbid, \$source_name, \$label, \$acc, \$desc);
  my $added = 0;
  while ( $sth->fetch() ) {

    my ($description,$junk) = split("[[]Source:",$desc);

    if(!defined($source{$source_name})){
      print STDERR "Could not find source_id for source $source_name for xref $acc\n";
      next;
    }
    my $xref_id = $self->get_xref($acc,$source{$source_name}, $species_id);
    if(!defined($xref_id)){
      $xref_id = $self->add_xref({ acc        => $acc,
				   label      => $label,
				   desc       => $description,
				   source_id  => $source{$source_name},
				   species_id => $species_id,
				   info_type  => "DIRECT"} );
      $added++;
    }
    my $transcript_id = $transcript_dbid;
    if(defined($transcript_stable_id) and $transcript_stable_id ne ""){
      $transcript_id = $transcript_stable_id;
    }
    my $gene_id = $gene_dbid;
    if(defined($gene_stable_id) and $gene_stable_id ne ""){
	$gene_id = $gene_stable_id;
    }
    
    #$self->add_direct_xref($xref_id, $transcript_id, "Transcript", "") if (defined($transcript_id));    
    $self->add_direct_xref($xref_id, $gene_id, "Gene", "") if (defined($gene_id));  
  }
  $sth->finish;

  print "Added $added Xrefs for EGs\n" if($verbose);
  return 0;
}

1;
