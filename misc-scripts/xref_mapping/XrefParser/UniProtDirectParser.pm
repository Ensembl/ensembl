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

package XrefParser::UniProtDirectParser;

use strict;
use warnings;
use Carp;
use DBI;

use base qw( XrefParser::BaseParser );
use XrefParser::Database;

# Parse file of Uniprot records and assign direct xrefs
# All assumed to be linked to translation


# --------------------------------------------------------------------------------

sub run_script {

 my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $file         = $ref_arg->{file};
  my $verbose      = $ref_arg->{verbose};
  my $dbi          = $ref_arg->{dbi};
  $dbi = $self->dbi unless $dbi;

  if((!defined $source_id) or (!defined $species_id) or (!defined $file) ){
    croak "Need to pass source_id, species_id and file as pairs";
  }
  $verbose |=0;

  my $user = "ensro";
  my $host;
  my $port;
  my $dbname;
  my $wget = "";

  if($file =~ /host[=][>](\S+?)[,]/){
    $host = $1;
  }
  if($file =~ /port[=][>](\S+?)[,]/){
    $port =  $1;
  }
  if($file =~ /dbname[=][>](\S+?)[,]/){
    $dbname = $1;
  }
  if($file =~ /wget[=][>](\S+?)[,]/){
    $wget = $1;
  }


  my $ua = LWP::UserAgent->new();
  $ua->timeout(10);
  $ua->env_proxy();

  my $response = $ua->get($wget);

  if ( !$response->is_success() ) {
    warn($response->status_line);
    return 1;
  }
 
  my $production_db =  XrefParser::Database->new({ host   => $host,
					     port   => $port,
					     user   => $user,
					     dbname => $dbname,
					     pass   => ""});
  my $prod_dbi = $production_db->dbi();

  if(!defined($prod_dbi)){
    return 1;
  }

  my ($prefix) = $prod_dbi->selectrow_array("SELECT species_prefix FROM species WHERE taxon = $species_id");

  my %prefix = ($species_id => $prefix);

  if(!defined($prefix{$species_id})){
    print "No prefix known for this species $species_id???\n";
    return 1;
  }

  my $parsed_count = 0;


  my %prot2ensembl;

  my $count = 0;

  my @lines = split(/\n/,$response->content);
  foreach my $line (@lines){
    my ($prot, $ens) = split /\s+/,$line;
    if($ens =~ /^$prefix{$species_id}P\d+$/){
      push @{$prot2ensembl{$prot}}, $ens;
   }
  }

  my $sw_source_id =  $self->get_source_id_for_source_name("uniprot/swissprot","sequence_mapped", $dbi);
  if($sw_source_id < 1){
    die "Could not find source id for uniprot/swissprot ???\n";
  }
  else{
    print "Source_id = $sw_source_id\n";
  }
  my $get_desc_sth = $dbi->prepare("select xref_id, version, label, description from xref where source_id = $sw_source_id and accession = ?");


  my $get_dependents_sth = $dbi->prepare("select dependent_xref_id, linkage_annotation, linkage_source_id  from dependent_xref where master_xref_id = ?");

  my $add_dependent_xref_sth = $dbi->prepare("INSERT INTO dependent_xref (master_xref_id,dependent_xref_id,linkage_annotation, linkage_source_id) VALUES (?,?,?,?)");


  my $get_aliases_sth =  $dbi->prepare("select synonym from synonym where xref_id = ?");
  my $add_alias_sth   =  $dbi->prepare("INSERT IGNORE INTO synonym (xref_id, synonym) VALUES (?, ?)");



  my $err_count=0;
  foreach my $key (keys %prot2ensembl){

    #
    # get the descrptions etc for the uniprot entry
    #
    $get_desc_sth->execute($key);
    my ($old_xref_id, $version, $label, $description);
    $get_desc_sth->bind_columns(\$old_xref_id, \$version, \$label, \$description);
    $get_desc_sth->fetch;
    if(!defined($old_xref_id)){
      print "Could not find $key in the database\n" if ($err_count <10);
      $err_count++;
      next;
    }
    $count++;

    #
    # get the dependents
    #
    my %linkage_anotation=();
    my %linkage_source_id=();
    my ($dependent_xref_id, $linkage_annotation, $linkage_source_id);
    $get_dependents_sth->execute($old_xref_id);
    $get_dependents_sth->bind_columns(\$dependent_xref_id, \$linkage_annotation, \$linkage_source_id);
    while($get_dependents_sth->fetch){
      $linkage_anotation{$dependent_xref_id} =  $linkage_annotation;
      $linkage_source_id{$dependent_xref_id} =  $linkage_source_id;
    }

#    print $key."\t";
    #
    # Add the new xref
    #

    my $xref_id = $self->add_xref({ acc        => $key,
				    version    => $version,
				    label      => $label,
				    desc       => $description,
				    source_id  => $source_id,
				    species_id => $species_id,
                                    dbi        => $dbi,
				    info_type  => "DIRECT"} );


    #
    # Add the synonyms
    #
    my $synonym;
    $get_aliases_sth->execute($old_xref_id);
    $get_aliases_sth->bind_columns(\$synonym);
    while($get_aliases_sth->fetch()){
      $add_alias_sth->execute($xref_id, $synonym) || croak "Could not add synonym for $xref_id, $synonym";
    }


    foreach my $trans (@{$prot2ensembl{$key}}){
      #
      #add the direct xref entry
      #

      $self->add_direct_xref( $xref_id, $trans, "Translation", '', $dbi);
#      print ":".$trans;

      #
      #add the dependents
      #
      foreach my $dep (keys %linkage_anotation){
	$add_dependent_xref_sth->execute($xref_id, $dep, $linkage_anotation{$dep}, $linkage_source_id{$dep});
      }
    }
  }


  print $count." entrys added\n".$err_count." not found\n";
  return 0;
}


1;
