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

package XrefMapper::Interpro;
use strict;
use warnings;

use vars '@ISA';
@ISA = qw{ XrefMapper::BasicMapper };


use XrefMapper::BasicMapper;

sub new {
  my($class, $mapper) = @_;

  my $self ={};
  bless $self,$class;
  $self->core($mapper->core);
  $self->xref($mapper->xref);
  $self->verbose($mapper->verbose);
  return $self;
}


sub process{
  my $self = shift;

  print "Writing InterPro\n" if($self->verbose);
  my( $ipro_count, $xref_count, $oxref_count, $goxref_count ) = (0,0,0,0);
  
  my $object_xref_id;
  my $sth = $self->xref->dbc->prepare("select max(object_xref_id) from object_xref");
  $sth->execute();
  $sth->bind_columns(\$object_xref_id);
  $sth->fetch();
  $sth->finish;

  $object_xref_id++;



#
# Get the sql and sths required.
#

  my $add_ox_sql = (<<"AOX");
INSERT INTO object_xref
       (object_xref_id, ensembl_id,ensembl_object_type, xref_id,
        linkage_type, ox_status, master_xref_id )
       VALUES (?, ?, ?, ?, ?, "DUMP_OUT", ?)
AOX
  my $add_object_xref_sth = $self->xref->dbc->prepare($add_ox_sql);


  my $add_go_xref_sql = (<<"AGX");
INSERT INTO go_xref (object_xref_id, linkage_type)
  VALUES (?, ?)
AGX
  my $add_go_xref_sth = $self->xref->dbc->prepare($add_go_xref_sql);

  my $ins_ix_sql = (<<"IXS");
INSERT INTO identity_xref (object_xref_id, query_identity, target_identity)
  VALUES (?, 100, 100)
IXS
  my $ins_ix_sth = $self->xref->dbc->prepare($ins_ix_sql);

  my $dep_sql = (<<"DLS");
SELECT dependent_xref_id, linkage_annotation
  FROM dependent_xref
   WHERE master_xref_id = ?
DLS
  my $dep_sth    = $self->xref->dbc->prepare($dep_sql);
  
  local $add_object_xref_sth->{RaiseError}; #catch duplicates
  local $add_object_xref_sth->{PrintError}; # cut down on error messages


  # Get a mapping of protein domains to ensembl translations for
  # interpro dependent xrefs
  my $core_sql = 'SELECT hit_name, translation_id FROM protein_feature' ;
  my $core_sth = $self->core->dbc->prepare($core_sql);
  $core_sth->execute();
  my %domain_to_translation = ();
  my ($domain, $translation);
  $core_sth->bind_columns(\$domain, \$translation);
  while ($core_sth->fetch()) {
    $domain_to_translation{$domain} ||= [];
    push @{$domain_to_translation{$domain}}, $translation;
  }


  my $int_sql = (<<"ISQ");
SELECT ip.interpro, ip.pfam, x2.xref_id, x2.source_id,
       dx.linkage_annotation, dx.master_xref_id
  FROM interpro ip, xref x
    LEFT JOIN dependent_xref dx 
         ON x.xref_id=dx.master_xref_id
    LEFT JOIN xref x2 
         ON dx.dependent_xref_id=x2.xref_id
      WHERE ip.interpro = x.accession and ip.dbtype != 'PRINTS'
ISQ
  $sth = $self->xref->dbc->prepare($int_sql);


  my $rv = $sth->execute();

#  my %added;
  my $dup=0;
  while( my $row = $sth->fetchrow_arrayref() ){
    my ( $interpro, $pfam, $dx_xref_id, $dx_source_id, $go_linkage, $master_id ) = @$row;
    if( $dx_xref_id ){
      foreach my $ensembl_id( @{$domain_to_translation{$pfam}||[]} ){
        #...And the interpro domain maps to a translation
        $add_object_xref_sth->execute($object_xref_id,
				      $ensembl_id,
				      'Translation',
				      $dx_xref_id,
				      'DEPENDENT',
				      $master_id);
	if($add_object_xref_sth->err){
	  my $err = $add_object_xref_sth->errstr;
	  if($err =~ /Duplicate/){
	    $dup++;
	    next;
	  }
	  else{
	    die "Problem adding object xref for interpro data\n";
	  }
	}
	$ins_ix_sth->execute($object_xref_id);
	$oxref_count++;
	if($go_linkage){
	  $add_go_xref_sth->execute($object_xref_id, $go_linkage );
	  $goxref_count ++;
	}



	#
	# Also add dependents of the xref and its etc...!!!
	#
	my @master_xref_ids;
	push @master_xref_ids, $dx_xref_id;
	$object_xref_id++;
	while (my $new_master_id = pop(@master_xref_ids)){
	  $dep_sth->execute($new_master_id);
	  my $dep_xref_id;
	  my $link;
	  $dep_sth->bind_columns(\$dep_xref_id, \$link);
          while($dep_sth->fetch()){
	    $add_object_xref_sth->execute($object_xref_id,
					  $ensembl_id,
					  'Translation',
					  $dep_xref_id,
					  'DEPENDENT',
					  $new_master_id);
	    if(!$add_object_xref_sth->err){
	      push @master_xref_ids, $dep_xref_id;
	      if($link){
		$add_go_xref_sth->execute($object_xref_id, $link );
	      }
	      $ins_ix_sth->execute($object_xref_id);
	    }
	    $object_xref_id++;
	  }
	}



	$object_xref_id++;
      }
    }
  }
  $sth->finish();


  if($self->verbose){
    print "\n".$dup." already existed\n\n";
    print "  Wrote $ipro_count interpro table entries\n";
    print "    including $oxref_count object xrefs, \n";
    print "    and $goxref_count go xrefs\n";
  }
  return;
}

1;
