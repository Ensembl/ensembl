package XrefMapper::Interpro;

use vars '@ISA';
@ISA = qw{ XrefMapper::BasicMapper };

use strict;
use warnings;
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

  my $add_object_xref_sth = $self->xref->dbc->prepare('insert into object_xref (object_xref_id, ensembl_id,ensembl_object_type, xref_id, linkage_type, ox_status, master_xref_id ) values (?, ?, ?, ?, ?, "DUMP_OUT", ?)');

  local $add_object_xref_sth->{RaiseError}; #catch duplicates
  local $add_object_xref_sth->{PrintError}; # cut down on error messages
  
  my $add_go_xref_sth = $self->xref->dbc->prepare('insert into go_xref (object_xref_id, linkage_type) values (?, ?)'); 

  my $ins_ix_sth = $self->xref->dbc->prepare("insert into identity_xref (object_xref_id, query_identity, target_identity) values(?, 100, 100)");

  # Get a mapping of protein domains to ensembl translations for
  # interpro dependent xrefs
  my $core_sql = "SELECT hit_name, translation_id FROM protein_feature" ;
  my $core_sth = $self->core->dbc->prepare($core_sql);
  $core_sth->execute();
  my %domain_to_translation = ();
  my ($domain, $translation);
  $core_sth->bind_columns(\$domain, \$translation);
  while ($core_sth->fetch()) {
    $domain_to_translation{$domain} ||= [];
    push @{$domain_to_translation{$domain}}, $translation;
  }

  my $dep_sth    = $self->xref->dbc->prepare("select dependent_xref_id, linkage_annotation from dependent_xref where master_xref_id = ?");
  
  # Get a list of interpro data, including dependent xrefs if avail
  $sth = $self->xref->dbc->prepare("
    SELECT ip.interpro, ip.pfam, x2.xref_id, x2.source_id,
           dx.linkage_annotation, dx.master_xref_id
      FROM interpro ip, xref x
        LEFT JOIN dependent_xref dx ON x.xref_id=dx.master_xref_id
          LEFT JOIN xref x2 ON dx.dependent_xref_id=x2.xref_id
            WHERE ip.interpro = x.accession");
  my $rv = $sth->execute();
#  my %interpro_cache;
  my %added;
  my $dup=0;
  while( my $row = $sth->fetchrow_arrayref() ){
    my ( $interpro, $pfam, $dx_xref_id, $dx_source_id, $go_linkage, $master_id ) = @$row;
    if( $dx_xref_id ){
      foreach my $ensembl_id( @{$domain_to_translation{$pfam}||[]} ){
        #...And the interpro domain maps to a translation
	$add_object_xref_sth->execute($object_xref_id, $ensembl_id, 'Translation', $dx_xref_id, 'DEPENDENT', $master_id);	  
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
	$added{$dx_source_id}++;
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
	    $add_object_xref_sth->execute($object_xref_id, $ensembl_id, 'Translation', $dep_xref_id, 'DEPENDENT', $new_master_id);	  
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
  
  
  print "\n".$dup." already existed\n\n" if($self->verbose);

  print("  Wrote $ipro_count interpro table entries\n") if($self->verbose);
  print("    including $oxref_count object xrefs, \n") if($self->verbose);
  print("    and $goxref_count go xrefs\n") if($self->verbose);
#  foreach my $key (keys %added){
#    print "id= $key has ".$added{$key}. " object xrefs added\n";
#  }
}

1;
