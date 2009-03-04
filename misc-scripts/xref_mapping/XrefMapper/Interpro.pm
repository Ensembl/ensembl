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
  return $self;
}


sub process{
  my $self = shift;

  print "Writing InterPro\n";
  my( $ipro_count, $xref_count, $oxref_count, $goxref_count ) = (0,0,0,0);
  
  my $object_xref_id;
  my $sth = $self->xref->dbc->prepare("select max(object_xref_id) from object_xref");
  $sth->execute();
  $sth->bind_columns(\$object_xref_id);
  $sth->fetch();
  $sth->finish;

  $object_xref_id++;

  my $add_object_xref_sth = $self->xref->dbc->prepare('insert into object_xref (object_xref_id, ensembl_id,ensembl_object_type, xref_id, linkage_type, ox_status ) values (?, ?, ?, ?, ?, "DUMP_OUT")');
  local $add_object_xref_sth->{RaiseError}; #catch duplicates
  local $add_object_xref_sth->{PrintError}; # cut down on error messages
  
  my $add_go_xref_sth = $self->xref->dbc->prepare('insert into go_xref (object_xref_id, linkage_type) values (?, ?)'); 

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
  
  # Get a list of interpro data, including dependent xrefs if avail
  $sth = $self->xref->dbc->prepare("
    SELECT ip.interpro, ip.pfam, x2.xref_id, x2.source_id,
           dx.linkage_annotation
      FROM interpro ip, xref x
        LEFT JOIN dependent_xref dx ON x.xref_id=dx.master_xref_id
          LEFT JOIN xref x2 ON dx.dependent_xref_id=x2.xref_id
            WHERE ip.interpro = x.accession");
  my $rv = $sth->execute();
#  my %interpro_cache;
  my %added;
  my $dup=0;
  while( my $row = $sth->fetchrow_arrayref() ){
    my ( $interpro, $pfam, $dx_xref_id, $dx_source_id, $go_linkage ) = @$row;
    if( $dx_xref_id ){
      foreach my $ensembl_id( @{$domain_to_translation{$pfam}||[]} ){
        #...And the interpro domain maps to a translation
	$add_object_xref_sth->execute($object_xref_id, $ensembl_id, 'Translation', $dx_xref_id, 'DEPENDENT');	  
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
	$added{$dx_source_id}++;
	$oxref_count++;
	if($go_linkage){
	  $add_go_xref_sth->execute($object_xref_id, $go_linkage );
	  $goxref_count ++;
	}
	$object_xref_id++;
      }
    }
  }
  $sth->finish();
  
  
  print "\n".$dup." already existed\n\n";

  print("  Wrote $ipro_count interpro table entries\n");
  print("    including $oxref_count object xrefs, \n");
  print("    and $goxref_count go xrefs\n");
  foreach my $key (keys %added){
    print "id= $key has ".$added{$key}. " object xrefs added\n";
  }
}

1;
