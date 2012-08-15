=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::DBSQL::AttributeAdaptor - Provides database interaction for
Bio::EnsEMBL::Attribute objects.


=head1 SYNOPSIS

  # $db is a Bio::EnsEMBL::DBSQL::DBAdaptor object:
  $attribute_adaptor = $db->get_AttributeAdaptor();

  $attributes = $attribute_adaptor->fetch_all_by_MiscFeature($feature);

  $attributes = $attribute_adaptor->fetch_all_by_Slice($slice);

  $attribute_adaptor->store_on_Slice( $slice, \@attributes );

  $attribute_adaptor->store_on_MiscFeature( $misc_feature,
    \@attributes )

=head1 DESCRIPTION

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::AttributeAdaptor;

use strict;
use warnings;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Attribute;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);


=head2 new

  Arg [...]  : Superclass args.  See Bio::EnsEMBL::DBSQL::BaseAdaptor
  Description: Instantiates a Bio::EnsEMBL::DBSQL::AttributeAdaptor
  Returntype : Bio::EnsEMBL::AttributeAdaptor
  Exceptions : none
  Caller     : DBAdaptor
  Status     : Stable

=cut


sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);


  # cache creation could go here
  return $self;
}

use vars '$AUTOLOAD';

sub AUTOLOAD {
  my ($self,@args) = @_;
  my @array_return=();
  my $ref_return = undef;
  $AUTOLOAD =~ /^.*::(\w+_)+(\w+)$/ ;

  my $sub = $1;
  my $type = $2;



#  print STDERR "AUTO".$AUTOLOAD."\n";

#  print STDERR "AUTOLOAD reached with call to $sub of type $type\n";
  if($self->can($sub)){
    return $self->$sub($type,@args);
  }
  else{
    warn("In AttribAdaptor cannot call sub $sub$type\n");
  }
  return undef;
}



sub store_on_ {
  my $self = shift;
  my $type = shift;
  my $object = shift;
  my $attributes = shift;
  my $table;


  my $object_id;
  if($type =~ /[GT][er][na][en]/){
    if (!ref($object)){
      $object_id = $object;
    }
    else {
      $object_id = $object->dbID;
    }
    $table = lc($type);
#    $type = lc($type);
  }
  else{
    if(!ref($object) || !$object->isa('Bio::EnsEMBL::'.$type)) {
      throw("$type argument is required. but you passed $object");
    }
    if($type eq "Slice"){
      $object_id = $object->get_seq_region_id();
      $table = "seq_region"; 
      $type = "seq_region";
    }
    else{
      if($type eq "MiscFeature"){
	$type = "misc_feature";
	$table = "misc"; 
      }
      else{
	$table = lc($type);
      }
      
      $object_id = $object->dbID();
      my $db = $self->db();
      
      if(!$object->is_stored($db)) {
	throw("$type is not stored in this database.");
      }
      
    }
  }
  my $sth = $self->prepare( "INSERT into ".$table."_attrib ".
			    "SET ".$type."_id = ?, attrib_type_id = ?, ".
			    "value = ? " );

  my $undef_circular_cache = 0;
  for my $attrib ( @$attributes ) {
    if(!ref($attrib) && $attrib->isa('Bio::EnsEMBL::Attribute')) {
      throw("Reference to list of Bio::EnsEMBL::Attribute objects " .
            "argument expected.");
    }

#    next if ! $attrib->value;

    my $atid = $self->_store_type( $attrib );
    if ((defined $attrib->code) and ($attrib->code eq 'circular_seq')) {
	$undef_circular_cache = 1;
    }
    $sth->bind_param(1,$object_id,SQL_INTEGER);
    $sth->bind_param(2,$atid,SQL_INTEGER);
    $sth->bind_param(3,$attrib->value,SQL_VARCHAR);
    $sth->execute();
  }

  if($table eq "seq_region") {        
    if ($undef_circular_cache) {
	#the slice is circular
	$object->{'circular'} = 1;
	my $slice_adaptor = $object->adaptor();
        #undefine slice adaptor->is_circular and the circular slice cache
	if (defined $slice_adaptor) {
	    $slice_adaptor->{'is_circular'} = undef;
	    $slice_adaptor->{'circular_sr_id_cache'} = {};
	}
    }
  }

  return;
}


sub remove_from_{
  my $self   = shift;
  my $type   = shift;
  my $object = shift;
  my $code   = shift;
  my $table;

  if(!ref($object) || !$object->isa('Bio::EnsEMBL::'.$type)) {
    throw("$type argument is required or a attrib code. but you passed $object");
  }

  my $object_id;
  if($type eq "Slice"){
    $object_id = $object->get_seq_region_id();
    $table = "seq_region"; 
    $type = "seq_region";
    if ((defined $code) and ($code eq 'circular_seq')) {
	#undefine slice->is_circular, slice adaptor->is_circular and the circular slice cache
	$object->{'circular'} = undef;
	my $slice_adaptor = $object->adaptor();
	if (defined $slice_adaptor) {
	    $slice_adaptor->{'is_circular'} = undef;
	    $slice_adaptor->{'circular_sr_id_cache'} = {};
	}
    }
  }
  else{
    if($type eq "MiscFeature"){
      $type = "misc_feature";
      $table = "misc"; 
    }
    else{
      $table = lc($type);
    }

    $object_id = $object->dbID();
    my $db = $self->db();
    
    if(!$object->is_stored($db)) {
      throw("$type is not stored in this database.");
    }

  }

  if(!defined($object_id)) {
    throw("$type must have dbID.");
  }

  my $sth;
  if(defined($code)){
    $sth = $self->prepare("DELETE a FROM ".$table."_attrib a ,attrib_type at " .
                         "WHERE a.attrib_type_id = at.attrib_type_id AND ".
                         "a.".$type."_id = ? AND ".
			 "at.code like ?");
    $sth->bind_param(1,$object_id,SQL_INTEGER);
    $sth->bind_param(2,$code,SQL_VARCHAR);
  }
  else{
    $sth = $self->prepare("DELETE FROM ".$table."_attrib " .
                         "WHERE ".$type."_id = ?");
    $sth->bind_param(1,$object_id,SQL_INTEGER);
  }
  $sth->execute();

  $sth->finish();

  return;
}



sub fetch_all_by_{
  my $self     = shift;
  my $type     = shift;
  my $object   = shift;
  my $code     = shift;
  my $table =undef;

  if(!ref($object) || !$object->isa('Bio::EnsEMBL::'.$type)) {
    throw("$type argument is required. but you passed $object");
  }

  my $object_id;
  if($type eq "Slice"){
    $object_id = $object->get_seq_region_id();
    $table = "seq_region"; 
    $type = "seq_region";
  }
  else{
    if($type eq "MiscFeature"){
      $type = "misc_feature";
      $table = "misc"; 
    }
    else{
      $table = lc($type);
    }

    $object_id = $object->dbID();
  }

  if(!defined($object_id)) {
    throw("$type must have dbID.");
  }


  my $sql = "SELECT at.code, at.name, at.description, t.value " .
              "FROM ".($table||$type)."_attrib t, attrib_type at ".
                 "WHERE t.".$type."_id = ? " .
	         "AND   at.attrib_type_id = t.attrib_type_id ";

  if(defined($code)){
    $sql .= 'AND at.code like "'.$code.'" ';
  }
		   
  my $sth = $self->prepare($sql);

  $sth->bind_param(1,$object_id,SQL_INTEGER);
  $sth->execute();

  my $results = $self->_obj_from_sth($sth);

  $sth->finish();

  return $results;
  
}


sub DESTROY{
}


#
# _id_check
#
# backwards compatibility check:
# check if $ensID is an object; if so, return $obj->dbID
#

sub _id_check {
  my $self = shift;
  my $ensID = shift;

  if ($ensID =~ /^\d+$/) {
    return $ensID;
  
  } elsif (ref($ensID) eq 'Bio::EnsEMBL::Gene' or
      ref($ensID) eq 'Bio::EnsEMBL::Transcript' or
      ref($ensID) eq 'Bio::EnsEMBL::Translation') {

    warning("You should pass a dbID rather than an ensembl object to store the attribute on");

    if ($ensID->dbID) {
      return $ensID->dbID;
    } else {
      throw("Ensembl object ".$ensID->display_id." doesn't have a dbID, can't store attribute");
    }

  } else {
    throw("Invalid dbID");
  }

}


# _store_type

sub _store_type {
  my $self = shift;
  my $attrib = shift;

  my $sth1 = $self->prepare
    ("INSERT IGNORE INTO attrib_type set code = ?, name = ?, ".
     "description = ?" );


  $sth1->bind_param(1,$attrib->code,SQL_VARCHAR);
  $sth1->bind_param(2,$attrib->name,SQL_VARCHAR);
  $sth1->bind_param(3,$attrib->description,SQL_LONGVARCHAR);

  my $rows_inserted =  $sth1->execute();

  my $atid = $sth1->{'mysql_insertid'};

  if($rows_inserted == 0) {
    # the insert failed because the code is already stored
    my $sth2 = $self->prepare
      ("SELECT attrib_type_id FROM attrib_type " .
       "WHERE code = ?");
    $sth2->bind_param(1,$attrib->code,SQL_VARCHAR);
    $sth2->execute();
    ($atid) = $sth2->fetchrow_array();

    $sth2->finish();

    if(!$atid) {
      throw("Could not store or fetch attrib_type code [".$attrib->code."]\n" .
	    "Wrong database user/permissions?");
    }
  }

  $sth1->finish();

  return $atid;
}


sub _obj_from_sth {
  my $self = shift;
  my $sth = shift;

  my ($code, $name, $desc, $value);
  $sth->bind_columns(\$code, \$name, \$desc, \$value);

  my @results;
  while($sth->fetch()) {
    push @results, Bio::EnsEMBL::Attribute->new_fast
                                             ( {'code' => $code,
                                                'name' => $name,
                                                'description' => $desc,
                                                'value' => $value} );
  }

  return \@results;
}




1;

