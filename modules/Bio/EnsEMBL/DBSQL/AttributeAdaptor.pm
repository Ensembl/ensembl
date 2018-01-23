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



=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

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
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);

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



## Storing methods

## Batch store

sub store_batch_on_Object {
  my ($self, $table, $attributes, $batch_size) = @_;
  # $attributes is a hashref where the key is the object ID
  # and the value is an array ref of Attribute objects

  # maintain a hash of attrib type IDs by code so we don't have to keep looking them up...
  my $attrib_type_ids = {};

  # create an arrayref of the values to store
  my $rows = [];
  $batch_size ||= scalar(values(%$attributes));

  while (my ($obj_id, $attribs) = each %{$attributes}) {
    for my $attrib (@{$attribs}) {

      my $attrib_type_id = $attrib_type_ids->{$attrib->code()};
      if (!defined $attrib_type_id) {
        $attrib_type_id = $self->_store_type($attrib);
        $attrib_type_ids->{$attrib->code()} = $attrib_type_id;
      }

      push @$rows, sprintf('(%d, %d, %s)', $obj_id, $attrib_type_id, $self->dbc()->db_handle()->quote($attrib->value()));

      if (scalar(@$rows) == $batch_size) {
        $rows = $self->_store_batch_rows($table, $rows);
      }
    }
  }

  $rows = $self->_store_batch_rows($table, $rows);

  return;
}

sub _store_batch_rows {
  my ($self, $table, $rows) = @_;
  if (scalar(@$rows) > 0) {
        $self->dbc()->sql_helper()->execute_update(-SQL => 'INSERT INTO ' . $table . '_attrib(' . $table . '_id, attrib_type_id, value) VALUES ' . join(',', @$rows));
  }
  return [];
}


sub store_batch_on_MiscAttrib {
  my ($self, $attributes, $batch_size) = @_;

  $self->store_batch_on_Object('misc_feature', $attributes, $batch_size);

  return;
}

sub store_batch_on_Slice {
  my ($self, $attributes, $batch_size) = @_;

  $self->store_batch_on_Object('seq_region', $attributes, $batch_size);

  return;
}

sub store_batch_on_Gene {
  my ($self, $attributes, $batch_size) = @_;

  $self->store_batch_on_Object('gene', $attributes, $batch_size);

  return;
}

sub store_batch_on_Transcript {
  my ($self, $attributes, $batch_size) = @_;

  $self->store_batch_on_Object('transcript', $attributes, $batch_size);

  return;
}

sub store_batch_on_Translation {
  my ($self, $attributes, $batch_size) = @_;

  $self->store_batch_on_Object('translation', $attributes, $batch_size);

  return;
}

sub store_batch_on_DnaDnaAlignFeature {
  my ($self, $attributes, $batch_size) = @_;

  $self->store_batch_on_Object('dna_align_feature', $attributes, $batch_size);

  return;
}


## Single store

sub store_on_Object {
  my ($self, $object_id, $attributes, $table, $type) = @_;

  my $db = $self->db();
  if (!defined $type) {
    $type = $table;
  }

  my $insert_ignore = $self->insert_ignore_clause();
  my $sth = $self->prepare( "${insert_ignore} INTO ${table}_attrib (${type}_id, attrib_type_id, value)" .
                            "VALUES (?, ?, ?)" );

  for my $attrib ( @$attributes ) {
    if(!ref($attrib) && $attrib->isa('Bio::EnsEMBL::Attribute')) {
      throw("Reference to list of Bio::EnsEMBL::Attribute objects " .
            "argument expected.");
    }

    my $atid = $self->_store_type( $attrib );
    $sth->bind_param(1,$object_id,SQL_INTEGER);
    $sth->bind_param(2,$atid,SQL_INTEGER);
    $sth->bind_param(3,$attrib->value,SQL_VARCHAR);
    $sth->execute();
  }

  return;
}


sub store_on_MiscFeature {
  my ($self, $object, $attributes) = @_;

  assert_ref( $object, 'Bio::EnsEMBL::MiscFeature');

  my $object_id = $object->dbID();
  $self->store_on_Object($object_id, $attributes, 'misc', 'misc_feature');

  return;
}

sub store_on_Slice {
  my ($self, $object, $attributes) = @_;

  assert_ref( $object, 'Bio::EnsEMBL::Slice');

  my $object_id = $object->get_seq_region_id();
  $self->store_on_Object($object_id, $attributes, 'seq_region');

  my $undef_circular_cache = 0;
  for my $attrib ( @$attributes ) {
    if ((defined $attrib->code) and ($attrib->code eq 'circular_seq')) {
        $undef_circular_cache = 1;
    }
  }

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

  return;
}

sub store_on_Gene {
  my ($self, $object, $attributes) = @_;

  my $object_id;
  if (!ref($object)){      
    $object_id = $object;
  }
  else {
    $object_id = $object->dbID;    
  }

  $self->store_on_Object($object_id, $attributes, 'gene');
 
  return;
}

sub store_on_Transcript {
  my ($self, $object, $attributes) = @_;

  my $object_id;
  if (!ref($object)){      
    $object_id = $object;
  }
  else {
    $object_id = $object->dbID;    
  }

  $self->store_on_Object($object_id, $attributes, 'transcript');

  return;
}

sub store_on_Translation {
  my ($self, $object, $attributes) = @_;

  my $object_id;
  if (!ref($object)){
    $object_id = $object;
  }
  else {
    $object_id = $object->dbID;
  }

  $self->store_on_Object($object_id, $attributes, 'translation');

  return;
}

sub store_on_DnaDnaAlignFeature {
  my ($self, $object, $attributes) = @_;

  my $object_id;
  if (!ref($object)){
    $object_id = $object;
  }
  else {
    $object_id = $object->dbID;
  }

  $self->store_on_Object($object_id, $attributes, 'dna_align_feature');

  return;
}


## Remove methods

sub remove_from_Object {
  my ($self, $object_id, $table, $code, $type) = @_;

  my $db = $self->db();
  
  if(!defined($object_id)) {
    throw("$table must have dbID.");
  } 
  if (!defined($type)) {
    $type = $table;
  }
  
  my $sth;
  if (defined($code)) {
    if ($db->dbc->driver() eq 'mysql') {
      $sth = $self->prepare("DELETE a FROM " . $table . "_attrib a, attrib_type at " .
                            "WHERE a.attrib_type_id = at.attrib_type_id AND ".
                            "a." . $type . "_id = ? AND ".
                            "at.code like ?");
    } else {
      $sth = $self->prepare(qq{DELETE FROM ${table}_attrib
                                WHERE ${type}_id = ? AND
                                       attrib_type_id IN
                               (SELECT attrib_type_id
                                  FROM attrib_type
                                 WHERE code LIKE ? ) });
    }
    $sth->bind_param(1,$object_id,SQL_INTEGER);
    $sth->bind_param(2,$code,SQL_VARCHAR);
  }
  else{
    $sth = $self->prepare("DELETE FROM " . $table . "_attrib " .
                         "WHERE " . $type . "_id = ?");
    $sth->bind_param(1,$object_id,SQL_INTEGER);
  }
  $sth->execute();
  $sth->finish();
  
  return;

}


sub remove_from_MiscFeature {
  my ($self, $object, $code) = @_;

  assert_ref($object, 'Bio::EnsEMBL::MiscFeature');

  my $object_id = $object->dbID();
  $self->remove_from_Object($object_id, 'misc', $code, 'misc_feature');

  return;

}

sub remove_from_Slice {
  my ($self, $object, $code) = @_;

  assert_ref($object, 'Bio::EnsEMBL::Slice');

  my $object_id = $object->get_seq_region_id();
  $self->remove_from_Object($object_id, 'seq_region', $code);

  return;

}

sub remove_from_Gene {
  my ($self, $object, $code) = @_;

  assert_ref($object, 'Bio::EnsEMBL::Gene');

  my $object_id = $object->dbID();
  $self->remove_from_Object($object_id, 'gene', $code);

  return;

}

sub remove_from_Transcript {
  my ($self, $object, $code) = @_;

  assert_ref($object, 'Bio::EnsEMBL::Transcript');

  my $object_id = $object->dbID();
  $self->remove_from_Object($object_id, 'transcript', $code);
  
  return;

}

sub remove_from_Translation {
  my ($self, $object, $code) = @_;

  assert_ref($object, 'Bio::EnsEMBL::Translation');

  my $object_id = $object->dbID();
  $self->remove_from_Object($object_id, 'translation', $code);
  
  return;

}

sub remove_from_DnaDnaAlignFeature {
  my ($self, $object, $code) = @_;

  assert_ref($object, 'Bio::EnsEMBL::DnaDnaAlignFeature');

  my $object_id = $object->dbID();
  $self->remove_from_Object($object_id, 'dna_align_feature', $code);
  
  return;

}



## Fetch methods

sub fetch_all {
  throw("Use of method fetch_all not supported for attributes");
}

sub fetch_by_code {
  my ($self, $code) = @_;

  my $sql = "SELECT attrib_type_id, code, name, description " .
            "FROM attrib_type 
             WHERE code = '$code' ";
  my $sth = $self->prepare($sql);
  $sth->execute();

  my ($attrib_type_id, $name, $desc);
  $sth->bind_columns(\$attrib_type_id, \$code, \$name, \$desc);

  my @results = $sth->fetchrow_array;
  $sth->finish();

  return \@results;
}



sub fetch_all_by_Object {
  my ($self, $object_id, $table, $code, $type) = @_;

  if (!defined $type) {
    $type = $table;
  }

  my $sql = "SELECT at.code, at.name, at.description, a.value " .
            "FROM " . $table . "_attrib a, attrib_type at "  .
            "WHERE at.attrib_type_id = a.attrib_type_id ";

  if (defined($code)) {
        $sql .= 'AND at.code like "' . $code . '" ';
  } 
  if (defined($object_id)) {
        $sql .= "AND a." .$type . "_id = " . $object_id;
  }

  my $sth = $self->prepare($sql);
  $sth->execute();
  my $results = $self->_obj_from_sth($sth);
  $sth->finish();

  return $results;
}


sub fetch_all_by_MiscFeature {
  my ($self, $object, $code) = @_;

  my $object_id;

  if (defined($object)) {
    assert_ref($object, 'Bio::EnsEMBL::MiscFeature');
    $object_id = $object->dbID();
  }

  my $results = $self->fetch_all_by_Object($object_id, 'misc', $code, 'misc_feature');

  return $results;

}

sub fetch_all_by_Slice {
  my ($self, $object, $code) = @_;

  my $object_id;

  if (defined($object)) {
    assert_ref($object, 'Bio::EnsEMBL::Slice');
    $object_id = $object->get_seq_region_id();
  }

  my $results = $self->fetch_all_by_Object($object_id, 'seq_region', $code);

  return $results;

}

sub fetch_all_by_Gene {
  my ($self, $object, $code) = @_;

  my $object_id;

  if (defined($object)) {
    assert_ref($object, 'Bio::EnsEMBL::Gene');
    $object_id = $object->dbID();
  }

  my $results = $self->fetch_all_by_Object($object_id, 'gene', $code);

  return $results;

}

sub fetch_all_by_Transcript {
  my ($self, $object, $code) = @_;

  my $object_id;

  if (defined($object)) {
    assert_ref($object, 'Bio::EnsEMBL::Transcript');
    $object_id = $object->dbID();
  }

  my $results = $self->fetch_all_by_Object($object_id, 'transcript', $code);

  return $results;

}

sub fetch_all_by_Translation {
  my ($self, $object, $code) = @_;

  my $object_id;

  if (defined($object)) {
    assert_ref($object, 'Bio::EnsEMBL::Translation');
    $object_id = $object->dbID();
  }

  my $results = $self->fetch_all_by_Object($object_id, 'translation', $code);
  
  return $results;

}

sub fetch_all_by_DnaDnaAlignFeature {
  my ($self, $object, $code) = @_;

  my $object_id;

  if (defined($object)) {
    assert_ref($object, 'Bio::EnsEMBL::DnaDnaAlignFeature');
    $object_id = $object->dbID();
  }

  my $results = $self->fetch_all_by_Object($object_id, 'dna_align_feature', $code);
  
  return $results;

}




#
# _id_check
#
# backwards compatibility check:
# check if $ensID is an object; if so, return $obj->dbID
#

sub _id_check {
  my $self  = shift;
  my $ensID = shift;

  if ($ensID =~ /^\d+$/) {
	return $ensID;

  } elsif (   ref($ensID) eq 'Bio::EnsEMBL::Gene'
		   or ref($ensID) eq 'Bio::EnsEMBL::Transcript'
		   or ref($ensID) eq 'Bio::EnsEMBL::Translation')
  {

	warning("You should pass a dbID rather than an ensembl object to store the attribute on");

	if ($ensID->dbID) {
	  return $ensID->dbID;
	} else {
	  throw("Ensembl object " . $ensID->display_id . " doesn't have a dbID, can't store attribute");
	}

  } else {
	throw("Invalid dbID");
  }

} ## end sub _id_check

# _store_type

sub _store_type {
  my $self   = shift;
  my $attrib = shift;

  my $insert_ignore = $self->insert_ignore_clause();
  my $sth1 = $self->prepare
    ("${insert_ignore} INTO attrib_type (code, name, description) values (?, ?, ?)" );


  $sth1->bind_param(1, $attrib->code,        SQL_VARCHAR);
  $sth1->bind_param(2, $attrib->name,        SQL_VARCHAR);
  $sth1->bind_param(3, $attrib->description, SQL_LONGVARCHAR);

  my $rows_inserted = $sth1->execute();

  my $atid = $self->last_insert_id('attrib_type_id', undef, 'attrib_type');

  if ($rows_inserted == 0) {
	# the insert failed because the code is already stored
	my $sth2 = $self->prepare("SELECT attrib_type_id FROM attrib_type " . "WHERE code = ?");
	$sth2->bind_param(1, $attrib->code, SQL_VARCHAR);
	$sth2->execute();
	($atid) = $sth2->fetchrow_array();

	$sth2->finish();

	if (!$atid) {
	  throw("Could not store or fetch attrib_type code [" . $attrib->code . "]\n" . "Wrong database user/permissions?");
	}
  }

  $sth1->finish();

  return $atid;
} ## end sub _store_type

sub _obj_from_sth {
  my $self = shift;
  my $sth  = shift;

  my ($code, $name, $desc, $value);
  $sth->bind_columns(\$code, \$name, \$desc, \$value);

  my @results;
  while ($sth->fetch()) {
	push @results,
	  Bio::EnsEMBL::Attribute->new_fast({'code'        => $code,
										 'name'        => $name,
										 'description' => $desc,
										 'value'       => $value});
  }

  return \@results;
}

1;

