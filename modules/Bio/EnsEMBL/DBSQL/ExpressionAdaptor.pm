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

Bio::EnsEMBL::DBSQL::ExpressionAdaptor - Provides database interaction for
Bio::EnsEMBL::Expression objects.


=head1 SYNOPSIS

  # $db is a Bio::EnsEMBL::DBSQL::DBAdaptor object:
  $expression_adaptor = $db->get_ExpressionAdaptor();

  $attributes = $expression_adaptor->fetch_all_by_gene($tissue);

  $expression_adaptor->store_on_Gene( $gene, \@expressions);


=head1 DESCRIPTION

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::ExpressionAdaptor;

use strict;
use warnings;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Expression;

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

  if($self->can($sub)){
    return $self->$sub($type,@args);
  }
  else{
    warn("In ExpressionAdaptor cannot call sub $sub$type\n");
  }
  return undef;
}



sub store_on_ {
  my $self        = shift;
  my $type        = shift;
  my $object      = shift;
  my $expressions = shift;
  my $table;


  my $object_id = $object->dbID();
  $type = lc($type);

  my $sth = $self->prepare( "INSERT into ".$type."_expression ".
			    "SET ".$type."_id = ?, tissue_id = ?, ".
			    "value = ? " );

  for my $exp ( @$expressions) {

    if(!ref($exp) && $exp->isa('Bio::EnsEMBL::Expression')) {
      throw("Reference to list of Bio::EnsEMBL::Expression objects " .
            "argument expected.");
    }

    my $expid = $self->_store_type( $exp);
    $sth->bind_param(1,$object_id,SQL_INTEGER);
    $sth->bind_param(2,$expid,SQL_INTEGER);
    $sth->bind_param(3,$exp->value,SQL_VARCHAR);
    $sth->execute();
  }

  return;
}


sub remove_from_ {
  my $self   = shift;
  my $type   = shift;
  my $object = shift;
  my $name   = shift;
  my $table;

  if(!ref($object) || !$object->isa('Bio::EnsEMBL::'.$type)) {
    throw("$type argument is required but you passed $object");
  }

  my $object_id = $object->dbID();

  if(!defined($object_id)) {
    throw("$type must have dbID.");
  }

  $table = lc($type);
  my $sth;
  if(defined($name)){
    $sth = $self->prepare("DELETE e FROM ".$table."_expression e, tissue t " .
                         "WHERE t.tissue_id = e.tissue_id AND ".
                         "e.".$type."_id = ? AND ".
			 "t.name like ?");
    $sth->bind_param(1,$object_id,SQL_INTEGER);
    $sth->bind_param(2,$name,SQL_VARCHAR);
  }
  else{
    $sth = $self->prepare("DELETE FROM ".$table."_expression " .
                         "WHERE ".$type."_id = ?");
    $sth->bind_param(1,$object_id,SQL_INTEGER);
  }
  $sth->execute();

  $sth->finish();

  return;
}

sub fetch_all {
}


sub get_all_tissues {
  my $self  = shift;

  my $sql = "SELECT DISTINCT name FROM tissue";
  my $sth = $self->prepare($sql);
  $sth->execute();
  my ($name, @out);

  $sth->bind_columns(\$name);
  while ($sth->fetch()) {
    push @out, $name;
  }
  return \@out;
}



sub fetch_all_by_{
  my $self       = shift;
  my $type       = shift;
  my $object     = shift;
  my $name       = shift;
  my $logic_name = shift;
  my $cutoff     = shift;

  if (defined($object)){
    if(!ref($object) || !$object->isa('Bio::EnsEMBL::'.$type)) {
      throw("$type argument is required. but you passed $object");
    }
  }

  $type = lc($type);
  my $object_id;
  $object_id = $object->dbID() if defined $object;

  my $sql = "SELECT t.name, t.description, t.ontology, e.".$type."_id, e.value, e.analysis_id " .
              "FROM ".$type."_expression e, tissue t ".
                 "WHERE e.tissue_id = t.tissue_id";

  if(defined($name)){
    $sql .= ' AND t.name like "'.$name . '"';
  }
  if(defined($object_id)){
    $sql .= " AND e.".$type."_id = ".$object_id;
  }
  if (defined($logic_name)){
    my $aa = $self->db->get_AnalysisAdaptor();
    my $an = $aa->fetch_by_logic_name($logic_name);
    my $an_id = $an->dbID();
    $sql .= ' AND e.analysis_id = ' . $an_id;
  }
  if (defined($cutoff)){
    $sql .= " AND e.value > $cutoff";
  }
		   
  my $sth = $self->prepare($sql);
  $sth->execute();
  my ($desc, $ontology, $value, $analysis_id);
  $sth->bind_columns(\$name, \$desc, \$ontology, \$object_id, \$value, \$analysis_id);

  my $object_adaptor = "get_" . $type . "Adaptor";
  my $adaptor = $self->db->$object_adaptor();
  my @out;

  while ($sth->fetch()) {
    my $object = $adaptor->fetch_by_dbID($object_id);

    my $exp = Bio::EnsEMBL::Expression->new_fast
              ( {'name' => $name,
                 'description' => $desc,
                 'ontology' => $ontology,
                 'object' => $object,
                 'value' => $value} );
    push @out, $exp;
  }

  $sth->finish();

  return \@out;
  
}


sub DESTROY{
}



# _store_type

sub _store_type {
  my $self   = shift;
  my $tissue = shift;

  my $sth1 = $self->prepare
    ("INSERT IGNORE INTO tissue set name = ?, description = ?, ontology = ?" );

  $sth1->bind_param(1,$tissue->name,SQL_VARCHAR);
  $sth1->bind_param(2,$tissue->description,SQL_LONGVARCHAR);
  $sth1->bind_param(3,$tissue->ontology,SQL_VARCHAR);

  my $rows_inserted =  $sth1->execute();

  my $atid = $sth1->{'mysql_insertid'};

  if($rows_inserted == 0) {
    # the insert failed because the tissue is already stored
    my $sth2 = $self->prepare
      ("SELECT tissue_id FROM tissue " .
       "WHERE name = ?");
    $sth2->bind_param(1,$tissue->name,SQL_VARCHAR);
    $sth2->execute();
    ($atid) = $sth2->fetchrow_array();

    $sth2->finish();

    if(!$atid) {
      throw("Could not store or fetch tissue [".$tissue->name."]\n" .
	    "Wrong database user/permissions?");
    }
  }

  $sth1->finish();

  return $atid;
}



1;

