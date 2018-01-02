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

=head1 NAME

Bio::EnsEMBL::DBSQL::GenomeContainer - Encapsulates all access to 
genome related information

=head1 SYNOPSIS

  use Bio::EnsEMBL::Registry;

  Bio::EnsEMBL::Registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
  );

  $genome =
    Bio::EnsEMBL::Registry->get_adaptor( "human", "core", "GenomeContainer" );

  my $version = $genome->get_version;

  my $ref_length = $genome->get_ref_length;

  my $coord_systems = $genome->get_coord_systems;



=head1 DESCRIPTION

This module is responsible for fetching and storing genome-wide information.
Genome is an abstract object which contains information linking the species, the assembly and the ensembl annotation.

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::GenomeContainer;

use strict;
use warnings;

use Bio::EnsEMBL::Genome;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw( deprecate throw warning );
use Bio::EnsEMBL::Utils::Scalar qw( assert_ref );

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);



=head2 new

  Arg [...]  : Superclass args.  See Bio::EnsEMBL::DBSQL::BaseAdaptor
  Description: Instantiates a Bio::EnsEMBL::DBSQL::GenomeContainer
  Returntype : Bio::EnsEMBL::GenomeContainer
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


=head2 store

  Arg [1]    : Statistic
               The type of statistic to store
  Arg [2]    : Value
               The corresponding value for the statistic
  Arg [3]    : (optional) Attribute
               If more than one value exists for the statistics, it will be distinguished by its attribute
  Example    : $genome_adaptor->store('coding_cnt', 20769);
  Description: Stores a genome statistic in the database
  Returntype : The database identifier (dbID) of the newly stored genome statistic
  Exceptions : 
  Caller     : general
  Status     : Stable

=cut

sub store {
  my ($self, $statistic, $value, $attribute) = @_;

  my $stats_id = $self->fetch_by_statistic($statistic, $attribute)->dbID;

  if (defined $stats_id) {
    $self->update($statistic, $value, $attribute);
  } else {
    my $db = $self->db();
    my $species_id = $db->species_id();

    my $store_genome_sql = q{
    INSERT INTO genome_statistics 
       SET statistic = ?,
               value = ?,
          species_id = ?,
           timestamp = now()
    };

    if (defined $attribute) {
      $store_genome_sql .= ", attrib_type_id = ?";
    }

    my $sth = $self->prepare($store_genome_sql);
    $sth->bind_param(1,   $statistic,         SQL_VARCHAR);
    $sth->bind_param(2,   $value,             SQL_INTEGER);
    $sth->bind_param(3,   $species_id,        SQL_INTEGER);

    if (defined $attribute) {
      my $attribute_adaptor = $db->get_AttributeAdaptor();
      my $attribute_object = Bio::EnsEMBL::Attribute->new(-code => $attribute);
      my $attribute_type_id = $attribute_adaptor->_store_type($attribute_object);
      $sth->bind_param(4, $attribute_type_id, SQL_VARCHAR);
    }

    $sth->execute();
    $sth->finish();

    $stats_id = $sth->{'mysql_insertid'};
  }

  return $stats_id;
  
}


=head2 update

  Arg [1]    : Statistic
               The type of statistic to update
  Arg [2]    : Value
               The corresponding value for the statistic
  Arg [3]    : (optional) Attribute
               If more than one value exists for the statistics, it will be distinguished by its attribute
  Example    : $genome_adaptor->update('coding_cnt', 20769);
  Description: Updates an existing genome statistic in the database
  Returntype : none
  Exceptions :
  Caller     : general
  Status     : Stable

=cut

sub update {
  my ($self, $statistic, $value, $attribute) = @_;

  my $db = $self->db();

  my $update_genome_sql = q{
    UPDATE genome_statistics
       SET value = ?,
       timestamp = now() 
  };

  if (defined $attribute) {
    $update_genome_sql .= ', attrib_type_id = ?';
  }

  $update_genome_sql .= ' WHERE statistic = ? and species_id = ?';
  
  my $sth = $self->prepare($update_genome_sql);
  $sth->bind_param(1,   $value,     SQL_INTEGER);

  my $increment = 2;
  if (defined $attribute) {
    my $attribute_adaptor = $db->get_AttributeAdaptor();
    my $attribute_object = Bio::EnsEMBL::Attribute->new(-code => $attribute);
    my $attribute_type_id = $attribute_adaptor->_store_type($attribute_object);
    $sth->bind_param($increment, $attribute_type_id, SQL_VARCHAR);
    $increment++;
  }
  
  $sth->bind_param($increment++, $statistic, SQL_VARCHAR);
  $sth->bind_param($increment, $db->species_id(), SQL_INTEGER);

  $sth->execute();
  $sth->finish();
}


=head2 _meta_container

  Arg [1]    : none
  Example    : $meta_container = $genome->_meta_container();
  Description: Internal method to return a MetaContainer object for the genome
  Returntype : Bio::EnsEMBL::DBSQL::MetaContainer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub _meta_container {
  my $self = shift;
  return $self->db->get_adaptor('MetaContainer');
}


=head2 get_version

  Arg [1]    : (optional) assembly version
  Example    : $version = $genome->get_version();
  Description: Getter/Setter for the assembly version

  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_version {
  my ($self, $version) = @_;
  if (defined $version) {
    $self->{'version'} = $version;
  }
  if (!defined $self->{'version'}) {
    my $csa = $self->db()->get_adaptor('CoordSystem');
    $self->{'version'} = $csa->get_default_version;
  }
  return $self->{'version'};
}

=head2 get_accession

  Arg [1]    : (optional) assembly accession
  Example    : $accession = $genome->get_accession();
  Description: Getter/setter for the accession of the assembly currently used

  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_accession {
  my ($self, $accession) = @_;
  if (defined $accession) {
    $self->{'accession'} = $accession;
  }
  if (!defined $self->{'accession'}) {
    $self->{'accession'} = $self->_meta_container->single_value_by_key('assembly.accession');
  }
  return $self->{'accession'};
}


=head2 get_assembly_name

  Arg [1]    : (optional) assembly name
  Example    : $assembly_name = $genome->get_assembly_name();
  Description: Getter/setter for the name of the assembly currently used

  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_assembly_name {
  my ($self, $assembly_name) = @_;
  if (defined $assembly_name) {
    $self->{'assembly_name'} = $assembly_name;
  }
  if (!defined $self->{'assembly_name'}) {
    $self->{'assembly_name'} = $self->_meta_container->single_value_by_key('assembly.name');
  }
  return $self->{'assembly_name'};
}


=head2 get_assembly_date

  Arg [1]    : (optional) assembly date
  Example    : $assembly_date = $genome->get_assembly_date();
  Description: Getter/setter for the date of the assembly currently used

  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_assembly_date {
  my ($self, $assembly_date) = @_;
  if (defined $assembly_date) {
    $self->{'assembly_date'} = $assembly_date;
  }
  if (!defined $self->{'assembly_date'}) {
    $self->{'assembly_date'} = $self->_meta_container->single_value_by_key('assembly.date');
  }
  return $self->{'assembly_date'};
}


=head2 get_genebuild_start_date

  Arg [1]    : (optional) genebuild start date
  Example    : $genebuild_start_date = $genome->get_genebuild_start_date();
  Description: Getter/setter for the start date of the genebuild currently used

  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_genebuild_start_date {
  my ($self, $genebuild_start_date) = @_;
  if (defined $genebuild_start_date) {
    $self->{'genebuild_start_date'} = $genebuild_start_date;
  }
  if (!defined $self->{'genebuild_start_date'}) {
    $self->{'genebuild_start_date'} = $self->_meta_container->single_value_by_key('genebuild.start_date');
  }
  return $self->{'genebuild_start_date'};
}


=head2 get_genebuild_method

  Arg [1]    : (optional) genebuild start date
  Example    : $genebuild_method = $genome->get_genebuild_method();
  Description: Getter/setter for the method of the genebuild currently used

  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_genebuild_method {
  my ($self, $genebuild_method) = @_;
  if (defined $genebuild_method) {
    $self->{'genebuild_method'} = $genebuild_method;
  }
  if (!defined $self->{'genebuild_method'}) {
    $self->{'genebuild_method'} = $self->_meta_container->single_value_by_key('genebuild.method');
  }
  return $self->{'genebuild_method'};
}


=head2 get_genebuild_initial_release_date

  Arg [1]    : (optional) genebuild initial release date
  Example    : $genebuild_initial_release_date = $genome->get_initial_release_date();
  Description: Getter/setter for the initial release date of the genebuild currently used

  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_genebuild_initial_release_date {
  my ($self, $genebuild_initial_release_date) = @_;
  if (defined $genebuild_initial_release_date) {
    $self->{'genebuild_initial_release_date'} = $genebuild_initial_release_date;
  }
  if (!defined $self->{'genebuild_initial_release_date'}) {
    $self->{'genebuild_initial_release_date'} = $self->_meta_container->single_value_by_key('genebuild.initial_release_date');
  }
  return $self->{'genebuild_initial_release_date'};
}


=head2 get_genebuild_last_geneset_update

  Arg [1]    : (optional) genebuild last geneset update
  Example    : $genebuild_last_geneset_update = $genome->get_last_geneset_update();
  Description: Getter/setter for the last geneset update of the genebuild currently used

  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_genebuild_last_geneset_update {
  my ($self, $genebuild_last_geneset_update) = @_;
  if (defined $genebuild_last_geneset_update) {
    $self->{'genebuild_last_geneset_update'} = $genebuild_last_geneset_update;
  }
  if (!defined $self->{'genebuild_last_geneset_update'}) {
    $self->{'genebuild_last_geneset_update'} = $self->_meta_container->single_value_by_key('genebuild.last_geneset_update');
  }
  return $self->{'genebuild_last_geneset_update'};
}


=head2 _get_length

  Arg [1]    : none
  Example    : $length = $genome->_get_length('toplevel');
  Description: Internal method to return the length for a type of slices
  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub _get_length {
  my ($self, $cs_name) = @_;
  my $slice_adaptor = $self->db->get_adaptor('Slice');
  my $seqlevel = $slice_adaptor->fetch_all($cs_name);
  my $count;
  foreach my $seq (@$seqlevel) {
    $count += $seq->length();
  }
  return $count;
}



=head2 get_ref_length

  Arg [1]    : (optional) golden path length
  Example    : $ref_length = $genome->get_ref_length();
  Description: Getter/setter for the golden path of the assembly currently used
  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_ref_length {
  my ($self, $ref_length) = @_;
  if (defined $ref_length) {
    $self->{'ref_length'} = $ref_length;
  }
  if (!defined $self->{'ref_length'}) {
    $self->{'ref_length'} = $self->fetch_by_statistic('ref_length')->value();
  }
  return $self->{'ref_length'};
}


=head2 get_total_length

  Arg [1]    : (optional) base pair length
  Example    : $total_length = $genome->get_total_length();
  Description: Getter/setter for the total length (number of base pairs) for the assembly currently used

  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_total_length {
  my ($self, $total_length) = @_;
  if (defined $total_length) {
    $self->{'total_length'} = $total_length;
  }
  if (!defined $self->{'total_length'}) {
    $self->{'total_length'} = $self->fetch_by_statistic('total_length')->value;
  }
  return $self->{'total_length'};
}

=head2 get_toplevel

  Arg [1]    : none
  Example    : $toplevel = $genome->get_toplevel();
  Description: Returns the toplevel for the assembly currently used

  Returntype : ListRef of Bio::EnsEMBL::Slice
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_toplevel {
  my ($self) = @_;
  my $sa = $self->db->get_adaptor('Slice');
  $self->{'toplevel'} = $sa->fetch_all('toplevel', undef, undef, 1);
  return $self->{'toplevel'};
}


=head2 get_karyotype

  Arg [1]    : none
  Example    : $karyotype = $genome->get_karyotype();
  Description: Returns the karyotype for the assembly currently used

  Returntype : ListRef of Bio::EnsEMBL::Slice
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_karyotype {
  my ($self) = @_;
  my $sa = $self->db->get_adaptor('Slice');
  $self->{'karyotype'} = $sa->fetch_all_karyotype;
  return $self->{'karyotype'};
}

=head2 get_coord_systems

  Arg [1]    : none
  Example    : $coord_systems = $genome->get_coord_systems();
  Description: Returns the coord_systems for the assembly currently used

  Returntype : ListRef of Bio::EnsEMBL::CoordSystem
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_coord_systems {
  my ($self, $all) = @_;
  my $csa = $self->db->get_adaptor('CoordSystem');
  if (!$all) {
    my $version = $self->get_version();
    $self->{'coord_systems'} = $csa->fetch_all_by_version($version);
  } else {
    $self->{'coord_systems'} = $csa->fetch_all();
  }
  return $self->{'coord_systems'};
}

=head2 _get_count

  Arg [1]    : none
  Example    : $count = $genome->_get_count('coding_cnt');
  Description: Internal method to return a count for a given attribute code
  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub _get_count {
  my ($self, $code, $attribute) = @_;
  my $statistic = $self->fetch_by_statistic($code, $attribute);
  return $statistic->value();
}

=head2 get_count

  Arg [1]    : none
  Example    : $count = $genome->get_count('coding_cnt');
  Description: Retrieve a count for a given attribute code
  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_count {
  my ($self, $code, $attribute) = @_;
  my $statistic = $self->fetch_by_statistic($code, $attribute);
  return $statistic->value();
}

=head2 fetch_all_statistics

  Arg [1]    : none
  Example    : $list = $genome->fetch_all_statistics();
  Description: Retrieve all entries stored in the genome_statistics table
  Returntype : ArrayRef of Bio::EnsEMBL::Genome
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_statistics {
  my ($self) = @_;
  my $db = $self->db;
  my $species_id = $self->db->species_id();
  my @results;
  my $sql = q{
    SELECT genome_statistics_id, statistic, value, species_id, code, name, description
      FROM genome_statistics, attrib_type
     WHERE genome_statistics.attrib_type_id = attrib_type.attrib_type_id
       AND species_id=?
  };

  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $species_id, SQL_INTEGER);
  $sth->execute();
  my $results = $self->_obj_from_sth($sth);
  $sth->finish();

  return $results;
}


=head2 fetch_by_statistic

  Arg [1]    : string $statistic
  Example    : $results = $genome->fetch_by_statistic('coding_cnt');
  Description: Returns a Genome object for a given statistic
  Returntype : Bio::EnsEMBL::Genome
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_statistic {
  my ($self, $statistic_name, $attribute) = @_;
  my $db = $self->db;
  my $fetch_sql = q{
    SELECT genome_statistics_id, statistic, value, species_id, code, name, description
      FROM genome_statistics, attrib_type 
     WHERE genome_statistics.attrib_type_id = attrib_type.attrib_type_id
       AND statistic = ? AND species_id=?
  };
  if (defined $attribute) {
    $fetch_sql .= " AND code = ?";
  }

  my $sth = $self->prepare($fetch_sql);
  $sth->bind_param(1, $statistic_name, SQL_VARCHAR);
  $sth->bind_param(2, $self->db->species_id, SQL_INTEGER);
  if (defined $attribute) {
    $sth->bind_param(3, $attribute, SQL_VARCHAR);
  }
  $sth->execute();
  my ($dbID, $statistic, $value, $species_id, $code, $name, $desc);
  $sth->bind_columns(\$dbID, \$statistic, \$value, \$species_id, \$code, \$name, \$desc);

  my @results = $sth->fetchrow_array;
  $sth->finish();

  return Bio::EnsEMBL::Genome->new_fast({'dbID'        => $dbID,
                                         'statistic'   => $statistic,
                                         'code'        => $code,
                                         'name'        => $name,
                                         'description' => $desc,
                                         'value'       => $value});

}

=head2 is_empty

  Arg [1]    : none
  Example    : $results = $genome->is_empty;
  Description: Boolean to check if there is data in the genome container
  Returntype : Boolean
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub is_empty {
  my $self = shift;
  my $db = $self->db;
  my $species_id = $self->db->species_id();
  my $is_empty = 1;
  my $count_sql = q{
    SELECT count(*) FROM genome_statistics
  };

  my $sth = $self->prepare($count_sql);
  $sth->execute();
  if ($sth->fetchrow()) {
    $is_empty = 0;
  }
  $sth->finish();
  return $is_empty;
}

=head2 get_attrib

  Arg [1]    : statistic
  Example    : $results = $genome->_get_attrib('coding_cnt');
  Description: Returns the attribute object for a given statistic
  Returntype : Bio::EnsEMBL::Attrib
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_attrib {
  my ($self, $statistic) = @_;
  my $db = $self->db();
  my $attribute_adaptor = $db->get_adaptor('attribute');
  my @attribs = @{ $attribute_adaptor->fetch_by_code($statistic) };
  my $attrib = Bio::EnsEMBL::Attribute->new(
  -code => $attribs[1],
  -name => $attribs[2],
  -description => $attribs[3]
  );
  return $attrib;
}

=head2 get_coding_count

  Arg [1]    : (optional) coding count
  Example    : $coding_count = $genome->get_coding_count();
  Description: Getter/setter for the number of coding genes in the current build

  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_coding_count {
  my ($self, $coding_count) = @_;
  if (defined $coding_count) {
    $self->{'coding_count'} = $coding_count;
  }
  if (!defined $self->{'coding_count'}) {
    $self->{'coding_count'} = $self->_get_count('coding_cnt');
  }
  return $self->{'coding_count'};
}

=head2 get_rcoding_count

  Arg [1]    : (optional) readthrough coding count
  Example    : $rcoding_count = $genome->get_rcoding_count();
  Description: Getter/setter for the number of readthrough coding genes in the current build

  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_rcoding_count {
  my ($self, $rcoding_count) = @_;
  if (defined $rcoding_count) {
    $self->{'rcoding_count'} = $rcoding_count;
  }
  if (!defined $self->{'rcoding_count'}) {
    $self->{'rcoding_count'} = $self->_get_count('coding_rcnt');
  }
  return $self->{'rcoding_count'};
}


=head2 get_snoncoding_count

  Arg [1]    : (optional) short non coding count
  Example    : $snoncoding_count = $genome->get_snoncoding_count();
  Description: Getter/setter for the number of short non coding genes in the current build

  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_snoncoding_count {
  my ($self, $snoncoding_count) = @_;
  if (defined $snoncoding_count) {
    $self->{'snoncoding_count'} = $snoncoding_count;
  }
  if (!defined $self->{'snoncoding_count'}) {
    $self->{'snoncoding_count'} = $self->_get_count('noncoding_cnt_s');
  }
  return $self->{'snoncoding_count'};
}

=head2 get_rsnoncoding_count

  Arg [1]    : (optional) readthrough short non coding count
  Example    : $rsnoncoding_count = $genome->get_rsnoncoding_count();
  Description: Getter/setter for the number of readthrough short non coding genes in the current build

  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_rsnoncoding_count {
  my ($self, $rsnoncoding_count) = @_;
  if (defined $rsnoncoding_count) {
    $self->{'rsnoncoding_count'} = $rsnoncoding_count;
  }
  if (!defined $self->{'rsnoncoding_count'}) {
    $self->{'rsnoncoding_count'} = $self->_get_count('noncoding_rcnt_s');
  }
  return $self->{'rsnoncoding_count'};
}

=head2 get_mnoncoding_count

  Arg [1]    : (optional) miscellaneous non coding count
  Example    : $mnoncoding_count = $genome->get_mnoncoding_count();
  Description: Getter/setter for the number of miscellaneous non coding genes in the current build

  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_mnoncoding_count {
  my ($self, $mnoncoding_count) = @_;
  if (defined $mnoncoding_count) {
    $self->{'mnoncoding_count'} = $mnoncoding_count;
  }
  if (!defined $self->{'mnoncoding_count'}) {
    $self->{'mnoncoding_count'} = $self->_get_count('noncoding_cnt_m');
  }
  return $self->{'mnoncoding_count'};
}

=head2 get_rmnoncoding_count

  Arg [1]    : (optional) readthrough miscellaneous non coding count
  Example    : $rmnoncoding_count = $genome->get_rmnoncoding_count();
  Description: Getter/setter for the number of readthrough miscellaneous non coding genes in the current build

  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_rmnoncoding_count {
  my ($self, $rmnoncoding_count) = @_;
  if (defined $rmnoncoding_count) {
    $self->{'rmnoncoding_count'} = $rmnoncoding_count;
  }
  if (!defined $self->{'rmnoncoding_count'}) {
    $self->{'rmnoncoding_count'} = $self->_get_count('noncoding_rcnt_m');
  }
  return $self->{'rmnoncoding_count'};
}


=head2 get_lnoncoding_count

  Arg [1]    : (optional) long non coding count
  Example    : $lnoncoding_count = $genome->get_lnoncoding_count();
  Description: Getter/setter for the number of long non coding genes in the current build

  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_lnoncoding_count {
  my ($self, $lnoncoding_count) = @_;
  if (defined $lnoncoding_count) {
    $self->{'lnoncoding_count'} = $lnoncoding_count;
  }
  if (!defined $self->{'lnoncoding_count'}) {
    $self->{'lnoncoding_count'} = $self->_get_count('noncoding_cnt_l');
  }
  return $self->{'lnoncoding_count'};
}

=head2 get_rlnoncoding_count

  Arg [1]    : (optional) readthrough long non coding count
  Example    : $rlnoncoding_count = $genome->get_rlnoncoding_count();
  Description: Getter/setter for the number of readthrough long non coding genes in the current build

  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_rlnoncoding_count {
  my ($self, $rlnoncoding_count) = @_;
  if (defined $rlnoncoding_count) {
    $self->{'rlnoncoding_count'} = $rlnoncoding_count;
  }
  if (!defined $self->{'rlnoncoding_count'}) {
    $self->{'rlnoncoding_count'} = $self->_get_count('noncoding_rcnt_l');
  }
  return $self->{'rlnoncoding_count'};
}

=head2 get_pseudogene_count

  Arg [1]    : (optional) pseudogene count
  Example    : $pseudogene_count = $genome->get_pseudogene_count();
  Description: Getter/setter for the number of pseudogenes in the current build

  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut


sub get_pseudogene_count {
  my ($self, $pseudogene_count) = @_;
  if (defined $pseudogene_count) {
    $self->{'pseudogene_count'} = $pseudogene_count;
  }
  if (!defined $self->{'pseudogene_count'}) {
    $self->{'pseudogene_count'} = $self->_get_count('pseudogene_cnt');
  }
  return $self->{'pseudogene_count'};
}

=head2 get_rpseudogene_count

  Arg [1]    : (optional) readthrough pseudogene count
  Example    : $rpseudogene_count = $genome->get_rpseudogene_count();
  Description: Getter/setter for the number of readthrough pseudogenes in the current build

  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut


sub get_rpseudogene_count {
  my ($self, $rpseudogene_count) = @_;
  if (defined $rpseudogene_count) {
    $self->{'rpseudogene_count'} = $rpseudogene_count;
  }
  if (!defined $self->{'rpseudogene_count'}) {
    $self->{'rpseudogene_count'} = $self->_get_count('pseudogene_rcnt');
  }
  return $self->{'rpseudogene_count'};
}

=head2 get_alt_coding_count

  Arg [1]    : (optional) alt coding count
  Example    : $alt_coding_count = $genome->get_alt_coding_count();
  Description: Getter/setter for the number of coding genes on alternate sequences

  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_alt_coding_count {
  my ($self, $alt_coding_count) = @_;
  if (defined $alt_coding_count) {
    $self->{'alt_coding_count'} = $alt_coding_count;
  }
  if (!defined $self->{'alt_coding_count'}) {
    $self->{'alt_coding_count'} = $self->_get_count('coding_acnt');
  }
  return $self->{'alt_coding_count'};
}

=head2 get_alt_rcoding_count

  Arg [1]    : (optional) alt readthrough coding count
  Example    : $alt_rcoding_count = $genome->get_alt_rcoding_count();
  Description: Getter/setter for the number of readthrough coding genes on alternate sequences

  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_alt_rcoding_count {
  my ($self, $alt_rcoding_count) = @_;
  if (defined $alt_rcoding_count) {
    $self->{'alt_rcoding_count'} = $alt_rcoding_count;
  }
  if (!defined $self->{'alt_rcoding_count'}) {
    $self->{'alt_rcoding_count'} = $self->_get_count('coding_racnt');
  }
  return $self->{'alt_rcoding_count'};
}


=head2 get_alt_snoncoding_count

  Arg [1]    : (optional) alt short non coding count
  Example    : $alt_snoncoding_count = $genome->get_alt_snoncoding_count();
  Description: Getter/setter for the number of short non coding genes on alternate sequences

  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_alt_snoncoding_count {
  my ($self, $alt_snoncoding_count) = @_;
  if (defined $alt_snoncoding_count) {
    $self->{'alt_snoncoding_count'} = $alt_snoncoding_count;
  }
  if (!defined $self->{'alt_snoncoding_count'}) {
    $self->{'alt_snoncoding_count'} = $self->_get_count('noncoding_acnt_s');
  }
  return $self->{'alt_snoncoding_count'};
}

=head2 get_alt_rsnoncoding_count

  Arg [1]    : (optional) alt readthrough short non coding count
  Example    : $alt_rsnoncoding_count = $genome->get_alt_rsnoncoding_count();
  Description: Getter/setter for the number of readthrough short non coding genes on alternate sequences

  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_alt_rsnoncoding_count {
  my ($self, $alt_rsnoncoding_count) = @_;
  if (defined $alt_rsnoncoding_count) {
    $self->{'alt_rsnoncoding_count'} = $alt_rsnoncoding_count;
  }
  if (!defined $self->{'alt_rsnoncoding_count'}) {
    $self->{'alt_rsnoncoding_count'} = $self->_get_count('noncoding_racnt_s');
  }
  return $self->{'alt_rsnoncoding_count'};
}

=head2 get_alt_mnoncoding_count

  Arg [1]    : (optional) alt miscellaneous non coding count
  Example    : $alt_mnoncoding_count = $genome->get_alt_mnoncoding_count();
  Description: Getter/setter for the number of miscellaneous non coding genes on alternate sequences

  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_alt_mnoncoding_count {
  my ($self, $alt_mnoncoding_count) = @_;
  if (defined $alt_mnoncoding_count) {
    $self->{'alt_mnoncoding_count'} = $alt_mnoncoding_count;
  }
  if (!defined $self->{'alt_mnoncoding_count'}) {
    $self->{'alt_mnoncoding_count'} = $self->_get_count('noncoding_acnt_m');
  }
  return $self->{'alt_mnoncoding_count'};
}

=head2 get_alt_rmnoncoding_count

  Arg [1]    : (optional) alt readthrough miscellaneous non coding count
  Example    : $alt_rmnoncoding_count = $genome->get_alt_rmnoncoding_count();
  Description: Getter/setter for the number of readthrough miscellaneous non coding genes on alternate sequences

  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_alt_rmnoncoding_count {
  my ($self, $alt_rmnoncoding_count) = @_;
  if (defined $alt_rmnoncoding_count) {
    $self->{'alt_rmnoncoding_count'} = $alt_rmnoncoding_count;
  }
  if (!defined $self->{'alt_rmnoncoding_count'}) {
    $self->{'alt_rmnoncoding_count'} = $self->_get_count('noncoding_racnt_m');
  }
  return $self->{'alt_rmnoncoding_count'};
}


=head2 get_alt_lnoncoding_count

  Arg [1]    : (optional) alt long non coding count
  Example    : $alt_lnoncoding_count = $genome->get_alt_lnoncoding_count();
  Description: Getter/setter for the number of long non coding genes on alternate sequences

  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_alt_lnoncoding_count {
  my ($self, $alt_lnoncoding_count) = @_;
  if (defined $alt_lnoncoding_count) {
    $self->{'alt_lnoncoding_count'} = $alt_lnoncoding_count;
  }
  if (!defined $self->{'alt_lnoncoding_count'}) {
    $self->{'alt_lnoncoding_count'} = $self->_get_count('noncoding_acnt_l');
  }
  return $self->{'alt_lnoncoding_count'};
}

=head2 get_alt_rlnoncoding_count

  Arg [1]    : (optional) alt readthrough long non coding count
  Example    : $alt_lnoncoding_count = $genome->get_alt_lnoncoding_count();
  Description: Getter/setter for the number of readthrough long non coding genes on alternate sequences

  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_alt_rlnoncoding_count {
  my ($self, $alt_rlnoncoding_count) = @_;
  if (defined $alt_rlnoncoding_count) {
    $self->{'alt_rlnoncoding_count'} = $alt_rlnoncoding_count;
  }
  if (!defined $self->{'alt_rlnoncoding_count'}) {
    $self->{'alt_rlnoncoding_count'} = $self->_get_count('noncoding_racnt_l');
  }
  return $self->{'alt_rlnoncoding_count'};
}



=head2 get_alt_pseudogene_count

  Arg [1]    : (optional) alt pseudogene count
  Example    : $alt_pseudogene_count = $genome->get_alt_pseudogene_count();
  Description: Getter/setter for the number of pseudogenes on alternate sequences

  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_alt_pseudogene_count {
  my ($self, $alt_pseudogene_count) = @_;
  if (defined $alt_pseudogene_count) {
    $self->{'alt_pseudogene_count'} = $alt_pseudogene_count;
  }
  if (!defined $self->{'alt_pseudogene_count'}) {
    $self->{'alt_pseudogene_count'} = $self->_get_count('pseudogene_acnt');
  }
  return $self->{'alt_pseudogene_count'};
}

=head2 get_alt_rpseudogene_count

  Arg [1]    : (optional) alt readthrough pseudogene count
  Example    : $alt_rpseudogene_count = $genome->get_alt_pseudogene_count();
  Description: Getter/setter for the number of readthrough pseudogenes on alternate sequences

  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_alt_rpseudogene_count {
  my ($self, $alt_rpseudogene_count) = @_;
  if (defined $alt_rpseudogene_count) {
    $self->{'alt_rpseudogene_count'} = $alt_rpseudogene_count;
  }
  if (!defined $self->{'alt_rpseudogene_count'}) {
    $self->{'alt_rpseudogene_count'} = $self->_get_count('pseudogene_racnt');
  }
  return $self->{'alt_rpseudogene_count'};
}

=head2 get_short_variation_count

  Arg [1]    : (optional) short variation count
  Example    : $short_variation_count = $genome->get_short_variation_count();
  Description: Getter/setter for the number of short variants in the current build

  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_short_variation_count {
  my ($self, $short_variation_count) = @_;
  if (defined $short_variation_count) {
    $self->{'short_variation_count'} = $short_variation_count;
  }
  if (!defined $self->{'short_variation_count'}) {
    $self->{'short_variation_count'} = $self->_get_count('SNPCount');
  }
  return $self->{'short_variation_count'};
}


=head2 get_prediction_count

  Arg [1]    : (optional) logic_name
  Example    : $prediction_count = $genome->get_prediction_count();
  Description: Getter for the number of predicted genes in the current build
               Can be restricted to a given analysis

  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_prediction_count {
  my ($self, $logic_name) = @_;
  return $self->_get_count('PredictionTranscript', $logic_name);
}


=head2 get_structural_variation_count

  Arg [1]    : none
  Example    : $structural_variation_count = $genome->get_structural_variation_count();
  Description: Return the number of structural variations in the current build
  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_structural_variation_count {
  my ($self, $structural_variation_count) = @_;
  if (defined $structural_variation_count) {
    $self->{'structural_variation_count'} = $structural_variation_count;
  }
  if (!defined $self->{'structural_variation_count'}) {
    $self->{'structural_variation_count'} = $self->_get_count('StructuralVariation');
   }
  return $self->{'structural_variation_count'};
}

=head2 get_transcript_count

  Arg [1]    : (optional) transcript count
  Example    : $transcript_count = $genome->get_transcript_count();
  Description: Getter/setter for the number of transcripts in the current build
  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_transcript_count {
  my ($self, $transcript_count) = @_;
  if (defined $transcript_count) {
    $self->{'transcript_count'} = $transcript_count;
  }
  if (!defined $self->{'transcript_count'}) {
    $self->{'transcript_count'} = $self->_get_count('transcript');
  }
  return $self->{'transcript_count'};
}

=head2 get_alt_transcript_count

  Arg [1]    : (optional) alt transcript count
  Example    : $alt_transcript_count = $genome->get_alt_transcript_count();
  Description: Getter/setter for the number of transcripts on alternate sequences in the current build
  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_alt_transcript_count {
  my ($self, $alt_transcript_count) = @_;
  if (defined $alt_transcript_count) {
    $self->{'alt_transcript_count'} = $alt_transcript_count;
  }
  if (!defined $self->{'alt_transcript_count'}) {
    $self->{'alt_transcript_count'} = $self->_get_count('alt_transcript');
  }
  return $self->{'alt_transcript_count'};
}


=head2 has_karyotype

  Arg        : None
  Example    : $has_karyotype = $genome->has_karyotype();
  Description: Boolean indicating whether a genome has a karyotype (ie chromosomes)
               or not
  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub has_karyotype {
  my $self = shift;

  my $db = $self->db();
  my $slice_adaptor = $db->get_SliceAdaptor();
  my $karyotype = $slice_adaptor->fetch_all_karyotype;

  return 0 unless scalar(@$karyotype);

  return 1;
}


=head2 is_high_coverage

  Arg        : None
  Example    : $is_high_coverage = $genome->is_high_coverage();
  Description: Boolean indicating whether an assembly is high coverage
               or not
  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub is_high_coverage {
  my $self = shift;

  my $coverage_depth = $self->_meta_container->single_value_by_key('assembly.coverage_depth');

  return 0 if !$coverage_depth;
  $coverage_depth = lc($coverage_depth);

  if ($coverage_depth eq 'high') {
    return 1;
  } elsif (($coverage_depth eq 'low') or ($coverage_depth eq 'medium')) {
    return 0;
  } elsif ($coverage_depth =~ /^([0-9]+)x$/) {
    return $1<6 ? 0 : 1;
  }

  return 0;
}

=head2 is_polyploid

  Arg        : None
  Example    : $is_polyploid = $genome->is_polyploid();
  Description: Returns whether the genome is or is not polyploid.
  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub is_polyploid {
  my $self = shift;

  my $polyploid = $self->_meta_container->single_value_by_key('ploidy');

  # polyploid could be not defined, meta_key is optional
  return 0 unless defined $polyploid;

  return $polyploid > 2;
}

=head2 get_genome_components

  Arg        : None
  Example    : $components = $genome->get_genome_components();
  Description: Returns the list of (diploid) components, for a
               polyploid genome
  Returntype : Arrayref
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_genome_components {
  my $self = shift;

  my $sql_helper = $self->dbc->sql_helper;

  my $sql = 
    "SELECT DISTINCT value 
     FROM seq_region_attrib JOIN attrib_type 
     USING (attrib_type_id) WHERE attrib_type.code='genome_component'";

  return $sql_helper->execute_simple(-SQL => $sql);
}

sub _obj_from_sth {
  my $self = shift;
  my $sth  = shift;

  my ($dbID, $statistic, $value, $species_id, $code, $name, $desc);
  $sth->bind_columns(\$dbID, \$statistic, \$value, \$species_id, \$code, \$name, \$desc);

  my @results;
  while ($sth->fetch()) {
        push @results,
          Bio::EnsEMBL::Genome->new_fast({'dbID'        => $dbID,
                                          'statistic'   => $statistic,
                                          'code'        => $code,
                                          'name'        => $name,
                                          'description' => $desc,
                                          'value'       => $value});
  }

  return \@results;
}

1;
