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

package XrefParser::InterproFromCoreParser;
  
use strict;
use warnings;
use Carp;

use vars qw(@ISA);
@ISA = qw(XrefParser::BaseParser);
use Bio::EnsEMBL::Registry;

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

  my $project;

  if ($file =~ /project[=][>](\S+?)[,]/) {
    $project = $1;
  }
  
  my $registry = "Bio::EnsEMBL::Registry";

  if ($project eq 'ensembl') {
    $registry->load_registry_from_multiple_dbs(
      {
        '-host' => 'mysql-ensembl-mirror.ebi.ac.uk',
	'-port'    => 4240,
        '-user' => 'ensro',
      },
    );
  } elsif ($project eq 'ensemblgenomes') {
    $registry->load_registry_from_multiple_dbs(
      {
        '-host' => 'mysql-eg-staging-1.ebi.ac.uk',
        '-port' => 4160,
        '-user' => 'ensro',
      },
      {
        '-host' => 'mysql-eg-staging-2.ebi.ac.uk',
        '-port' => 4275,
        '-user' => 'ensro',
      },
    );
  } else {
    die("Missing or unsupported project value. Supported values: ensembl, ensemblgenomes");
  }

  my %id2name = $self->species_id2name;
  my $species_name = $id2name{$species_id}[0];

  my %interpro = $self->get_core_interpro($registry, $species_name);

  my $add_interpro_sth = $self->dbi()->prepare(
    "INSERT IGNORE INTO interpro (interpro, pfam, dbtype) VALUES(?,?,?)"
  );

  my $add_xref_sth = $self->dbi()->prepare(
    "INSERT IGNORE INTO xref ".
    "(accession, label, description, source_id, species_id, info_type) ".
    "VALUES(?,?,?,?,?,?)"
  );
  
  # The InterproScan pipeline uses additional sources for Interpro
  # links (e.g. Gene3D, Panther), so in order to replicate that in the xref
  # database, remove the restriction on the column contents.
  $self->dbi()->do("ALTER TABLE interpro MODIFY COLUMN dbtype VARCHAR(25);");

  foreach my $interpro_id (sort keys %interpro) {
    foreach my $db_type (sort keys %{$interpro{$interpro_id}}) {
      foreach my $id (sort keys %{$interpro{$interpro_id}{$db_type}}) {
        my $added = 
          $add_xref_sth->execute(
            $interpro_id,
            $interpro{$interpro_id}{$db_type}{$id}{'short_name'},
            $interpro{$interpro_id}{$db_type}{$id}{'name'},
            $source_id,
            $species_id,
            $interpro{$interpro_id}{$db_type}{$id}{'info_type'},
          );

        if ( !$added ) {
          print STDERR "Problem adding '$interpro_id'\n";
          return 1;    # 1 is an error
        }

        $added =
          $add_interpro_sth->execute(
            $interpro_id,
            $id,
            $db_type
          );

        if ( !$added ) {
          print STDERR "Problem adding '$interpro_id'/".$interpro{$interpro_id}{$db_type}{'id'}."\n";
          return 1;    # 1 is an error
        }

      }
    }
  }

  return 0;
}

sub get_core_interpro {
  my ($self, $registry, $species_name) = @_;

  my $dba = $registry->get_DBAdaptor($species_name, "core");

  # Get interpro terms and related information
  my %interpro;
  my $sql =
    'select distinct '.
    'i.interpro_ac, i.id, '.
    'x.display_label, x.description, x.info_type, '.
    'a.logic_name '.
    'from xref x '.
    'inner join interpro i on x.dbprimary_acc = i.interpro_ac '.
    'inner join protein_feature pf on i.id = pf.hit_name '.
    'inner join analysis a on pf.analysis_id = a.analysis_id;';
  my $sth = $dba->dbc()->prepare($sql);
  $sth->execute();

  # Ensembl analysis logic names don't match with the terms that
  # Interpro uses, but mapping is easy enough.
  my %dbtypes = (
    'gene3d'      => 'GENE3D',
    'hmmpanther'  => 'PANTHER',
    'pfam'        => 'PFAM',
    'pfscan'      => 'PROFILE',
    'pirsf'       => 'PIRSF',
    'prints'      => 'PRINTS',
    'scanprosite' => 'PROSITE',
    'smart'       => 'SMART',
    'superfamily' => 'SSF',
    'tigrfam'     => 'TIGRFAMs',
  );

  while (my @row = $sth->fetchrow_array()) {
    my $interpro_id = $row[0];
    my $db_type = $dbtypes{$row[5]};
    my $id = $row[1];

    if (defined $db_type) {
      $interpro{$interpro_id}{$db_type}{$id}{'short_name'} = $row[2];
      $interpro{$interpro_id}{$db_type}{$id}{'name'} = $row[3];
      $interpro{$interpro_id}{$db_type}{$id}{'info_type'} = $row[4];
    }
  }

  print "Retrieved ".scalar(keys %interpro)." interpro ids.\n";
  return %interpro;
}


1;
