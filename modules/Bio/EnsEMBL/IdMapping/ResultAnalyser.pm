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

Bio::EnsEMBL::IdMapping::ResultAnalyser - analyse stable Id mapping results

=head1 SYNOPSIS

  # get a result analyser
  my $analyser = Bio::EnsEMBL::IdMapping::ResultAnalyser->new(
    -LOGGER => $logger,
    -CONF   => $conf,
    -CACHE  => $cache
  );

  # analyse results
  $analyser->analyse( $gene_mappings,
    $stable_id_mapper->get_all_stable_id_events('similarity') );

  # write results to file
  $analyser->write_results_to_file;

  # create click lists
  $analyser->create_clicklist;

  # mapping_summary
  $analyser->create_mapping_summary;

=head1 DESCRIPTION

This is a utility module which analyses the stable Id mapping results
by providing various sorts of mapping statistics. It also creates
clicklists and a mapping summary.

=head1 METHODS

  analyse
  analyse_db
  classify_source_genes_by_type
  classify_genes_by_mapping_simple
  classify_genes_by_mapping
  add
  get
  get_all_by_subclass
  get_all_by_class
  get_count_by_subclass
  get_count_by_class
  get_all_classes
  class_key
  write_results_to_file
  create_clicklist
  create_mapping_summary
  read_from_file

=cut


package Bio::EnsEMBL::IdMapping::ResultAnalyser;

use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::IdMapping::BaseObject;
our @ISA = qw(Bio::EnsEMBL::IdMapping::BaseObject);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::ScriptUtils qw(path_append);


=head2 analyse

  Arg[1]      : Bio::EnsEMBL::IdMapping::MappingList $gene_mappings - the gene
                mappings to analyse
  Arg[2]      : Arrayref of Strings - similarity events
  Example     : $analyser->analyse($gene_mappings,
                  $stable_id_mapper->get_all_stable_id_events('similarity'));
  Description : Analyses the results of a stable Id mapping run.
  Return type : none
  Exceptions  : thrown on wrong or missing arguments
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub analyse {
  my $self = shift;
  my $gene_mappings = shift;
  my $similarity_events = shift;
  
  # argument check
  unless ($gene_mappings and
          $gene_mappings->isa('Bio::EnsEMBL::IdMapping::MappingList')) {
    throw("Need a Bio::EnsEMBL::IdMapping::MappingList of genes.");
  }

  unless ($similarity_events and ref($similarity_events) eq 'ARRAY') {
    throw("Need a list of similarity events.");
  }

  # classify source genes by type (logic_name-biotype)
  $self->classify_source_genes_by_type;

  # classify source genes by mapping status
  $self->classify_genes_by_mapping($gene_mappings, $similarity_events);
}


=head2 classify_source_genes_by_type

  Example     : $analyser->classify_source_genes_by_type;
  Description : Classifies source genes by type and adds them to the internal
                datastructure. For the format of the classification string see
                class_key().
  Return type : none
  Exceptions  : none
  Caller      : internal
  Status      : At Risk
              : under development

=cut

sub classify_source_genes_by_type {
  my $self = shift;

  foreach my $s_gene (values %{ $self->cache->get_by_name('genes_by_id', 'source') }) {
    $self->add('source', $self->class_key($s_gene), 'all', $s_gene->stable_id);
  }
}


=head2 classify_genes_by_mapping_simple

  Arg[1]      : Bio::EnsEMBL::IdMapping::MapppingList $gene_mappings - gene
                mappings to classify
  Example     : $analyser->classify_genes_by_mapping_simple;
  Description : Classifies target genes by mapping ('mapped' or 'unmapped').
  Return type : none
  Exceptions  : thrown on wrong or missing argument
  Caller      : This method is not in use at the momen.
  Status      : At Risk
              : under development

=cut

sub classify_genes_by_mapping_simple {
  my $self = shift;
  my $gene_mappings = shift;

  # argument check
  unless ($gene_mappings and
          $gene_mappings->isa('Bio::EnsEMBL::IdMapping::MappingList')) {
    throw("Need a Bio::EnsEMBL::IdMapping::MappingList of genes.");
  }

  my %result = ();
  
  # firrst, create a lookup hash of source genes by target internal ID
  my %source_genes_by_target = ();
  foreach my $e (@{ $gene_mappings->get_all_Entries }) {
    my $s_gene = $self->cache->get_by_key('genes_by_id', 'source', $e->source);
    my $t_gene = $self->cache->get_by_key('genes_by_id', 'target', $e->target);
    $source_genes_by_target{$t_gene->id} = $s_gene;
  }

  # now loop over target genes
  foreach my $t_gene (values %{ $self->cache->get_by_name('genes_by_id', 'target') }) {
  
    # check if target gene has all required properties set
    unless ($t_gene->logic_name and $t_gene->biotype) {
      $self->logger->warning("Missing data for target gene: ".
        $t_gene->to_string."\n", 1);
    }

    my $class = $self->class_key($t_gene);

    # classify as '1' if mapped (using source gene's stable ID), otherwise '0'
    if (my $s_gene = $source_genes_by_target{$t_gene->id}) {
      $self->add('target', $class, 'mapped', $s_gene->stable_id);
    } else {
      $self->add('target', $class, 'unmapped', $t_gene->stable_id);
    }

  }
}


=head2 classify_genes_by_mapping

  Arg[1]      : Bio::EnsEMBL::IdMapping::MapppingList $gene_mappings - gene
                mappings to classify
  Arg[2]      : Arrayref of Strings - similarity events
  Example     : $analyser->classify_genes_by_mapping;
  Description : Classifies genes by mapping. Status is
                  'mapped' => stable Id was mapped
                  'lost_similar' => stable Id not mapped, but there is a
                                    similarity entry for the source Id
                  'lost_definite' => not mapped and no similarity
  Return type : none
  Exceptions  : thrown on wrong or missing argument
  Caller      : This method is not in use at the momen.
  Status      : At Risk
              : under development

=cut

sub classify_genes_by_mapping {
  my $self = shift;
  my $gene_mappings = shift;
  my $similarity_events = shift;
  
  # argument check
  unless ($gene_mappings and
          $gene_mappings->isa('Bio::EnsEMBL::IdMapping::MappingList')) {
    throw("Need a Bio::EnsEMBL::IdMapping::MappingList of genes.");
  }

  unless ($similarity_events and ref($similarity_events) eq 'ARRAY') {
    throw("Need a list of similarity events.");
  }

  # mapped genes
  foreach my $e (@{ $gene_mappings->get_all_Entries }) {
    my $s_gene = $self->cache->get_by_key('genes_by_id', 'source', $e->source);
    $self->add('source', $self->class_key($s_gene), 'mapped',
      $s_gene->stable_id);
  }

  # lookup hash for similarities
  my %similar = ();
  foreach my $event (@{ $similarity_events }) {
    my ($stable_id) = split("\t", $event);
    $similar{$stable_id} = 1;
  }
  
  # deleted genes
  foreach my $s_gene (values %{ $self->cache->get_by_name('genes_by_id', 'source') }) {
    
    my $stable_id = $s_gene->stable_id;
    my $class = $self->class_key($s_gene);

    unless ($self->get('source', $class, 'mapped', $stable_id)) {

      # sub-classify as 'lost_similar' or 'lost_definite'
      if ($similar{$stable_id}) {
        $self->add('source', $class, 'lost_similar', $stable_id);
      } else {
        $self->add('source', $class, 'lost_definite', $stable_id);
      }

    }
  }
  
}


=head2 add

  Arg[1]      : String $dbtype - db type ('source' or 'target')
  Arg[2]      : String $class - key identifying a gene type (see class_key())
  Arg[3]      : String $subclass - status identifier (e.g. 'mapped', 'lost')
  Arg[4]      : String $stable_id - gene stable Id
  Arg[5]      : String $val - value (usually 0 or 1)
  Example     : $analyser->add('source', 'KNOWN-ensembl-protein_coding',
                  'mapped', 'ENSG00002342', 1);
  Description : Add a stable Id / property pair to a name/dbtype lookup hash.
  
                The datastructure is a bit of a bloat, but is general enough to
                be used as a lookup hash and to generate statistics (counts by
                type) and debug lists (dump by type).
  Return type : String - the added value
  Exceptions  : none
  Caller      : internal
  Status      : At Risk
              : under development

=cut

sub add {
  my ($self, $dbtype, $class, $subclass, $stable_id, $val) = @_;

  # private method, so no argument check done for performance reasons

  # default to a value of '1'
  $val = 1 unless (defined($val));

  $self->{$dbtype}->{$class}->{$subclass}->{$stable_id} = $val;
}


=head2 get

  Arg[1]      : String $dbtype - db type ('source' or 'target')
  Arg[2]      : String $class - key identifying a gene type (see class_key())
  Arg[3]      : String $subclass - status identifier (e.g. 'mapped', 'lost')
  Arg[4]      : String $stable_id - gene stable Id
  Example     : my $mapping_status = $analyser->get('source',
                  'KNOWN-ensembl-protein_coding', 'mapped', 'ENSG00002342');
  Description : Gets a stable Id mapping status from the internal datastructure.
  Return type : String
  Exceptions  : none
  Caller      : internal
  Status      : At Risk
              : under development

=cut

sub get {
  my ($self, $dbtype, $class, $subclass, $stable_id) = @_;

  # private method, so no argument check done for performance reasons

  return $self->{$dbtype}->{$class}->{$subclass}->{$stable_id};
}


=head2 get_all_by_subclass

  Arg[1]      : String $dbtype - db type ('source' or 'target')
  Arg[2]      : String $class - key identifying a gene type (see class_key())
  Arg[3]      : String $subclass - status identifier (e.g. 'mapped', 'lost')
  Example     : my @mapped_stable_ids = @{
                  $analyser->get_all_by_subclass(
                    'source', 'KNOWN-ensembl-protein_coding',
                    'mapped'
                  ) };
  Description : Gets a list of stable Id for a given subclass.
  Return type : Arrayref of String (stable Ids)
  Exceptions  : thrown on missing arguments
  Caller      : internal
  Status      : At Risk
              : under development

=cut

sub get_all_by_subclass {
  my ($self, $dbtype, $class, $subclass) = @_;

  # argument check
  throw("Need a dbtype (source|target).") unless ($dbtype);
  throw("Need a class.") unless ($class);
  throw("Need a subclass.") unless ($subclass);

  return [ keys %{ $self->{$dbtype}->{$class}->{$subclass} || {} } ];
}


=head2 get_all_by_class

  Arg[1]      : String $dbtype - db type ('source' or 'target')
  Arg[2]      : String $class - key identifying a gene type (see class_key())
  Example     : my @stable_ids = @{
                  $analyser->get_all_by_class( 'source',
                    'KNOWN-ensembl-protein_coding' ) };
  Description : Gets a list of stable Id for a given class.
  Return type : Arrayref of String (stable Ids)
  Exceptions  : thrown on missing arguments
  Caller      : internal
  Status      : At Risk
              : under development

=cut

sub get_all_by_class {
  my ($self, $dbtype, $class) = @_;

  # argument check
  throw("Need a dbtype (source|target).") unless ($dbtype);
  throw("Need a class.") unless ($class);

  my %merged = ();

  foreach my $subclass (keys %{ $self->{$dbtype}->{$class} || {} }) {
    while (my ($key, $val) = each(%{ $self->{$dbtype}->{$class}->{$subclass} })) {
      $merged{$key} = $val;
    }
  }

  return [ keys %merged ];
}


=head2 get_count_by_subclass

  Arg[1]      : String $dbtype - db type ('source' or 'target')
  Arg[2]      : String $class - key identifying a gene type (see class_key())
  Arg[3]      : String $subclass - status identifier (e.g. 'mapped', 'lost')
  Example     : my $num_mapped = $analyser->get_count_by_subclass('source',
                  'KNOWN-ensembl-protein_coding', 'mapped');
  Description : Gets the number of stable Ids for a given subclass.
  Return type : Int
  Exceptions  : thrown on missing arguments
  Caller      : internal
  Status      : At Risk
              : under development

=cut

sub get_count_by_subclass {
  my ($self, $dbtype, $class, $subclass) = @_;

  # argument check
  throw("Need a dbtype (source|target).") unless ($dbtype);
  throw("Need a class.") unless ($class);
  throw("Need a subclass.") unless ($subclass);

  return scalar(keys %{ $self->{$dbtype}->{$class}->{$subclass} || {} });
}


=head2 get_count_by_class

  Arg[1]      : String $dbtype - db type ('source' or 'target')
  Arg[2]      : String $class - key identifying a gene type (see class_key())
  Example     : my $num_mapped = $analyser->get_count_by_class('source',
                  'KNOWN-ensembl-protein_coding');
  Description : Gets the number of stable Ids for a given class.
  Return type : Int
  Exceptions  : thrown on missing arguments
  Caller      : internal
  Status      : At Risk
              : under development

=cut

sub get_count_by_class {
  my ($self, $dbtype, $class) = @_;

  # argument check
  throw("Need a dbtype (source|target).") unless ($dbtype);
  throw("Need a class.") unless ($class);

  return scalar(@{ $self->get_all_by_class($dbtype, $class) });
}


=head2 get_all_classes

  Arg[1]      : String $dbtype - db type ('source' or 'target')
  Example     : foreach my $class (@{ $analyser->get_all_classes('source') }) {
                  print "$class\n";
                }
  Description : Gets a list of classes in the ResultAnalyser.
  Return type : Arrayref of String
  Exceptions  : thrown on missing argument
  Caller      : internal
  Status      : At Risk
              : under development

=cut

sub get_all_classes {
  my ($self, $dbtype) = @_;

  # argument check
  throw("Need a dbtype (source|target).") unless ($dbtype);

  return [ sort keys %{ $self->{$dbtype} || {} } ];
}


=head2 class_key

  Arg[1]      : Bio::EnsEMBL::IdMapping::TinyGene $gene - a gene object
  Example     : my $class = $analyser->class_key($gene);
  Description : Generates a key identifying a gene class. This identifier is 
                composed from the gene's logic naame, and biotye.
  Return type : String
  Exceptions  : none
  Caller      : internal
  Status      : At Risk
              : under development

=cut

sub class_key {
  my ($self, $gene) = @_;
  return join('-', map { $gene->$_ } qw(logic_name biotype));
}


=head2 write_results_to_file

  Example     : $analyser->write_results_to_file;
  Description : Writes the results of the result analysis to a file. This is a 
                human-readable text detailing the mapping statistics.
  Return type : none
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub write_results_to_file { 
  my $self = shift;

  my $fh = $self->get_filehandle('gene_detailed_mapping_stats.txt', 'stats');

  my $fmt1 = "%-60s%-16s%-16s%-16s\n";
  my $fmt2 = "%-60s%5.0f (%7s) %5.0f (%7s) %5.0f (%7s)\n";
  my $fmt3 = "%3.2f%%";

  print $fh "Gene detailed mapping results:\n\n";

  print $fh sprintf($fmt1, "Gene type", "mapped", "lost (similar)",
    "lost (definite)");

  print $fh ('-'x108), "\n";

  foreach my $class (@{ $self->get_all_classes('source') }) {
    next if ($class eq 'all');

    my $total = $self->get_count_by_class('source', $class);

    # avoid division by zero error
    unless ($total) {
      $self->logger->warning("No count found for $class.\n", 1);
      next;
    }
    
    my $mapped = $self->get_count_by_subclass('source', $class, 'mapped');
    my $similar = $self->get_count_by_subclass('source', $class,
      'lost_similar');
    my $lost = $self->get_count_by_subclass('source', $class, 'lost_definite');

    print $fh sprintf($fmt2,
                      $class,
                      $mapped,  sprintf($fmt3, $mapped/$total*100),
                      $similar, sprintf($fmt3, $similar/$total*100),
                      $lost,    sprintf($fmt3, $lost/$total*100));
  }

  close($fh);
}


=head2 create_clicklist

  Example     : $analyser->create_clicklist;
  Description : Writes an html file which contains a list of all lost genes,
                with hyperlinks to the appropriate archive website. This is to
                manually check lost genes.
  Return type : none
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub create_clicklist {
  my $self = shift;

  my $fh = $self->get_filehandle('genes_lost.html', 'stats');

  # start html output
  print $fh qq(<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">\n);
  print $fh qq(<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en-gb"  lang="en-gb">);
  print $fh "<head>\n";
  print $fh "<title>Lost genes ";
  print $fh $self->conf->param('sourcedbname'), ' -&gt; ',
            $self->conf->param('targetdbname');
  print $fh "</title>\n";
  print $fh "</head>\n<body>\n";

  my $prefix = $self->conf->param('urlprefix');
  unless ($prefix) {
    $self->logger->warning("No urlprefix set, clicklists might not be useable.\n", 1);
  }

  my $navigation;
  my $clicklist;

  foreach my $class (@{ $self->get_all_classes('source') }) {
    next if ($class eq 'all');

    $navigation .= "$class ";
    $clicklist .= "<h1>$class</h1>\n";
    
    foreach my $subclass (qw(lost_similar lost_definite)) {

      # navigation
      $navigation .= qq(<a href="#${class}-$subclass">$subclass</a> );
      
      # clicklist
      $clicklist .= "<h2>$subclass</h2>\n";

      foreach my $stable_id (@{ $self->get_all_by_subclass('source', $class, $subclass) }) {
        $clicklist .= qq(<a href="${prefix}$stable_id">$stable_id</a><br />\n);
      }
      
    }

    $navigation .= "<br />\n";
  }

  # print navigation and clicklist
  print $fh "$navigation\n\n";
  print $fh "$clicklist\n\n";

  # html footer
  print $fh "</body></html>\n";

  close($fh);
}


=head2 create_mapping_summary

  Example     : $analyser->create_mapping_summary();
  Description : Writes a text file containing a summary of the mapping stats.
                This will be emailed to the genebuilder for evaluation (you will
                have to manually send the email, using the text in
                "mapping_summary.txt" as the template).
  Return type : none
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub create_mapping_summary {
  my $self = shift;
  
  my $fh = $self->get_filehandle('mapping_summary.txt');

  #
  # title
  # 
  print $fh qq(Stable ID mapping results\n);
  print $fh qq(=========================\n\n);

  #
  # timing
  #
  print $fh "Run at:  ".localtime()."\n";
  print $fh "Runtime: ";
  print $fh $self->logger->runtime, "\n\n";

  #
  # parameters used for this run
  #
  print $fh $self->conf->list_param_values;
  print $fh "\n";

  #
  # mapping stats
  #
  foreach my $type (qw(exon transcript translation gene gene_detailed)) {
    my $filename = "${type}_mapping_stats.txt";
    
    if ($self->file_exists($filename, 'stats')) {
      print $fh $self->read_from_file($filename, 'stats');
      print $fh "\n\n";
    } else {
      print $fh "No mapping stats found for $type.\n\n";
    }
  }

  #
  # db uploads
  #
  my @uploads = (
    ['stable_ids'  => 'Stable IDs'],
    ['events'      => 'Stable ID events and mapping session'],
    ['archive'     => 'Gene and peptide archive'],
  );
  
  my $fmt1 = "%-40s%-20s\n";

  print $fh qq(Data uploaded to db:\n);
  print $fh qq(====================\n\n);

  if ($self->conf->param('dry_run')) {
   
    print $fh "None (dry run).\n";
  
  } else {
  
    foreach my $u (@uploads) {
      my $uploaded = 'no';
      $uploaded = 'yes' if ($self->conf->is_true("upload_".$u->[0]));
      print $fh sprintf($fmt1, $u->[1], $uploaded);
    }
    
  }

  print $fh "\n";

  #
  # stats and clicklist
  #
  my @output = (
    ['stats'    => 'statistics (including clicklists of deleted IDs)'],
    ['debug'    => 'detailed mapping output for debugging'],
    ['tables'   => 'data files for db upload'],
  );
  
  my $fmt2 = "%-20s%-50s\n";

  print $fh qq(\nOutput directories:\n);
  print $fh qq(===================\n\n);

  print $fh sprintf($fmt2, qw(DIRECTORY DESCRIPTION));
  print $fh ('-'x72), "\n";

  print $fh sprintf($fmt2, 'basedir', $self->conf->param('basedir'));

  foreach my $o (@output) {
    print $fh sprintf($fmt2, '$basedir/'.$o->[0], $o->[1]);
  }

  print $fh "\n";

  #
  # clicklist of first 10 deleted genes
  #
  print $fh qq(\nFirst 10 deleted known genes:\n);
  print $fh qq(=============================\n\n);

  my $in_fh = $self->get_filehandle('genes_lost.txt', 'debug', '<');
  my $prefix = $self->conf->param('urlprefix');
  my $i;
  
  while (<$in_fh>) {
    last if (++$i > 10);
    
    chomp;
    my ($stable_id, $type) = split(/\s+/);
    
    next unless ($type eq 'known');

    print $fh sprintf($fmt2, $stable_id, "${prefix}$stable_id");
  }

  close($in_fh);
  close($fh);
}


=head2 read_from_file

  Arg[1]      : String $filename - name of file to read
  Arg[2]      : (optional) String $append - directory name to append to basedir
  Example     : my $stats_text = $analyser->read_from_file('gene_mapping_stats',
                  'stats');
  Description : Reads mapping stats from a file.
  Return type : String
  Exceptions  : none
  Caller      : internal
  Status      : At Risk
              : under development

=cut

sub read_from_file {
  my $self = shift;
  my $filename = shift;
  my $append = shift;

  my $in_fh = $self->get_filehandle($filename, $append, '<');

  my $txt;

  while (<$in_fh>) {
    $txt .= $_;
  }

  return $txt;
}

1;

