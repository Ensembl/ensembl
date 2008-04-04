package Bio::EnsEMBL::IdMapping::ResultAnalyser;

=head1 NAME


=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS


=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http:#www.ensembl.org/info/about/code_licence.html for details.

=head1 AUTHOR

Patrick Meidl <meidl@ebi.ac.uk>, Ensembl core API team

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<ensembl-dev@ebi.ac.uk>

=cut


use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::IdMapping::BaseObject;
our @ISA = qw(Bio::EnsEMBL::IdMapping::BaseObject);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::ScriptUtils qw(path_append);


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

  # classify source genes by type (status-logic_name-biotype)
  $self->classify_source_genes_by_type;

  # classify source genes by mapping status
  $self->classify_genes_by_mapping($gene_mappings, $similarity_events);
}

#
# Analyse stable ID data from two existing dbs.
# This is for potential stand-alone use of this module.
#
# [todo]
sub analyse_db {
}


sub classify_source_genes_by_type {
  my $self = shift;

  foreach my $s_gene (values %{ $self->cache->get_by_name('genes_by_id', 'source') }) {
    $self->add('source', $self->class_key($s_gene), 'all', $s_gene->stable_id);
  }
}


sub classify_target_genes_by_type {
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
    unless ($t_gene->status and $t_gene->logic_name and $t_gene->biotype) {
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


#
# genes will be classified as:
#   - mapped
#   - deleted
#     - lost_similar
#     - lost_definite
#
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


#
# Add a stable ID / property pair to a name/dbtype lookup hash.
#
# This datastructure is a bloat for some applications, but is general enough to
# be used as a lookup hash and to generate statistics (counts by type) and
# debug lists (dump by type).
#
sub add {
  my ($self, $dbtype, $class, $subclass, $stable_id, $val) = @_;

  # private method, so no argument check done for performance reasons

  # default to a value of '1'
  $val = 1 unless (defined($val));

  $self->{$dbtype}->{$class}->{$subclass}->{$stable_id} = $val;
}


sub get {
  my ($self, $dbtype, $class, $subclass, $stable_id) = @_;

  # private method, so no argument check done for performance reasons

  return $self->{$dbtype}->{$class}->{$subclass}->{$stable_id};
}


sub get_all_by_subclass {
  my ($self, $dbtype, $class, $subclass) = @_;

  # argument check
  throw("Need a dbtype (source|target).") unless ($dbtype);
  throw("Need a class.") unless ($class);
  throw("Need a subclass.") unless ($subclass);

  return [ keys %{ $self->{$dbtype}->{$class}->{$subclass} || {} } ];
}


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


sub get_count_by_subclass {
  my ($self, $dbtype, $class, $subclass) = @_;

  # argument check
  throw("Need a dbtype (source|target).") unless ($dbtype);
  throw("Need a class.") unless ($class);
  throw("Need a subclass.") unless ($subclass);

  return scalar(keys %{ $self->{$dbtype}->{$class}->{$subclass} || {} });
}


sub get_count_by_class {
  my ($self, $dbtype, $class) = @_;

  # argument check
  throw("Need a dbtype (source|target).") unless ($dbtype);
  throw("Need a class.") unless ($class);

  return scalar(@{ $self->get_all_by_class($dbtype, $class) });
}


sub get_all_classes {
  my ($self, $dbtype) = @_;

  # argument check
  throw("Need a dbtype (source|target).") unless ($dbtype);

  return [ sort keys %{ $self->{$dbtype} || {} } ];
}


sub class_key {
  my ($self, $gene) = @_;
  return join('-', map { $gene->$_ } qw(status logic_name biotype));
}


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


sub create_summary_email {
  my $self = shift;
  
  my $fh = $self->get_filehandle('summary_email.txt');

  #
  # title
  # 
  print $fh qq(Stable ID mapping results\n);
  print $fh qq(=========================\n\n);

  #
  # timing
  #
  print $fh "Run at:  ".localtime;
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

