=pod

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

=head1 NAME

Bio::EnsEMBL::Pipeline::FASTA::SCPBlast

=head1 DESCRIPTION

Performs a find in the Blast index directory, for the given species and copies
them to the specified target servers.

Allowed parameters are:

=over 8

=item no_scp - If true then we will not run SCP but still finish cleanly without error

=item type - The type of dump to copy. Required parameter

=item genomic_dir - Needed if you are copying DNA genomic files

=item genes_dir - Needed if you are copying DNA gene files

=item target_servers - The servers to copy to. Expects to be an array

=item species - Species to work with

=item scp_user - The user to scp as. Defaults to the current user

=item scp_identity - Give an identity file to use during ssh commands 
                     (useful when you are not scping as yourself)

=item base_path - The base of the dumps. The source blast directory is 
                  constructed from this path

=back

=cut

package Bio::EnsEMBL::Pipeline::FASTA::SCPBlast;

use strict;
use warnings;

use base qw/Bio::EnsEMBL::Pipeline::FASTA::Base/;

use Bio::EnsEMBL::Utils::Scalar qw/check_ref/;
use File::Spec;

sub param_defaults {
  my ($self) = @_;
  return {
    no_scp => 0,
#    genomic_dir => '',
#    genes_dir => '',
#    target_servers => ['srv1', 'srv2'],

    scp_user => $ENV{USER}, #defaults to the current user
#    scp_identity => '', 

#    type => 'genes'/'genomic',
#    species => '',     
  };
}

sub fetch_input {
  my ($self) = @_;
  if($self->param('no_scp')) {
    $self->info('Skipping as no_scp has been specified');
    return;
  }
  
  my $servers = $self->param('target_servers');
  
  if(!check_ref($servers, 'ARRAY') || ! @{$servers}) {
    my $msg = 'Will not perform copy as we have no servers';
    my $is_error = 0;
    $self->db()->get_JobMessageAdaptor()->register_message(
      $self->input_job()->dbID(), $msg, $is_error
    );
    $self->info($msg);
    return;
  }
  
  foreach my $key (qw/type species/) {
    $self->throw("Key $key is required") unless $self->param($key);
  }
  my $type = $self->param('type');
  if($type ne 'genomic' && $type ne 'genes') {
    $self->throw("param 'type' must be set to 'genomic' or 'genes'");
  }
  $self->target_dir(); #prodding for fetch's sake
  return;
}

sub run {
  my ($self) = @_;
  if($self->param('no_scp')) {
    $self->info('Skipping as no_scp has been specified');
    return;
  }
  my $servers = $self->param('target_servers');
  return unless @{$servers};
  my $files = $self->get_files();
  foreach my $server (@{$servers}) {
    $self->info('Copying files to %s for species %s', $server, $self->param('species'));
    $self->copy_to_server($files, $server);
  }
  return;
}

sub write_output {
  my ($self) = @_;
  $self->cleanup_DBAdaptor();
  return;
}

sub copy_to_server {
  my ($self, $files, $server) = @_;
  my $target_dir = $self->target_dir();
  $self->check_remote_dir($target_dir, $server);
  my $user = $self->param('scp_user');
  my $identity = $self->identity_param();
  foreach my $file (@{$files}) {
    my ($volume, $directory, $filename) = File::Spec->splitpath($file);
    my $target_path = File::Spec->catfile($target_dir, $filename);
    my $cmd = sprintf('scp %s %s %s@%s:%s', $identity, $file, $user, $server, $target_path);
    $self->fine('Executing %s', $cmd);
    system($cmd) and $self->throw(sprintf("Cannot run command '%s'. RC %d", $cmd, ($?>>8)));
  }
  return;
}

sub get_files {
  my ($self) = @_;
  my $species = $self->web_name();
  my $filter = sub {
    my ($filename) = @_;
    return ($filename =~ /^$species.+fa.+$/) ? 1 : 0;
  };
  my $files = $self->find_files($self->blast_dir(), $filter);
  $self->info('Found %d file(s) to copy', scalar(@{$files}));
  return $files;
}

sub blast_dir {
  my ($self) = @_;
  return $self->get_dir('blast', $self->param('type'));
}

sub target_dir {
  my ($self) = @_;
  my $t = $self->param('type');
  my $key = "${t}_dir";
  my $dir = $self->param($key);
  $self->throw("Cannot locate the parameter $key. We expect to do so") unless $dir;
  return $dir;
}

sub check_remote_dir {
  my ($self, $remote_dir, $server) = @_;
  my ($echo_rc) = $self->ssh_cmd($server, "echo -n 1");
  $self->throw("Cannot connect to $server") if $echo_rc; #1 means fail
  my ($exists_rc) = $self->ssh_cmd($server, "test -d $remote_dir");
  if($exists_rc == 1) {
    $self->info('Directory %s does not exist on %s. Will create it');
    my ($mkdir_rc, $mkdir_out) = $self->ssh_cmd($server, "mkdir -p $remote_dir");
    if($mkdir_rc == 1) {
      $self->throw("Cannot create the directory $remote_dir on $server. Output from cmd was $mkdir_out. Check and rerun");
    }
  }
  return;
}

sub ssh_cmd {
  my ($self, $server, $cmd) = @_;
  my $user = $self->param('scp_user');
  my $identity = $self->identity_param();
  $self->fine("Running command '%s' on '%s' as user '%s'", $cmd, $server, $user);
  my $ssh_cmd = sprintf('ssh %s %s@%s "%s"', $identity, $user, $server, $cmd);
  my $output = `$ssh_cmd`;
  my $rc = $? >> 8;
  return ($rc, $output);
}

sub identity_param {
  my ($self) = @_;
  return ($self->param('scp_identity')) ? '-i '.$self->param('scp_identity') : q{};
}

1;
