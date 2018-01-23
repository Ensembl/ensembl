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

package Bio::EnsEMBL::DataFile;

use strict;
use warnings;

use base qw/Bio::EnsEMBL::Storable/;

use Bio::EnsEMBL::ApiVersion;
use Bio::EnsEMBL::Utils::Argument qw/rearrange/;
use Bio::EnsEMBL::Utils::Exception qw/throw warning/;
use Bio::EnsEMBL::Utils::Scalar qw/:assert/;
use Bio::EnsEMBL::Utils::URI qw/is_uri/;
use File::Spec;

=head2 new

  Arg [-ADAPTOR]      : Bio::EnsEMBL::DBSQL::DataFileAdaptor
  Arg [-DBID]         : Integer $dbID
  Arg [-COORD_SYSTEM] : Bio::EnsEMBL::CoordSystem $coord_system
  Arg [-ANALYSIS]     : Bio::EnsEMBL::Analysis $analysis
  Arg [-NAME]         : String $name
  Arg [-VERSION_LOCK] : Boolean $version_lock
  Arg [-ABSOLUTE]     : Boolean $absolute
  Arg [-URL]          : String $url
  Arg [-FILE_TYPE]    : String $file_type
  Example			      : Bio::EnsEMBL::DataFile->new();
  Description	      : Returns a new instance of this object
  Returntype 	      : Bio::EnsEMBL::DataFile
  Exceptions 	      : Thrown if data is not as expected

=cut

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  my ($coord_system, $analysis, $name, $version_lock, $absolute, $url, $file_type) = 
    rearrange([qw/coord_system analysis name version_lock absolute url file_type/], @args);
  
  $self->coord_system($coord_system);
  $self->analysis($analysis);
  $self->name($name);
  $self->version_lock($version_lock);
  $self->absolute($absolute);
  $self->url($url);
  $self->file_type($file_type);
  
  return $self;
}


=head2 get_ExternalAdaptor

  Arg[1]     	: Scalar; (optional) base path. Uses defaults if not given 
  Arg[2]        : Scalar; (optional) file type
  Example       : my $ea = $df->get_ExternalAdaptor('/base/path');
  Description	: Delegates to the parent adaptor to retrieve the external 
                  adaptor for this data type
  Returntype 	: Adaptor; will be an adaptor that can read the given data file
  Exceptions 	: Thrown if there is no attached adaptor. 

=cut

sub get_ExternalAdaptor {
  my ($self, $base_path, $requested_type) = @_;
  my $adaptor = $self->adaptor();
  throw "No DataFileAdaptor found in this object. Cannot request ExternalAdaptor" if ! $adaptor;
  return $adaptor->DataFile_to_adaptor($self, $base_path, $requested_type);
}

=head2 path

  Arg[1]      : Scalar base of the path to use. Can be ignored if the instance
                already represents a canonical path 
  Example     : my $f = $df->path();
  Description : Used to generate the path to the file resource. Can return a
                path to the file or a URL but it is up to the using code to
                know how to interprate the different returned forms.
                
                If the data file url is canonical then this is just returned. 
                If not then a path is generated of the form 
                B</base/path/production_name/coord_system_version/[software_version]/db_group/name.ext> 
                
  Returntype  : Scalar the absolute path/url to the given resource
  Exceptions  : Thrown if the linked Coordinate System lacks a version and the
                current database also lacks a default version
  Caller      : public

=cut


sub path {
  my ($self, $base) = @_;
  my $all_paths = $self->get_all_paths($base);
  return $all_paths->[0];
}

sub get_all_paths {
  my ($self, $base) = @_;
    
  return [$self->url()] if $self->absolute();
  
  my @all_paths;
  
  $base = $self->adaptor()->get_base_path($base) if ! $base;

  my $production_name = $self->adaptor()->db()->get_MetaContainer()->get_production_name();
  my $cs_version = $self->coord_system()->version();
  if(! $cs_version) {
    my ($highest_cs) = @{$self->adaptor()->db()->get_CoordSystemAdaptor()->fetch_all()};
    $cs_version = $highest_cs->version();
  }
  if(!$cs_version) {
    my $name = $self->name();
    throw "The file '${name}' in species '${$production_name} is attached to a CoordinateSystem lacking a version and has no default assembly. Please fix";
  }
  
  my @portions;
  push(@portions, $production_name);
  push(@portions, $cs_version);
  push(@portions, software_version()) if $self->version_lock();
  push(@portions, $self->adaptor()->db()->group());
  
  #Targets are the files to generate
  my @targets;
  #If URL is populated we assume we need to add this onto the end but removing the /
  if($self->url()) {
    my @split = split(/\//, $self->url());
    push(@targets, [@split]);
  }
  else {
    my $extensions = $self->adaptor()->DataFile_to_extensions($self);
    foreach my $ext (@{$extensions}) {
      my $filename = sprintf(q{%s.%s}, $self->name(), $ext);
      push(@targets, [$filename]);
    }
  }
  
  my $is_uri = is_uri($base);
  foreach my $t (@targets) {
    my $path;
    if($is_uri) {
      $path = join(q{/}, $base, @portions, @{$t});
    }
    else {
      $path = File::Spec->catfile($base, @portions, @{$t});
    }
    push(@all_paths, $path);
  }
  return \@all_paths;
}

=head2 coord_system

  Arg[1]      : Bio::EnsEMBL::CoordSystem Optional setter  
  Description : Mutator for the coord system field. All files are linked to one
  Returntype  : Bio::EnsEMBL::CoordSystem
  Exceptions  : Thrown if not of the expected type

=cut


sub coord_system {
  my ($self, $coord_system) = @_;
  if(defined $coord_system) {
    assert_ref($coord_system, 'Bio::EnsEMBL::CoordSystem', 'coord_system');
  	$self->{'coord_system'} = $coord_system;
  }
  return $self->{'coord_system'};
}

=head2 analysis

  Arg[1]      : Bio::EnsEMBL::Analysis Optional setter  
  Description : Mutator for the analysis field. All files are linked to one
  Returntype  : Bio::EnsEMBL::Analysis
  Exceptions  : Thrown if not of the expected type

=cut

sub analysis {
  my ($self, $analysis) = @_;
  if(defined $analysis) {
    assert_ref($analysis, 'Bio::EnsEMBL::Analysis', 'analysis');
  	$self->{'analysis'} = $analysis;
  }
  return $self->{'analysis'};
}

=head2 name

  Arg[1]      : String Optional setter  
  Description : Mutator for the name of the file. Can be used in file location
                generation
  Returntype  : String

=cut

sub name {
  my ($self, $name) = @_;
  if(defined $name) {
  	$self->{'name'} = $name;
  }
  return $self->{'name'};
}

=head2 version_lock

  Arg[1]      : Boolean Optional setter  
  Description : Boolean indicating if the file is linked to the version of the
                database it was found in.
  Returntype  : Boolean

=cut

sub version_lock {
  my ($self, $version_lock) = @_;
  if(defined $version_lock) {
    assert_boolean($version_lock, 'version_lock');
    $self->{'version_lock'} = $version_lock;
  }
  return $self->{'version_lock'};
}

=head2 absolute

  Arg[1]      : Boolean Optional setter  
  Description : Indicates if the URL of this file is an absolute one i.e. 
                should be used verbatim or not.
  Returntype  : Boolean

=cut

sub absolute {
  my ($self, $absolute) = @_;
  if(defined $absolute) {
    assert_boolean($absolute, 'absolute');
  	$self->{'absolute'} = $absolute;
  }
  return $self->{'absolute'};
}

=head2 url

  Arg[1]      : String Optional setter  
  Description : Location of the file. Can be optional and if set means once
                we are in an automatic location use this value to locate
                the file.
  Returntype  : String

=cut

sub url {
  my ($self, $url) = @_;
  $self->{'url'} = $url if defined $url;
  return $self->{'url'};
}

=head2 file_type

  Arg[1]      : String Optional setter  
  Description : The type of file we are working with. Can be used to generate
                a file name.
  Returntype  : String

=cut

sub file_type {
  my ($self, $file_type) = @_;
  $self->{'file_type'} = $file_type if defined $file_type;
  return $self->{'file_type'};
}

#=head2 files
#
#  Args       	: 
#  Example			: my $files = @{$df->files()};
#  Description	: Returns all the file names we expect to cover for a flat file
#  Returntype 	: type return_description
#  Exceptions 	: 
#  Caller     	: caller
#  Status     	: status
#
#=cut
#
#
#sub files {
#  my ($self) = @_;
#  
#}

1;
