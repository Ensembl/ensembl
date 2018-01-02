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

Bio::EnsEMBL::DBFile::FileAdaptor - Base Adaptor for direct file access

=head1 DESCRIPTION

Basic wrapper class to provide access to file based data.

This is primarily aimed at indexed Collection(.col) files which are optimised for Slice 
based queries. Collections store fixed width width/windowed data as BLOBS.  This makes 
it possible to seek to the a required location given slice coordinate and read the only 
the required amount of data covering the slice.

Currently only works as hybrid DBAdaptor e.g. ResultFeatureAdaptor which inherits both from 
here and BaseFeatureAdaptor.

=cut



package Bio::EnsEMBL::DBFile::FileAdaptor;

use Bio::EnsEMBL::Utils::Exception qw(throw warning deprecate);
use strict;
use warnings;


=head2 get_filehandle

  Arg[1]     : string     - filepath
  Arg[2]     : HASHREF    - Optional params, see open_file
  Example    : my $fh     = $self->get_filehandle($filepath, 1);
  Description: Gets and caches a simple file handle.
  Returntype : GLOB/undef - filehandle
  Exceptions : warns if cache entry exists but is not defined 
  Caller     : general
  Status     : at risk

=cut

sub get_filehandle{
  my ($self, $filepath, $params_hash) = @_;

  my $file_op = '<';

  if(exists $params_hash->{-file_operator}){
	$file_op = $params_hash->{-file_operator};
  }else{
	$params_hash->{-file_operator} = $file_op;
  }

  if(! exists $self->{file_cache}{$filepath}{filehandle}){
	my $fh = $self->Bio::EnsEMBL::DBFile::FileAdaptor::open_file($filepath, $params_hash);

	if(defined $fh){
	  $self->{file_cache}{$filepath}{filehandle} = $fh;
	  #$self->initialise_filehandle($filepath) if $self->can('initialise_filehandle');
	  $self->initialise_filehandle($filepath) if($file_op eq '<');
	}
  }
  elsif(! defined $self->{file_cache}{$filepath}{filehandle}){
	#This maybe one of several read/seek errors which will have already been warned
	warn "Encountered and error with file handle for $filepath\n";
  }
  #else
  # check against cache file op
  # to make sure we aren't trying to open an already open fh with a different operator

 
  return $self->{file_cache}{$filepath}{filehandle};
}


=head2 open_file

  Arg[1]     : string     - filepath
  Arg[2]     : HASHREF    - Optional params:
                          -binmode       => 0|1,   # Boolean i.e. treat file as binary
                          -file_operator => '>'    # Default is '<'
                         #-perms_octal   =>  # Requires FileHandle
  Example    : my $fh     = $self->open_file($filepath, {-binmode = > 1, -file_operator => '>'});
  Description: Opens a file for reading or writing.
  Returntype : GLOB/undef - filehandle
  Exceptions : warns if file open fails
               warns if file operator unsupported
               warns if failed to set binmode
  Caller     : general
  Status     : at risk

=cut

sub open_file{
  my ($self, $filepath, $params_hash) = @_;

  #Validate params_hash? 
  #rearrange? Will not warn/throw for invalid keys?
  #perms octal, requires FileHandle? See EFGUtils::open_file



  my $file_op = $params_hash->{-file_operator} || '<';

  if(($file_op ne '<') &&
	 ($file_op ne '>') &&
	 ($file_op ne '>>')){
	throw("Cannot perform open with unsupported operator:\t${file_op}${filepath}");
  }

  my $fh;
  my $success = open $fh, $file_op, $filepath;
  #$fh will be still be GLOB on fail
  
  #These warn instead of throw/die to allow
  #open_file to be used to test a file
  #this prevents throws/die when an attempting to access an absent file (good for webcode)
  #could alternatively change to throw/die and eval where required
  #prevents need to catch everywhere else and potential double reporting of error

  if(! $success){
	#undef $fh;
	throw("Failed to open:\t$filepath\n$!\n");
  }
  elsif($params_hash->{-binmode}){
	$success = binmode $fh;
	  
	if(! $success){
	  throw("Failed to set binmode:\t$filepath\n$!");
	  #undef $fh;
	}
  }

  return $fh;
}


=head2 validate_file_length

  Arg[1]     : string  - filepath
  Arg[2]     : int     - expected length in bytes
  Example    : $self->validate_file_length($filepath, $expected_length);
  Description: Utility method which can be used during file creation
  Returntype : None
  Exceptions : warns if file open fails
               throws if file is not expected length
  Caller     : general
  Status     : at risk - change to seek to accounts for 'logical characters'

=cut

sub validate_file_length{
  my ($self, $filepath, $expected_length, $binmode) = @_;

  #Currently not using cache as we rarely want to 
  #use the file handle afterwards


  #THIS WAS USING EFGUtils::open_file imported in the Collector::ResultFeature!!!!
  #which is just a sub not a class method, and is in a parallel inheritance path
  #No warnings about redefining method :(
  #Force use of FileAdaptor::open_file

  my $fh = $self->Bio::EnsEMBL::DBFile::FileAdaptor::open_file($filepath, {-binmode => $binmode});


  #sysseek always returns length in bytes, change to seek which 
  #uses logical characters i.e. actual encoding?
  #Does seek use bytes in binmode and chars in non-binmode?

  my $seeked_bytes = sysseek($fh, 0, 2);# 2 is SEEK_END
  #There is no systell function. Use sysseek(FH, 0, 1) for that.

  if($seeked_bytes < $expected_length){
	throw("File is shorter($seeked_bytes) than expected($expected_length):\t$filepath\n");
  }
  elsif($seeked_bytes > $expected_length){
	throw("File is longer($seeked_bytes) than expected($expected_length):\t$filepath\n");
  }
 
  return;
}





### STUBB/TEMPLATE METHODS ###
#
#   If required hese should be over-ridden in the 
#   descendant FileAdaptor e.g. CollectionAdaptor
#   Listed here rather for visibility (rather than 
#   using 'can')


sub initialise_filehandle{
  return;
}



1;
