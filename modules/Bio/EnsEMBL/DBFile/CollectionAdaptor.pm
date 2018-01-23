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

Bio::EnsEMBL::DBFile::CollectionAdaptor

=head1 SYNOPSIS

For use with a Bio::EnsEMBL::Collector e.g.

    package Bio::EnsEMBL::Funcgen::DBSQL::ResultFeatureAdaptor;

    @ISA = qw(Bio::EnsEMBL::Funcgen::DBSQL::BaseFeatureAdaptor 
              Bio::EnsEMBL::Funcgen::Collector::ResultFeature 
              Bio::EnsEMBL::DBFile::CollectionAdaptor);
    #DBSQL and DBFile inheritance here due to dynamic nature of ResultFeatureAdaptor


Fetch wrapper methods access file based data via read_collection_blob:

    sub _fetch_from_file_by_Slice_ResultSet{

	    #define filepath/config

        my $packed_scores =  $self->read_collection_blob(
		    										   $filepath,
			    									   $efg_sr_id,
				    								   $conf->{$window_size}{'byte_offset'},
					    							   $conf->{$window_size}{'byte_length'},
						    						  );

        #Do unpacking and object creation here

    }

=head1 DESCRIPTION

Adaptor for direct collection(.col) file access, which are binary compressed fixed 
width format files providing window based values across the genome. Collection files
integrate an index block which contains seq_region byte off set values.

NOTE: By default all collection files are generated and packed using little endian encoding. 
Due to the lack of standards of float encoding(wrt to endianess) perl packs using the 
implicit endianess of the underlying architecture. This means that accessing float
collection files located on a big endian architecture will produce unexpected results.

# endian issues will disappear with knetfile xsubs

=head1 SEE ALSO

Bio::EnsEMBL::DBFile::FileAdaptor

=cut



package Bio::EnsEMBL::DBFile::CollectionAdaptor;

use strict;
use warnings;

use Bio::EnsEMBL::DBFile::FileAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning deprecate);
use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::DBFile::FileAdaptor);


=head2 initialise_filehandle

  Arg[1]     : string  - filepath
  Example    : $self->initialise_filehandle($filepath);
  Description: Initialises the filehandle for use, in this case reads 
               the index (seq_region offsets)
  Returntype : None
  Exceptions : warns if read fails
  Caller     : Bio::EnsEMBL::DBFile::FileAdaptor::get_filehandle
  Status     : at risk

=cut

sub initialise_filehandle{
  my ($self, $filepath) = @_;
  my $fh = $self->{file_cache}{$filepath}{filehandle};
  
  #offsets include the length of the complete index block
  my ($index_size, $read_bytes, $index, $num_keys, %offset_index);
  
  ### INDEX FORMAT ###
  #First block of the index the index size in bytes(not inc size block).
  #
  #Rest of index is a hash of sr_id(v 2 bytes) key offset(V 4 bytes) value pairs
  #V (long) is 4 bytes(via sys/read), which is actually an Config{intsize} i.e. i? 
  #long is 8 bytes according to Config{longsize}!

  #read uses logical characters not necessarily in bytes
  #altho this does seem to read bytes, maybe due to binmode?
  #seek is in bytes
  #Changed to sysread/read which both use bytes explicitly
  #Can't mix sysread/seek due to I/O buffering differences

  
  #Read index_size first encoded as v(2 bytes)
  $read_bytes = sysread($fh, $index_size, 2);
    
  if(! ((defined $read_bytes) && ($read_bytes == 2))){
	#! defined is error 0 is end of file
	warn "Failed to read index size from $filepath\n$!";

	#Delete fh as it is useless/unsafe to retry
	undef $self->{file_cache}{$filepath}{filehandle};
  }
  else{	#Read index
	($index_size) = unpack('v', $index_size);
	$read_bytes = sysread($fh, $index, $index_size);  #Now read index proper
	
	if(! ((defined $read_bytes) && ($read_bytes == $index_size))){
	  #! defined is error 0 is end of file
	  warn "Failed to read index from $filepath\n$!";

	  #Delete fh as it is useless/unsafe to retry
	  undef $self->{file_cache}{$filepath}{filehandle};
	}
	else{
	  #Number of key-value pairs => $index_size /(size of key(v 2bytes) + size of offset(V 4bytes))
	  $num_keys        = $index_size/6;
	  my $unpack_template = '(vV)'.$num_keys,;
	  
	  %offset_index = unpack($unpack_template, $index);
	  $self->{file_cache}{$filepath}{off_sets} = \%offset_index;
	}
  }

  return $self->{file_cache}{$filepath}{off_sets};
}


=head2 read_collection_blob

  Arg[1]     : string - filepath
  Arg[2]     : int    - seq_region_id
  Arg[3]     : int    - seq_region offset. The byte offset required to
                        locate the required start position
  Arg[4]     : int    - byte length to read
  Example    : my $blob_substr = $self->read_collection_blob($filepath,
                                                             $sr_key,
                                                             $sr_offset,
                                                             $byte_length);
  Description: Reads bytes from file given a seq_region_key, byte offset and byte length.
               Sets filehandle to undef if read fails.
  Returntype : string - packed binary data
  Exceptions : warns if seek or read errors
  Caller     : general e.g. fetch_from_file_by_Slice_ResultSet
  Status     : at risk

=cut

# We could change this to take a Slice, hence we could check 
# whether an EOF error is because the slice is out of range 
# and undef only if it is in range i.e. the index/file is corrupt
# overkill?
# This is something the Slice API should warn about
# but will still cause undef'd filehandle here
# Index should also contain ends, so we can validate whether the slice is out of range???


sub read_collection_blob{
  my($self, $filepath, $sr_key, $sr_offset, $byte_length) = @_;
	
  my $blob_substr;
  my $fh = $self->get_filehandle($filepath, {-binmode => 1});

  if(defined $fh){
	#Return from query cache here?
	#cache key = "$filepath:$key:$sr_offset:$byte_length"

	#define total offset

	#if(! exists $self->{file_cache}{$filepath}{off_sets}{$sr_key}){
	#  #warn "sr_key($sr_key) is not part of index for $filepath\n";
	#}
	#else{

	if(exists $self->{file_cache}{$filepath}{off_sets}{$sr_key}){

 	  my $total_offset = $self->{file_cache}{$filepath}{off_sets}{$sr_key} + $sr_offset;
	  my $seeked = sysseek($fh, $total_offset, 0);#0(whence) is SEEK_SET.

	  if(! $seeked){
		warn("Failed to seek to byte $total_offset in $filepath");
		#Don't undef fh here as this valid Slice maybe out of range
		#and we don't want to kill a valid fh
		#i.e. Slice start/end is past end of seq_region
	  }
	  else{
		my $read_bytes = sysread($fh, $blob_substr, $byte_length);
		
		if(! ((defined $read_bytes) && ($read_bytes == $byte_length))){
		  #! defined is error 0 is end of file
		  warn "Failed to read from $filepath\n$!";

		  if($read_bytes == 0){
			#This maybe because the slice is out of range!
			#The API gives no warning about this
						
			warn "End Of File encountered\n";
			warn "Total offset:\t".$self->{file_cache}{$filepath}{off_sets}{$sr_key}.
			  "  key($sr_key)  + $sr_offset = $total_offset\n";

			#add some checks against the theoretical/true length of the file?
		  }
		  else{  #Delete fh as it is useless/unsafe to retry
			undef $self->{file_cache}{$filepath}{filehandle};
			#$blob_substr is now set to empty string by read
			undef $blob_substr;
		  }
		}		
	  }
	}	
  }

  return $blob_substr;
}


1;
