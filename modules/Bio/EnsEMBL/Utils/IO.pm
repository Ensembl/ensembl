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

package Bio::EnsEMBL::Utils::IO;

=pod


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=pod

=head1 NAME

Bio::EnsEMBL::Utils::IO

=head1 SYNOPSIS

	use Bio::EnsEMBL::Utils::IO qw/slurp work_with_file slurp_to_array fh_to_array/;
	#or
	# use Bio::EnsEMBL::Utils::IO qw/:slurp/; # brings in any method starting with slurp
	# use Bio::EnsEMBL::Utils::IO qw/:array/; # brings in any method which ends with _array
	# use Bio::EnsEMBL::Utils::IO qw/:gz/;    # brings all methods which start with gz_
	# use Bio::EnsEMBL::Utils::IO qw/:bz/;    # brings all methods which start with bz_
	# use Bio::EnsEMBL::Utils::IO qw/:zip/;   # brings all methods which start with zip_
	# use Bio::EnsEMBL::Utils::IO qw/:all/;   # brings all methods in
	
  # As a scalar
  my $file_contents = slurp('/my/file/location.txt');
  print length($file_contents);
  
  # As a ref
  my $file_contents_ref = slurp('/my/file/location.txt', 1);
  print length($$file_contents_ref);
  
  # Sending it to an array
  my $array = slurp_to_array('/my/location');
  work_with_file('/my/location', 'r', sub {
    $array = process_to_array($_[0], sub {
      #Gives us input line by line
      return "INPUT: $_";
    });
  });
  
  # Simplified vesion but without the post processing
  $array = fh_to_array($fh);
  
  # Sending this back out to another file
  work_with_file('/my/file/newlocation.txt', 'w', sub {
    my ($fh) = @_;
    print $fh $$file_contents_ref;
    return;
  });
  
  # Gzipping the data to another file
  gz_work_with_file('/my/file.gz', 'w', sub {
    my ($fh) = @_;
    print $fh $$file_contents_ref;
    return;
  });
  
  # Working with a set of lines manually
  work_with_file('/my/file', 'r', sub {
    my ($fh) = @_;
    iterate_lines($fh, sub {
      my ($line) = @_;
      print $line; #Send the line in the file back out
      return;
    });
    return;
  });
  
  # Doing the same in one go
  iterate_file('/my/file', sub {
    my ($line) = @_;
    print $line; #Send the line in the file back out
    return;
  });
  
  # Move all data from one file handle to another. Bit like a copy
  move_data($src_fh, $trg_fh);
  	
=head1 DESCRIPTION

A collection of subroutines aimed to helping IO based operations

=head1 METHODS

See subroutines.

=head1 MAINTAINER

$Author$

=head1 VERSION

$Revision$

=cut

use strict;
use warnings;

use base qw(Exporter);

our $GZIP_OK = 0;
our $BZIP2_OK = 0;
our $ZIP_OK = 0;

our @EXPORT_OK = qw/slurp slurp_to_array fh_to_array process_to_array work_with_file gz_slurp gz_slurp_to_array gz_work_with_file bz_slurp bz_slurp_to_array bz_work_with_file zip_slurp zip_slurp_to_array zip_work_with_file spurt filter_dir iterate_file iterate_lines move_data/;
our %EXPORT_TAGS = (
  all     => [@EXPORT_OK],
  slurp   => [qw/slurp slurp_to_array gz_slurp gz_slurp_to_array/],
  spurt   => [qw/spurt/],
  array   => [qw/fh_to_array process_to_array slurp_to_array gz_slurp_to_array/],
  gz      => [qw/gz_slurp gz_slurp_to_array gz_work_with_file/],
  bz      => [qw/bz_slurp bz_slurp_to_array bz_work_with_file/],
  zip     => [qw/zip_slurp zip_slurp_to_array zip_work_with_file/],
  iterate => [qw/iterate_file iterate_lines/],
);
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Utils::Scalar qw(:assert);
use IO::File;

eval {
  require IO::Compress::Gzip;
  require IO::Uncompress::Gunzip;
  $GZIP_OK = 1;
};

eval {
  require IO::Compress::Bzip2;
  require IO::Uncompress::Bunzip2;
  $BZIP2_OK = 1;
};

eval {
  require IO::Compress::Zip;
  require IO::Uncompress::Unzip;
  $ZIP_OK = 1;
};

=head2 slurp()

  Arg [1]     : string $file
  Arg [2]     : boolean; $want_ref
  Arg [3]     : boolean; $binary
                Indicates if we want to return a scalar reference
  Description : Forces the contents of a file into a scalar. This is the 
                fastest way to get a file into memory in Perl. You can also
                get a scalar reference back to avoid copying the file contents
                in Scalar references. If the input file is binary then specify
                with the binary flag
  Returntype  : Scalar or reference of the file contents depending on arg 2
  Example     : my $contents = slurp('/tmp/file.txt');
  Exceptions  : If the file did not exist or was not readable
  Status      : Stable

=cut

sub slurp {
	my ($file, $want_ref, $binary) = @_;
	my $contents = q{};
	work_with_file($file, 'r', sub {
	  my ($fh) = @_;
	  binmode($fh) if $binary;
    my $size_left = -s $file;
    while( $size_left > 0 ) {
      my $read_cnt = sysread($fh, $contents, $size_left, length($contents));
      unless( $read_cnt ) {
        throw "read error in file $file: $!" ;
        last;
      }
      $size_left -= $read_cnt ;
    }
	  return;
	});
	return ($want_ref) ? \$contents : $contents;
}

=head2 spurt()

  Arg [1]     : string $file
  Arg [2]     : string $contents
  Arg [3]     : boolean; $append
  Arg [4]     : boolean; $binary
  Description : Convenient method to safely open a file and dump some content into it.
                $append can be set to append to the file instead of resetting it first.
                $binary can be set if the content you are printing is not plain-text.
  Returntype  : None
  Example     : spurt('/tmp/file.txt', $contents);
  Exceptions  : If the file could not be created or was not writable
  Status      : Stable

=cut

sub spurt {
  my ( $file, $contents, $append, $binary ) = @_;
  work_with_file(
    $file,
    $append ? 'a' : 'w',
    sub {
      my ($fh) = @_;
      binmode($fh) if $binary;
      syswrite( $fh, $contents );
    } );
}

=head2 gz_slurp

  Arg [1]     : string $file
  Arg [2]     : boolean; $want_ref Indicates if we want to return a scalar reference
  Arg [3]     : boolean; $binary
  Arg [4]     : HashRef arguments to pass into IO compression layers
  Description : Forces the contents of a file into a scalar. This is the 
                fastest way to get a file into memory in Perl. You can also
                get a scalar reference back to avoid copying the file contents
                in Scalar references. If the input file is binary then specify
                with the binary flag
  Returntype  : Scalar or reference of the file contents depending on arg 2
  Example     : my $contents = slurp('/tmp/file.txt.gz');
  Exceptions  : If the file did not exist or was not readable
  Status      : Stable

=cut

sub gz_slurp {
  my ($file, $want_ref, $binary, $args) = @_;
  my $contents;
  gz_work_with_file($file, 'r', sub {
    my ($fh) = @_;
    local $/ = undef;
    binmode($fh) if $binary;
    $contents = <$fh>;
    return;
  }, $args);
  return ($want_ref) ? \$contents : $contents;
}

=head2 bz_slurp

  Arg [1]     : string $file
  Arg [2]     : boolean; $want_ref Indicates if we want to return a scalar reference
  Arg [3]     : boolean; $binary
  Arg [4]     : HashRef arguments to pass into IO compression layers
  Description : Forces the contents of a file into a scalar. This is the 
                fastest way to get a file into memory in Perl. You can also
                get a scalar reference back to avoid copying the file contents
                in Scalar references. If the input file is binary then specify
                with the binary flag
  Returntype  : Scalar or reference of the file contents depending on arg 2
  Example     : my $contents = slurp('/tmp/file.txt.bz2');
  Exceptions  : If the file did not exist or was not readable
  Status      : Stable

=cut

sub bz_slurp {
  my ($file, $want_ref, $binary, $args) = @_;
  my $contents;
  bz_work_with_file($file, 'r', sub {
    my ($fh) = @_;
    local $/ = undef;
    binmode($fh) if $binary;
    $contents = <$fh>;
    return;
  }, $args);
  return ($want_ref) ? \$contents : $contents;
}

=head2 zip_slurp

  Arg [1]     : string $file
  Arg [2]     : boolean; $want_ref Indicates if we want to return a scalar reference
  Arg [3]     : boolean; $binary
  Arg [4]     : HashRef arguments to pass into IO compression layers
  Description : Forces the contents of a file into a scalar. This is the 
                fastest way to get a file into memory in Perl. You can also
                get a scalar reference back to avoid copying the file contents
                in Scalar references. If the input file is binary then specify
                with the binary flag
  Returntype  : Scalar or reference of the file contents depending on arg 2
  Example     : my $contents = slurp('/tmp/file.txt.zip');
  Exceptions  : If the file did not exist or was not readable
  Status      : Stable

=cut

sub zip_slurp {
  my ($file, $want_ref, $binary, $args) = @_;
  my $contents;
  zip_work_with_file($file, 'r', sub {
    my ($fh) = @_;
    local $/ = undef;
    binmode($fh) if $binary;
    $contents = <$fh>;
    return;
  }, $args);
  return ($want_ref) ? \$contents : $contents;
}


=head2 slurp_to_array

  Arg [1]     : string $file
  Arg [2]     : boolean $chomp
  Description : Sends the contents of the given file into an ArrayRef
  Returntype  : ArrayRef
  Example     : my $contents_array = slurp_to_array('/tmp/file.txt');
  Exceptions  : If the file did not exist or was not readable
  Status      : Stable

=cut

sub slurp_to_array {
  my ($file, $chomp) = @_;
  my $contents;
  work_with_file($file, 'r', sub {
	  my ($fh) = @_;
	  $contents = fh_to_array($fh, $chomp);
	  return;
	});
	return $contents;
}

=head2 gz_slurp_to_array

  Arg [1]     : string $file
  Arg [2]     : boolean $chomp
  Arg [3]     : HashRef arguments to pass into IO compression layers
  Description : Sends the contents of the given gzipped file into an ArrayRef
  Returntype  : ArrayRef
  Example     : my $contents_array = gz_slurp_to_array('/tmp/file.txt.gz');
  Exceptions  : If the file did not exist or was not readable
  Status      : Stable

=cut

sub gz_slurp_to_array {
  my ($file, $chomp, $args) = @_;
  my $contents;
  gz_work_with_file($file, 'r', sub {
    my ($fh) = @_;
    $contents = fh_to_array($fh, $chomp);
    return;
  }, $args);
  return $contents;
}

=head2 bz_slurp_to_array

  Arg [1]     : string $file
  Arg [2]     : boolean $chomp
  Arg [3]     : HashRef arguments to pass into IO compression layers
  Description : Sends the contents of the given bzipped file into an ArrayRef
  Returntype  : ArrayRef
  Example     : my $contents_array = bz_slurp_to_array('/tmp/file.txt.bz2');
  Exceptions  : If the file did not exist or was not readable
  Status      : Stable

=cut

sub bz_slurp_to_array {
  my ($file, $chomp, $args) = @_;
  my $contents;
  bz_work_with_file($file, 'r', sub {
    my ($fh) = @_;
    $contents = fh_to_array($fh, $chomp);
    return;
  }, $args);
  return $contents;
}

=head2 zip_slurp_to_array

  Arg [1]     : string $file
  Arg [2]     : boolean $chomp
  Arg [3]     : HashRef arguments to pass into IO compression layers
  Description : Sends the contents of the given zipped file into an ArrayRef
  Returntype  : ArrayRef
  Example     : my $contents_array = zip_slurp_to_array('/tmp/file.txt.zip');
  Exceptions  : If the file did not exist or was not readable
  Status      : Stable

=cut

sub zip_slurp_to_array {
  my ($file, $chomp, $args) = @_;
  my $contents;
  zip_work_with_file($file, 'r', sub {
    my ($fh) = @_;
    $contents = fh_to_array($fh, $chomp);
    return;
  }, $args);
  return $contents;
}

=head2 fh_to_array

  Arg [1]     : Glob/IO::Handle $fh
  Arg [2]     : boolean $chomp
  Description : Sends the contents of the given filehandle into an ArrayRef. 
                Will perform chomp on each line if specified. If you require
                any more advanced line based processing then see 
                L<process_to_array>.
  Returntype  : ArrayRef
  Example     : my $contents_array = fh_to_array($fh);
  Exceptions  : None
  Status      : Stable

=cut

sub fh_to_array {
  my ($fh, $chomp) = @_;
  if($chomp) {
    return process_to_array($fh, sub {
      my ($line) = @_;
      chomp($line);
      return $line;
    });
  }
  my @contents = <$fh>;
  return \@contents;
}

=head2 process_to_array

  Arg [1]     : Glob/IO::Handle $fh
  Arg [2]     : CodeRef $callback
  Description : Sends the contents of the given file handle into an ArrayRef
                via the processing callback. Assumes line based input.
  Returntype  : ArrayRef
  Example     : my $array = process_to_array($fh, sub { return "INPUT: $_"; });
  Exceptions  : If the fh did not exist or if a callback was not given.
  Status      : Stable

=cut

sub process_to_array {
  my ($fh, $callback) = @_;
  assert_file_handle($fh, 'FileHandle');
  assert_ref($callback, 'CODE', 'callback');
  my @contents;
  iterate_lines($fh, sub {
    my ($line) = @_;
    push(@contents, $callback->($line));
    return;
  });
  return \@contents;
}

=head2 iterate_lines

  Arg [1]     : Glob/IO::Handle $fh
  Arg [2]     : CodeRef $callback
  Description : Iterates through each line from the given file handle and
                hands them to the callback one by one
  Returntype  : None
  Example     : iterate_lines($fh, sub { print "INPUT: $_"; });
  Exceptions  : If the fh did not exist or if a callback was not given.
  Status      : Stable

=cut

sub iterate_lines {
  my ($fh, $callback) = @_;
  assert_file_handle($fh, 'FileHandle');
  assert_ref($callback, 'CODE', 'callback');
  while( my $line = <$fh> ) {
    $callback->($line);
  }
  return;
}

=head2 iterate_file

  Arg [1]     : string $file
  Arg [3]     : CodeRef the callback which is used to iterate the lines in
                the file
  Description : Iterates through each line from the given file and
                hands them to the callback one by one
  Returntype  : None
  Example     : iterate_file('/my/file', sub { print "INPUT: $_"; });
  Exceptions  : If the file did not exist or if a callback was not given.
  Status      : Stable

=cut


sub iterate_file {
  my ($file, $callback) = @_;
  work_with_file($file, 'r', sub {
    my ($fh) = @_;
    iterate_lines($fh, $callback);
    return;
  });
  return;
}



=head2 work_with_file

  Arg [1]     : string $file
  Arg [2]     : string; $mode 
                Supports all modes specified by the C<open()> function as well as those 
                supported by IO::File
  Arg [3]     : CodeRef the callback which is given the open file handle as
                its only argument
  Description : Performs the nitty gritty of checking if a file handle is open
                and closing the resulting filehandle down.
  Returntype  : None
  Example     : work_with_file('/tmp/out.txt', 'w', sub { 
                  my ($fh) = @_; 
                  print $fh 'hello'; 
                  return;
                });
  Exceptions  : If we could not work with the file due to permissions
  Status      : Stable

=cut

sub work_with_file {
  my ($file, $mode, $callback) = @_;
  throw "We need a file name to open" if ! $file;
  throw "We need a mode to open the requested file with" if ! $mode;
  assert_ref($callback, 'CODE', 'callback');
  my $fh = IO::File->new($file, $mode) or
    throw "Cannot open '${file}' in  mode '${mode}': $!";
  $callback->($fh);
  close($fh) or throw "Cannot close FH from ${file}: $!";
  return;
}

=head2 gz_work_with_file

  Arg [1]     : string $file
  Arg [2]     : string; $mode 
                Supports modes like C<r>, C<w>, C<\>> and C<\<>
  Arg [3]     : CodeRef the callback which is given the open file handle as
                its only argument
  Arg [4]     : HashRef used to pass options into the IO 
                compression/uncompression modules
  Description : Performs the nitty gritty of checking if a file handle is open
                and closing the resulting filehandle down.
  Returntype  : None
  Example     : gz_work_with_file('/tmp/out.txt.gz', 'w', sub { 
                  my ($fh) = @_; 
                  print $fh 'hello'; 
                  return;
                });
  Exceptions  : If we could not work with the file due to permissions
  Status      : Stable

=cut

sub gz_work_with_file {
  my ($file, $mode, $callback, $args) = @_;
  throw "IO::Compress was not available" if ! $GZIP_OK;
  throw "We need a file name to open" if ! $file;
  throw "We need a mode to open the requested file with" if ! $mode;
  assert_ref($callback, 'CODE', 'callback');
  $args ||= {};
  
  my $fh;
  {
    no warnings qw/once/;
    if($mode =~ '>$' || $mode eq 'w') {
      $args->{Append} = 1 if $mode =~ />>$/;
      $fh = IO::Compress::Gzip->new($file, %$args) or throw "Cannot open '$file' for writing: $IO::Compress::Gzip::GzipError";
    }
    elsif($mode eq '<' || $mode eq 'r') {
      $fh = IO::Uncompress::Gunzip->new($file, %$args) or throw "Cannot open '$file' for writing: $IO::Uncompress::Gunzip::GunzipError";
    }
    else {
      throw "Could not decipher a mode from '$mode'";
    }
  };
  $callback->($fh);
  close($fh) or throw "Cannot close FH from ${file}: $!";
  return;
}

=head2 bz_work_with_file

  Arg [1]     : string $file
  Arg [2]     : string; $mode 
                Supports modes like C<r>, C<w>, C<\>> and C<\<>
  Arg [3]     : CodeRef the callback which is given the open file handle as
                its only argument
  Arg [4]     : HashRef used to pass options into the IO 
                compression/uncompression modules
  Description : Performs the nitty gritty of checking if a file handle is open
                and closing the resulting filehandle down.
  Returntype  : None
  Example     : bz_work_with_file('/tmp/out.txt.bz2', 'w', sub { 
                  my ($fh) = @_; 
                  print $fh 'hello'; 
                  return;
                });
  Exceptions  : If we could not work with the file due to permissions
  Status      : Stable

=cut

sub bz_work_with_file {
  my ($file, $mode, $callback, $args) = @_;
  throw "IO::Compress was not available" if ! $BZIP2_OK;
  throw "We need a file name to open" if ! $file;
  throw "We need a mode to open the requested file with" if ! $mode;
  assert_ref($callback, 'CODE', 'callback');
  $args ||= {};
  
  my $fh;
  {
    no warnings qw/once/;
    if($mode =~ '>$' || $mode eq 'w') {
      $args->{Append} = 1 if $mode =~ />>$/;
      $fh = IO::Compress::Bzip2->new($file, %$args) or throw "Cannot open '$file' for writing: $IO::Compress::Bzip2::Bzip2Error";
    }
    elsif($mode eq '<' || $mode eq 'r') {
      $fh = IO::Uncompress::Bunzip2->new($file, %$args) or throw "Cannot open '$file' for writing: $IO::Uncompress::Bunzip2::Bunzip2Error";
    }
    else {
      throw "Could not decipher a mode from '$mode'";
    }
  };
  $callback->($fh);
  close($fh) or throw "Cannot close FH from ${file}: $!";
  return;
}

=head2 zip_work_with_file

  Arg [1]     : string $file
  Arg [2]     : string; $mode 
                Supports modes like C<r>, C<w>, C<\>> and C<\<>
  Arg [3]     : CodeRef the callback which is given the open file handle as
                its only argument
  Arg [4]     : HashRef used to pass options into the IO 
                compression/uncompression modules
  Description : Performs the nitty gritty of checking if a file handle is open
                and closing the resulting filehandle down.
  Returntype  : None
  Example     : zip_work_with_file('/tmp/out.txt.zip', 'w', sub { 
                  my ($fh) = @_; 
                  print $fh 'hello'; 
                  return;
                });
  Exceptions  : If we could not work with the file due to permissions
  Status      : Stable

=cut

sub zip_work_with_file {
  my ($file, $mode, $callback, $args) = @_;
  throw "IO::Compress was not available" if ! $ZIP_OK;
  throw "We need a file name to open" if ! $file;
  throw "We need a mode to open the requested file with" if ! $mode;
  assert_ref($callback, 'CODE', 'callback');
  $args ||= {};
  
  my $fh;
  {
    no warnings qw/once/;
    if($mode =~ '>$' || $mode eq 'w') {
      $args->{Append} = 1 if $mode =~ />>$/;
      $fh = IO::Compress::Zip->new($file, %$args) or throw "Cannot open '$file' for writing: $IO::Compress::Zip::ZipError";
    }
    elsif($mode eq '<' || $mode eq 'r') {
      $fh = IO::Uncompress::Unzip->new($file, %$args) or throw "Cannot open '$file' for writing: $IO::Uncompress::Unzip::UnzipError";
    }
    else {
      throw "Could not decipher a mode from '$mode'";
    }
  };
  $callback->($fh);
  close($fh) or throw "Cannot close FH from ${file}: $!";
  return;
}

=head2 filter_dir

  Arg [1]     : String; directory
  Arg [2]     : CodeRef; the callback which is given a file in the
                directory as its only argument
  Description : Return the lexicographically sorted content of a directory.
                The callback allows to specify the criteria an entry in
                the directory must satisfy in order to appear in the content.
  Returntype  : Arrayref; list with the filtered files/directory
  Example     : filter_dir('/tmp', sub { 
                  my $file = shift;

                  # select perl scripts in the directory
                  return $file if $file =~ /\.pl$/; 
                });
  Exceptions  : If the directory cannot be opened or its handle
                cannot be closed
  Status      : Stable

=cut

sub filter_dir {
  my ($dir, $callback) = @_;

  assert_ref($callback, 'CODE', 'callback');

  opendir(my $dh, $dir) or throw "Cannot open directory $dir";
  my @files = sort grep { $callback->($_) } readdir($dh);
  closedir($dh) or throw "Cannot close directory $dir";

  return \@files;
}


=head2 move_data

  Arg [1]     : FileHandle $src_fh
  Arg [2]     : FileHandle $trg_fh
  Arg [3]     : int $buffer. Defaults to 8KB
  Description : Moves data from the given source filehandle to the target one
                using a 8KB buffer or user specified buffer
  Returntype  : None
  Example     : move_data($src_fh, $trg_fh, 16*1024); # copy in 16KB chunks
  Exceptions  : If inputs were not as expected

=cut

sub move_data {
  my ($src_fh, $trg_fh, $buffer_size) = @_;
  assert_file_handle($src_fh, 'SourceFileHandle');
  assert_file_handle($trg_fh, 'TargetFileHandle');
  
  $buffer_size ||= 8192; #Default 8KB
  my $buffer;
  while(1) {
    my $read = sysread($src_fh, $buffer, $buffer_size);
    if(! defined $read) {
      throw "Error whilst reading from filehandle: $!";
    }
    if($read == 0) {
      last;
    }
    my $written = syswrite($trg_fh, $buffer);
    if(!defined $written) {
      throw "Error whilst writing to filehandle: $!";
    }
  }
  return;
}

1;
