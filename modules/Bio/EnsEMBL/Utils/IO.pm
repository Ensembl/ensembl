package Bio::EnsEMBL::Utils::IO;

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

=cut

=pod

=head1 NAME

Bio::EnsEMBL::Utils::IO

=head1 SYNOPSIS

	use Bio::EnsEMBL::Utils::IO qw/slurp work_with_file slurp_to_array fh_to_array/;
	#or
	# use Bio::EnsEMBL::Utils::IO qw/:slurp/; #brings in any method starting with slurp
	# use Bio::EnsEMBL::Utils::IO qw/:array/; #brings in any method which ends with _array
	# use Bio::EnsEMBL::Utils::IO qw/:gz/;    #brings all methods which start with gz_
	# use Bio::EnsEMBL::Utils::IO qw/:all/;   #brings all methods in
	
	#As a scalar
  my $file_contents = slurp('/my/file/location.txt');
  print length($file_contents);
  
  #As a ref
  my $file_contents_ref = slurp('/my/file/location.txt', 1);
  print length($$file_contents_ref);
  
  #Sending it to an array
  my $array = slurp_to_array('/my/location');
  work_with_file('/my/location', 'r', sub {
    $array = process_to_array($_[0], sub {
      #Gives us input line by line
      return "INPUT: $_";
    });
  });
  
  #Simplified vesion but without the post processing
  $array = fh_to_array($fh);
  
  #Sending this back out to another file
  work_with_file('/my/file/newlocation.txt', 'w', sub {
    my ($fh) = @_;
    print $fh $$file_contents_ref;
    return;
  });
  
  #Gzipping the data to another file
  gz_work_with_file('/my/file.gz', 'w', sub {
    my ($fh) = @_;
    print $fh $$file_contents_ref;
    return;
  });
  
  #Working with a set of lines manually
  work_with_file('/my/file', 'r', sub {
    my ($fh) = @_;
    iterate_lines($fh, sub {
      my ($line) = @_;
      print $line; #Send the line in the file back out
      return;
    });
    return;
  });
  
  #Doing the same in one go
  iterate_file('/my/file', sub {
    my ($line) = @_;
    print $line; #Send the line in the file back out
    return;
  });
  
  #Move all data from one file handle to another. Bit like a copy
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
our @EXPORT_OK = qw/slurp slurp_to_array fh_to_array process_to_array work_with_file gz_slurp gz_slurp_to_array gz_work_with_file iterate_file iterate_lines move_data/;
our %EXPORT_TAGS = (
  all     => [@EXPORT_OK],
  slurp   => [qw/slurp slurp_to_array gz_slurp gz_slurp_to_array/],
  array   => [qw/fh_to_array process_to_array slurp_to_array gz_slurp_to_array/],
  gz      => [qw/gz_slurp gz_slurp_to_array gz_work_with_file/],
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
    my $size_left = -s $fh;
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

=head2 gz_slurp()

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
  Example     : my $contents = slurp('/tmp/file.txt.gz');
  Exceptions  : If the file did not exist or was not readable
  Status      : Stable

=cut

sub gz_slurp {
  my ($file, $want_ref, $binary) = @_;
  my $contents;
  gz_work_with_file($file, 'r', sub {
    my ($fh) = @_;
    local $/ = undef;
    binmode($fh) if $binary;
    $contents = <$fh>;
    return;
  });
  return ($want_ref) ? \$contents : $contents;
}

=head2 slurp_to_array()

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

=head2 gz_slurp_to_array()

  Arg [1]     : string $file
  Arg [2]     : boolean $chomp
  Description : Sends the contents of the given gzipped file into an ArrayRef
  Returntype  : ArrayRef
  Example     : my $contents_array = slurp_to_array('/tmp/file.txt.gz');
  Exceptions  : If the file did not exist or was not readable
  Status      : Stable

=cut

sub gz_slurp_to_array {
  my ($file, $chomp) = @_;
  my $contents;
  gz_work_with_file($file, 'r', sub {
    my ($fh) = @_;
    $contents = fh_to_array($fh, $chomp);
    return;
  });
  return $contents;
}

=head2 fh_to_array()

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



=head2 work_with_file()

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

=head2 gz_work_with_file()

  Arg [1]     : string $file
  Arg [2]     : string; $mode 
                Supports modes like C<r>, C<w>, C<\>> and C<\<>
  Arg [3]     : CodeRef the callback which is given the open file handle as
                its only argument
  Description : Performs the nitty gritty of checking if a file handle is open
                and closing the resulting filehandle down.
  Returntype  : None
  Example     : work_with_file('/tmp/out.txt.gz', 'w', sub { 
                  my ($fh) = @_; 
                  print $fh 'hello'; 
                  return;
                });
  Exceptions  : If we could not work with the file due to permissions
  Status      : Stable

=cut

sub gz_work_with_file {
  my ($file, $mode, $callback) = @_;
  throw "IO::Compress was not available"if ! $GZIP_OK;
  throw "We need a file name to open" if ! $file;
  throw "We need a mode to open the requested file with" if ! $mode;
  assert_ref($callback, 'CODE', 'callback');
  my $fh;
  {
    no warnings qw/once/;
    if($mode =~ '>$' || $mode eq 'w') {
      my $append = ($mode =~ />>$/) ? 1 : 0;
      $fh = IO::Compress::Gzip->new($file, Append => $append) or throw "Cannot open '$file' for writing: $IO::Compress::Gzip::GzipError";
    }
    elsif($mode eq '<' || $mode eq 'r') {
      $fh = IO::Uncompress::Gunzip->new($file) or throw "Cannot open '$file' for writing: $IO::Uncompress::Gunzip::GunzipError";
    }
    else {
      throw "Could not decipher a mode from '$mode'";
    }
  };
  $callback->($fh);
  close($fh) or throw "Cannot close FH from ${file}: $!";
  return;
  return;
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
