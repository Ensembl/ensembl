package XrefParser::FetchFiles;

use strict;
use warnings;

# Given one or several FTP or HTTP URIs, download them.  If an URI is
# for a file or MySQL connection, then these will be ignored.  For
# FTP, standard shell file name globbing is allowed (but not regular
# expressions).  HTTP does not allow file name globbing.  The routine
# returns a list of successfully downloaded local files or an empty list
# if there was an error.


use Carp;
use DBI;
use Digest::MD5 qw(md5_hex);
use Getopt::Long;
use POSIX qw(strftime);

use File::Basename;
use File::Spec::Functions;
use IO::File;
use Net::FTP;
use URI;
use URI::file;
use Text::Glob qw( match_glob );
use LWP::UserAgent;


my $base_dir = File::Spec->curdir();

sub new {
    my ($proto) = @_;

    my $class = ref $proto || $proto;
    return bless {}, $class;
}

sub fetch_files {
  my ($self, $arg_ref) = @_;


  my $dest_dir         = $arg_ref->{dest_dir};
  my $user_uris        = $arg_ref->{user_uris};
  my $deletedownloaded = $arg_ref->{del_down};
  my $checkdownload    = $arg_ref->{chk_down};
  my $verbose          = $arg_ref->{verbose} ;

  my @processed_files;

  foreach my $user_uri (@$user_uris) {
    # Change old-style 'LOCAL:' URIs into 'file:'.
    $user_uri =~ s/^LOCAL:/file:/ix;
    my $uri = URI->new($user_uri);

    if ( $uri->scheme() eq 'script' ) {
      push( @processed_files, $user_uri );
    } elsif ( $uri->scheme() eq 'file' ) {

      # Deal with local files.

      $user_uri =~ s/file://x;
      if ( -s $user_uri ) {
        push( @processed_files, $user_uri );
      } else {
        printf( "==> Can not find file '%s' (or it is empty)\n",
                $user_uri );
        return ();
      }
    } elsif ( $uri->scheme() eq 'ftp' ) {
      # Deal with FTP files.

      my $file_path = catfile( $dest_dir, basename( $uri->path() ) );

      if ( $deletedownloaded && -e $file_path ) {
        if ($verbose) {
          printf( "Deleting '%s'\n", $file_path );
        }
        unlink($file_path);
      }

      if ( $checkdownload && -s $file_path ) {
        # The file is already there, no need to connect to a FTP
        # server.  This also means no file name globbing was
        # used (for globbing FTP URIs, we always need to connect
        # to a FTP site to see what files are there).

        if ($verbose) {
          printf( "File '%s' already exists\n", $file_path );
        }
        push( @processed_files, $file_path );
        next;
      }

      if ( -e $file_path ) { unlink($file_path) }

      if ($verbose) {
        printf( "Connecting to FTP host '%s' for file '%s' \n",
                $uri->host(), $file_path );
      }

      my $ftp = Net::FTP->new( $uri->host(), 'Debug' => 0);
      if ( !defined($ftp) ) {
        printf( "==> Can not open FTP connection: %s\n", $@ );
        return ();
      }

      if ( !$ftp->login( 'anonymous', '-anonymous@' ) ) {
        printf( "==> Can not log in on FTP host: %s\n",
                $ftp->message() );
        return ();
      }

      if ( !$ftp->cwd( dirname( $uri->path() ) ) ) {
        printf( "== Can not change directory to '%s': %s\n",
                dirname( $uri->path() ), $ftp->message() );
        return ();
      }

      $ftp->binary();

      if(!$ftp->ls()){
	$ftp->pasv(1);
      }
      foreach my $remote_file ( ( @{ $ftp->ls() } ) ) {
        if ( !match_glob( basename( $uri->path() ), $remote_file ) ) {
          next;
        }

        $file_path = catfile( $dest_dir, basename($remote_file) );

        if ( $deletedownloaded && -e $file_path ) {
          if ($verbose) {
            printf( "Deleting '%s'\n", $file_path );
          }
          unlink($file_path);
        }

        if ( $checkdownload && -s $file_path ) {
          if ($verbose) {
            printf( "File '%s' already exists\n", $file_path );
          }
        } else {

          if ( -e $file_path ) { unlink($file_path) }

          if ( !-d dirname($file_path) ) {
            if ($verbose) {
              printf( "Creating directory '%s'\n",
                      dirname($file_path) );
            }
            if ( !mkdir( dirname($file_path) ) ) {
              printf( "==> Can not create directory '%s': %s",
                      dirname($file_path), $! );
              return ();
            }
          }

          if ($verbose) {
            printf( "Fetching '%s' (size = %s)\n",
                    $remote_file,
                    $ftp->size($remote_file) || '(unknown)' );
            printf( "Local file is '%s'\n", $file_path );
          }

          if ( !$ftp->get( $remote_file, $file_path ) ) {
            printf( "==> Could not get '%s': %s\n",
                    basename( $uri->path() ), $ftp->message() );
            return ();
          }
        } ## end else [ if ( $checkdownload &&...)]

        if ( $file_path =~ /\.(gz|Z)$/x ) {
          # Read from zcat pipe
          #
          my $cmd = "gzip -t $file_path";
          if ( system($cmd) != 0 ) {
            printf( "system command '%s' failed: %s - "
                      . "Checking of gzip file failed - "
                      . "FILE CORRUPTED ?\n\n",
                    $cmd, $? );

            if ( -e $file_path ) {
              if ($verbose) {
                printf( "Deleting '%s'\n", $file_path );
              }
              unlink($file_path);
            }
            return ();
          } else {
            if ($verbose) {
              printf( "'%s' passed (gzip -t) corruption test.\n",
                      $file_path );
            }
          }
        }
        push( @processed_files, $file_path );

      } ## end foreach my $remote_file ( (...))

    } elsif ( $uri->scheme() eq 'http' ) {
      # Deal with HTTP files.

      my $file_path = catfile( $dest_dir, basename( $uri->path() ) );

      if ( $deletedownloaded && -e $file_path ) {
        if ($verbose) {
          printf( "Deleting '%s'\n", $file_path );
        }
        unlink($file_path);
      }

      if ( $checkdownload && -s $file_path ) {
        # The file is already there, no need to connect to a
        # HTTP server.

        if ($verbose) {
          printf( "File '%s' already exists\n", $file_path );
        }
        push( @processed_files, $file_path );
        next;
      }

      if ( -e $file_path ) { unlink($file_path) }

      if ( !-d dirname($file_path) ) {
        if ($verbose) {
          printf( "Creating directory '%s'\n", dirname($file_path) );
        }
        if ( !mkdir( dirname($file_path) ) ) {
          printf( "==> Can not create directory '%s': %s",
                  dirname($file_path), $! );
          return ();
        }
      }

      if ($verbose) {
        printf( "Connecting to HTTP host '%s'\n", $uri->host() );
        printf( "Fetching '%s'\n",                $uri->path() );
      }

      if ( $checkdownload && -s $file_path ) {
        if ($verbose) {
          printf( "File '%s' already exists\n", $file_path );
        }
      } else {

        if ($verbose) {
          printf( "Local file is '%s'\n", $file_path );
        }

        if ( -e $file_path ) { unlink($file_path) }

        my $ua = LWP::UserAgent->new();
        $ua->env_proxy();

        my $response =
          $ua->get( $uri->as_string(), ':content_file' => $file_path );

        if ( !$response->is_success() ) {
          printf( "==> Could not get '%s': %s\n",
                  basename( $uri->path() ), $response->content() );
          return ();
        }
      }

      push( @processed_files, $file_path );

    } elsif ( $uri->scheme() eq 'mysql' ) {
      # Just leave MySQL data untouched for now.
      push( @processed_files, $user_uri );
    } else {
      printf( "==> Unknown URI scheme '%s' in URI '%s'\n",
              $uri->scheme(), $uri->as_string() );
      return ();
    }
  } ## end foreach my $user_uri (@user_uris)

  return @processed_files;
} ## end sub fetch_files

1;
