#!/usr/local/ensembl/bin/perl -w

# this script will lokk through all listed databases for 
# database names which match given regexp.
# This version should be ok if one or more of the given database
# servers is not available.


use strict;
use DBI;

# seconds before we give up on database
my $timeout = 5;

my @all_locators = qw {
		     ecs2:3361:ensro::
		     ecs2:3362:ensro::
		     ecs2:3363:ensro::
		     ecs2:3364:ensro::
		     ecs2:3365:ensro::
		     ecs2:3366:ensro::

		     ecs4:3350:ensro::
		     ecs4:3351:ensro::
		     ecs4:3352:ensro::
		     ecs4:3353:ensro::

		     ecs1a:3306:ensro::
		     ecs1b:3306:ensro::
		     ecs1c:3306:ensro::
		     ecs1d:3306:ensro::
		     ecs1e:3306:ensro::
		     ecs1f:3306:ensro::
		     ecs1g:3306:ensro::
		     ecs1h:3306:ensro::

		     ecs3:3304:ensro::

		     ia64g:3306:ensro::
		     ia64e:3306:ensro::
		     ia64f:3306:ensro::
		      };


my $pattern=shift;
my $size = 0;
my %sizes;

if( @ARGV ) {
  $size = 1;
}

if( !$pattern ) {
  printf( "You need to supply a regexp as argument.\n" );
  printf( "Use -list if you want to see which databases are configured.\n" );
  printf( "If you have any parameters after the regexp, the program will print\n" );
  printf( " database sizes as well.\n" );
  printf( "\nExample: search_dbs.pl \"^stabenau\" -sizes | sort -knr3\n" ); 
  exit;
}

my $list = 0;

if( $pattern eq "-list" ) {
  $list = 1;
}


for my $loc ( @all_locators ) {
  my @elems = split( ":", $loc );
  my @dbnames = ();
  %sizes = ();
  my $dsn = sprintf( "DBI:mysql:host=%s;port=%s", $elems[0], $elems[1] );

  $SIG{ALRM} = sub{die("timeout");};

  eval {
    alarm $timeout;
    my $db = DBI->connect( $dsn, $elems[2], $elems[3], { RaiseError => 1 } );
    my $res = $db->selectall_arrayref( "show databases" );
    @dbnames = map { $_->[0] } @$res;
    $db->disconnect();
    alarm 0;
  };

  if( !@dbnames ) {
    print STDERR "$loc NOT OK\n";
  } else {
    if( $size ) {
      for my $dbname ( @dbnames ) {
	if( $dbname =~ /$pattern/ ) {
	  eval {
	    alarm $timeout;
	    my $db = DBI->connect( $dsn, $elems[2], $elems[3], { RaiseError => 1, PrintError=> 0 } );
	    $db->do( "use $dbname" );
	    my $t_status = $db->selectall_arrayref( "show table status" );
	    my $size = 0;
	    map { $size += $_->[6]; $size += $_->[8] } @$t_status;
	    print "$loc $dbname $size\n";
	    $db->disconnect();
	    alarm 0;
	  };
	  if( $@ ) {
	    print( "Problem on $loc $dbname.\n ", $@ );
	  }
	}
      }
    } else {
      if( $list ) {
	print STDERR "$loc ok\n";
      } else {
	for my $dbname ( @dbnames ) {
	  if( $dbname =~ /$pattern/ ) {
	    print "$loc $dbname\n";
	  }
	}
      }
    }
  }
}


