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

# Module to ease mysql schema conversion
# Usage:

# Make a schema_converter with new ( source_dbh, target_dbh )
# source database should be filled, target is empty schema


# For each target table there will be a transfer
#  Either with a self specified transfer function
#  with a custom select
#  from a renamed source table
#  from the same name source table   

# configure the transfer (and the order) with

# table_rename( "oldname", "newname" ) 
# table_skip( "tablename" )
# do_first( "newtable1", "newtable2", "newtable3" )


# Each standard table transfer 
#  Either do a custom select statement
#  or transfer columns with same name or renamed into each other
#  specify columns to omit in target or get error

# column_rename( "tablename", "oldcolname", "newcolname" )
# column_skip( $targetdb, "table", "column" ) 
# custom_select( $targetdb, "tablename", 

# Each row may be modified (custom select or standard select) 
#  specify a row_modifier function for the target table
#  It takes a list ref and returns a list ref with the modified values
#  ( you have to know the order of columns which come in have to go out db )

# set_row_modifier( "tablename", function_reference ) 
#
# potentially large (target) tables should be declared with
# big_table( "newtablename" );

package SchemaConverter;

use strict;
use DBI;
use Data::Dumper;


sub new {
  my ( $class, @args ) = @_;
  
  my $self = {};
  bless $self, $class;

  $self->source_dbh( $args[0] );
  $self->target_dbh( $args[1] );

  $self->read_dbs();

  return $self;
}

sub tmp_dir {
  my ( $self, $arg ) = @_;
  
  ( defined $arg ) &&
    ( $self->{'tmp_dir'} = $arg );

  return $self->{'tmp_dir'};
}

sub source_dbh {
  my ( $self, $arg ) = @_;
  
  ( defined $arg ) &&
    ( $self->{'source_dbh'} = $arg );

  return $self->{'source_dbh'};
}

sub target_dbh {
  my ( $self, $arg ) = @_;
  
  ( defined $arg ) &&
    ( $self->{'target_dbh'} = $arg );

  return $self->{'target_dbh'};
}

sub close_dbh {
  my $self = shift;

  $self->source_dbh()->disconnect();
  $self->target_dbh()->disconnect();
}


sub transfer {
  my $self = shift;
  
  local *FH;
  my $tmpdir;

  if( ! defined $self->tmp_dir() ) {
    $self->close_dbh();
    die( "No tmp_dir specified" );
  } else {
    $tmpdir = $self->tmp_dir();
  }

  # first we should check if all standard trnasfers can go
  # otherwise testing this is a pain
  $self->check_possible_transfer();
  my @ordered_tables = $self->{targetdb}{dofirst};
  my %all_tables = %{$self->{targetdb}{tables}};
  for my $tablename ( @ordered_tables ) {
    delete $all_tables{$tablename};
  }
  for my $tablename ( keys %all_tables ) {
    push( @ordered_tables, $tablename );
  }

  for my $tablename ( @ordered_tables ) {
    my $skip = 0;
    print STDERR "Transfer $tablename ";

    open( FH, ">$tmpdir/$tablename.txt" ) or die "cant open dumpfile";
    
    if( exists $self->{targetdb}{tables}{$tablename}{transfer} ) {
      my $transfunc = $self->{targetdb}{tables}{$tablename}{transfer};
      &$transfunc( $self->source_dbh(), $self->target_dbh(), $tablename, \*FH );
    } else {
      my $sourcetable;
      if( exists $self->{targetdb}{tables}{$tablename}{link} ) {
	$sourcetable = $self->{targetdb}{tables}{$tablename}{link};
	if( $sourcetable eq "" ) {
	  # skip this table
	  $skip = 1;
	}
      } else {
	# find the sourcetable
	if( exists $self->{targetdb}{tables}{$tablename}{select} ) {

	  # if we have custom select, sourcetable doesnt make sense
	  $sourcetable = undef;

	} elsif( ! exists $self->{sourcedb}{tables}{$tablename} ) {

	  die "Couldnt find source for $tablename. Enter empty sourcetable.";

	} else {
	  $sourcetable = $tablename;
	}
      }
      if( ! $skip ) {
	$self->standard_table_transfer( $sourcetable, $tablename, \*FH );
      }
    }

    close FH;

    if( ! $skip ) {
      $self->target_dbh->do( "load data infile '$tmpdir/$tablename.txt' into table $tablename" );
    }

    # upload ?
    unlink "$tmpdir/$tablename.txt";
    print STDERR " finished\n";
  }

  # close databases
  
}


# this function checks if a standard transfer is possible.
# custom selects and custom transfer functions are not checked
# it will only return if it can do the transfer otherwise die

sub check_possible_transfer {
  my $self = shift;

  for my $tablename ( keys %{$self->{targetdb}{tables}} ) {
    my $sourcetable;

    if( exists $self->{targetdb}{tables}{$tablename}{transfer} ) {
      next;
    }

    if( exists $self->{targetdb}{tables}{$tablename}{link} ) {
      $sourcetable = $self->{targetdb}{tables}{$tablename}{link};
      if( $sourcetable eq "" ) {
	# skip this table
	next;
      }
    } else {
      # find the sourcetable
      if( exists $self->{targetdb}{tables}{$tablename}{select} ) {

	# custom select, no check
	next;
      } elsif( ! exists $self->{sourcedb}{tables}{$tablename} ) {

	die "Couldnt find source for $tablename. Enter empty sourcetable.";

      } else {
	$sourcetable = $tablename;
      }
    }

    my @newcols = @{$self->{targetdb}{tables}{$tablename}{columns}};
    my @oldcols = @{$self->{sourcedb}{tables}{$sourcetable}{columns}};
    
    my %rename;

    if( exists $self->{targetdb}{tables}{$tablename}{columnrename} ) {
      %rename = %{$self->{targetdb}{tables}{$tablename}{columnrename}};
    } else {
      %rename = ();
    }

    # find all source columns and build select statement 

    for my $colname ( @newcols ) {
      my $selname;
      
      if( exists $rename{$colname} ) {
	$selname = $rename{$colname};
	if( $selname eq "" ) {
	  $selname = "NULL";
	}
      } else {
	my $colExists = 0;
	for my $oldcol ( @oldcols ) {
	  if(  $oldcol eq $colname ) {
	    $selname = $colname;
	    $colExists = 1;
	    last;
	  }
	}
	if( ! $colExists ) {
	  die "Couldnt fill $tablename.$colname\n";
	}
      }
    }
  }
}


sub standard_table_transfer {
  my ( $self, $sourcetable, $targettable, $tmpfile ) = @_;
  
  my $sourcedb = $self->source_dbh();
  my $targetdb = $self->target_dbh();

  # look for custom select
  my $select = "";
  if( exists $self->{targetdb}{tables}{$targettable}{select} ) {
    $select = $self->{targetdb}{tables}{$targettable}{select};
  } else {
    # check if all columns have matching names
    my @newcols = @{$self->{targetdb}{tables}{$targettable}{columns}};
    my @oldcols = @{$self->{sourcedb}{tables}{$sourcetable}{columns}};
    
    my %rename;

    if( exists $self->{targetdb}{tables}{$targettable}{columnrename} ) {
      %rename = %{$self->{targetdb}{tables}{$targettable}{columnrename}};
    } else {
      %rename = ();
    }

    # find all source columns and build select statement 

    for my $colname ( @newcols ) {
      my $selname;

      if( exists $rename{$colname} ) {
	$selname = $rename{$colname};
	if( $selname eq "" ) {
	  $selname = "NULL";
	}
      } else {
	my $colExists = 0;
	for my $oldcol ( @oldcols ) {
	  if(  $oldcol eq $colname ) {
	    $selname = $colname;
	    $colExists = 1;
	    last;
	  }
	}
	if( ! $colExists ) {
	  die "Couldnt fill $targettable.$colname\n";
	}
      }

      $select .= " $selname,";
    }
    chop( $select );
    $select = "SELECT $select from $sourcetable";
  }    
		   
  # MySQL specific ...
  # DBD doesn't use cursors; loads whole table during execute()
  # problem for large tables - "mysql_use_result" gets around this
  my $sth;
  if (defined $self->{targetdb}{tables}{$targettable}{its_a_big_un}) {
    $sth = $self->source_dbh()->prepare( $select , { mysql_use_result => 1 } );
  }
  else {
    $sth = $self->source_dbh()->prepare( $select );
  }
  $sth->execute();
  
  my $row;
  if( exists $self->{targetdb}{tables}{$targettable}{row_modify} ) {
    my $rowmod = $self->{targetdb}{tables}{$targettable}{row_modify};

    while( my $arref = $sth->fetchrow_arrayref() ) {

      $row = &{$rowmod}($arref);

      print $tmpfile ( join( "\t",@{$row} ),"\n" );
    }
  } else {
    while( my $arref = $sth->fetchrow_arrayref() ) {
      $row = join( "\t", @$arref );
      print $tmpfile "$row\n";
    }
  }
}



sub read_dbs {
  my $self = shift;

  my $dbh;

  for my $db_name ('targetdb', 'sourcedb' ) {
    if( $db_name eq 'targetdb' ) {
      $dbh = $self->target_dbh();
    } else {
      $dbh = $self->source_dbh();
    }

    my $sth = $dbh->prepare( "show tables" );
    $sth->execute();

    while( my $arref = $sth->fetchrow_arrayref() ) {
      $self->{$db_name}{tables}{$arref->[0]} = {};
    }

    my @tables = keys %{$self->{$db_name}{tables}};
    for my $table ( @tables ) {
      $sth = $dbh->prepare( "show columns from $table" );
      $sth->execute();
      while( my $arref = $sth->fetchrow_arrayref () ) {
	push( @{$self->{$db_name}{tables}{$table}{columns}}, $arref->[0] );
      }
    }
  }
}


sub table_rename {
  my ( $self, $oldtable, $newtable ) = @_;
  $self->{targetdb}{tables}{$newtable}{link} = $oldtable;
}

sub table_skip {
  my ( $self, $newtable ) = @_;
  $self->{targetdb}{tables}{$newtable}{link} = "";
}

sub column_rename {
  my ( $self, $newtable, $oldcol, $newcol ) = @_;
  $self->{targetdb}{tables}{$newtable}{columnrename}{$newcol} = $oldcol;
}

sub column_skip {
  my ( $self, $newtable,  $newcol ) = @_;
  $self->{targetdb}{tables}{$newtable}{columnrename}{$newcol} = "";
}

sub custom_select {
  my ( $self, $newtable, $select ) = @_;
  $self->{targetdb}{tables}{$newtable}{select} = $select;
}

sub set_row_modifier {
  my ( $self, $newtable, $row_modifier ) = @_;
  $self->{targetdb}{tables}{$newtable}{row_modify} = $row_modifier;
}
  
sub clear_target {
  my $self = shift;
  for my $tablename ( keys %{$self->{targetdb}{tables}} ) {
    $self->target_dbh()->do( "delete from $tablename" );
  }
}

sub big_table {
  my ( $self, $table ) = @_;

  $self->{targetdb}{tables}{$table}{its_a_big_un} = 1;
}


sub do_first {
  my $self = shift;
  my @ordered_table_list= @_;
  $self->{targetdb}{dofirst} = \@ordered_table_list;
}

1;
