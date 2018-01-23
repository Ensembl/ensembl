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

Bio::EnsEMBL::Utils::SqlHelper

=head1 VERSION

$Revision$

=head1 SYNOPSIS

  use Bio::EnsEMBL::Utils::SqlHelper;

  my $helper =
    Bio::EnsEMBL::Utils::SqlHelper->new( -DB_CONNECTION => $dbc );

  my $arr_ref = $helper->execute(
    -SQL      => 'select name, age from tab where col =?',
    -CALLBACK => sub {
      my @row = @{ shift @_ };
      return { name => $row[0], age => $row[1] };
    },
    -PARAMS => ['A'] );

  use Data::Dumper;
  print Dumper($arr_ref), "\n";
  # Prints out [name=>'name', age=>1] maybe ....


  # For transactional work; only works if your MySQL table
  # engine/database supports transactional work (such as InnoDB)

  $helper->transaction(
    -CALLBACK => sub {
      if ( $helper->execute_single_result(
                                      -SQL => 'select count(*) from tab'
           ) )
      {
        return $helper->execute_update('delete from tab');
      } else {
        return
          $helper->batch( -SQL  => 'insert into tab (?,?)',
                          -DATA => [ [ 1, 2 ], [ 1, 3 ], [ 1, 4 ] ] );
      }
    } );

=head1 DESCRIPTION

Easier database interaction

=cut

package Bio::EnsEMBL::Utils::SqlHelper;

use warnings;
use strict;

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref check_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Utils::Iterator;
use English qw( -no_match_vars ); #Used for $PROCESS_ID
use Scalar::Util qw(weaken); #Used to not hold a strong ref to DBConnection
use Time::HiRes;

=pod

=head2 new()

  Arg [DB_CONNECTION] : Bio::EnsEMBL::DBSQL::DBConnection $db_connection
  Returntype          : Instance of helper
  Exceptions          : If the object given as a DBConnection is not one or it
                        was undefined
  Status              : Stable
  Description         : Creates a new instance of this object.
  Example             :
  my $dba = get_dba('mydb');    # New DBAdaptor from somewhere
  my $helper = Bio::EnsEMBL::Utils::SqlHelper->new(
    -DB_CONNECTION => $dba->dbc() 
  );

  $helper->execute_update( -SQL    => 'update tab set flag=?',
                           -PARAMS => [1] );

=cut

sub new {
	my ( $class, @args ) = @_;
	
	my ($db_connection) = rearrange([qw(db_connection)], @args);
	 
	my $self = bless( {}, ref($class) || $class );
	throw('-DB_CONNECTION construction parameter was undefined.') 
	 unless defined $db_connection;
	$self->db_connection($db_connection);
	
	return $self;
}

=pod

=head2 db_connection()

  Arg [1]     : Bio::EnsEMBL::DBSQL::DBConnection $db_connection
  Description : Sets and retrieves the DBConnection 
  Returntype  : Bio::EnsEMBL::DBSQL::DBConnection
  Exceptions  : If the object given as a DBConnection is not one or if an 
                attempt is made to set the value more than once
  Status      : Stable

=cut

sub db_connection {
  my ($self, $db_connection) = @_;
  if(defined $db_connection) {
    if(exists $self->{db_connection}) {
      throw('Cannot reset the DBConnection object; already defined ');
    }
    assert_ref($db_connection, 'Bio::EnsEMBL::DBSQL::DBConnection', 'db_connection');
    $self->{db_connection} = $db_connection;
    weaken $self->{db_connection};
  }
  return $self->{db_connection};
}

# --------- SQL Methods

=pod

=head2 execute() - Execute a SQL statement with a custom row handler

  Arg [SQL]             : string SQL to execute
  Arg [CALLBACK]        : CodeRef; The callback to use for mapping a row to a data  
                          point; leave blank for a default mapping to a 2D array
  Arg [USE_HASHREFS]    : boolean If set to true will cause HashRefs to be returned 
                          to the callback & not ArrayRefs
  Arg [PARAMS]          : ArrayRef The binding parameters to the SQL statement
  Arg [PREPARE_PARAMS]  : boolean Parameters to be passed onto the Statement Handle 
                          prepare call
  Arg [ITERATOR]        : boolean Request a Bio::EnsEMBL::Utils::Iterator
                          rather than a 2D array
  Returntype :  ArrayRef or Bio::EnsEMBL::Utils::Iterator
  Exceptions :  If errors occur in the execution of the SQL
  Status     :  Stable
  Example    :
  my $arr_ref = $helper->execute(
    -SQL      => 'select a,b,c from tab where col =?',
    -CALLBACK => sub {
      my @row = @{ shift @_ };
      return { A => $row[0], B => $row[1], C => $row[2] };
    },
    -PARAMS => ['A'] );

  #Or with hashrefs
  my $arr_ref = $helper->execute(
    -SQL          => 'select a,b,c from tab where col =?',
    -USE_HASHREFS => 1,
    -CALLBACK     => sub {
      my $row = shift @_;
      return { A => $row->{a}, B => $row->{b}, C => $row->{c} };
    },
    -PARAMS => ['A'] );

  Description:
Uses a callback defined by the sub decalaration. Here we specify how
the calling code will deal with each row of a database's result set. The
sub can return any type of Object/hash/data structure you require.

Should you not specify a callback then a basic one will be assigned to
you which will return a 2D array structure e.g.

  my $arr_ref = $helper->execute(
                           -SQL => 'select a,b,c from tab where col =?',
                           -PARAMS => ['A'] );

This is equivalent to DBI's selectall_arrayref() subroutine.

As an extension to this method you can write a closure subroutine which
takes in two parameters. The first is the array/hash reference & the
second is the statement handle used to execute. 99% of the time you will
not need it but there are occasions where you do need it. An example of
usage would be:

  my $conn = get_conn();    #From somwewhere
  my $arr_ref = $conn->execute(
    -SQL          => 'select a,b,c from tab where col =?',
    -USE_HASHREFS => 1,
    -CALLBACK     => sub {
      my ( $row, $sth ) = @_;
      #Then do something with sth
      return { A => $row->[0], B => $row->[1], C => $row->[2] };
    },
    -PARAMS => ['A'] );

Any arguments to bind to the incoming statement. This can be a set of scalars
or a 2D array if you need to specify any kind of types of sql objects i.e.

  use DBI qw(:sql_types);

  my $conn = get_conn();
  my $arr_ref = $conn->execute(
    -SQL =>
      'select a,b,c from tab where col =? and num_col=? and other=?',
    -USE_HASHREFS => 1,
    -CALLBACK     => sub {
      my @row = @{ shift @_ };
      return { A => $row[0], B => $row[1], C => $row[2] };
    },
    -PARAMS => [ '1', SQL_VARCHAR ],
    [ 2, SQL_INTEGER ],
    'hello' );

Here we import DBI's sql types into our package and then pass in
multiple anonymous array references as parameters. Each param is
tested in the input and if it is detected to be an ARRAY reference we
dereference the array and run DBI's bind_param method. In fact you can
see each part of the incoming paramaters array as the contents to call
C<bind_param> with. The only difference is the package tracks the bind
position for you.

We can get back a L<Bio::EnsEMBL::Utils::Iterator> object which can be used
to iterate over the results set without first materializing the data into 
memory. An example would be:

   my $iterator = $helper->execute(
                           -SQL => 'select a,b,c from tab where col =?',
                           -PARAMS => ['A'] 
                           -ITERATOR => 1);
   while($iterator->has_next()) {
     my $row = $iterator->next();
     #Do something
   }

This is very useful for very large datasets.

=cut

sub execute {
	my ( $self, @args ) = @_;
	my ($sql, $callback, $use_hashrefs, $params, $prepare_params, $iterator) = 
	 rearrange([qw(sql callback use_hashrefs params prepare_params iterator)], @args);
	my $has_return = 1;
	
	#If no callback then we execute using a default one which returns a 2D array
	if(!defined $callback) {
    if($use_hashrefs) {
      $callback = $self->_mappers()->{hash_ref};
    }
    else {
      $callback = $self->_mappers()->{array_ref};
    }    
	}
	
	return $self->_execute( $sql, $callback, $has_return, $use_hashrefs, $params, $prepare_params, $iterator );
}

=pod

=head2 execute_simple()

  Arg [SQL]           : string $sql
  Arg [PARAMS]        : ArrayRef $params
  Arg [CALLBACK]      : CodeRef $callback
  Returntype : ArrayRef of 1D elements
  Exceptions : If errors occur in the execution of the SQL
  Status     : Stable
  Example : my $classification =
    $helper->execute_simple(
       -SQL => 'select meta_val from meta where meta_key =? order by meta_id',
       -PARAMS => ['species.classification'] 
    );
  Description: Similar to execute() but without a sub-routine reference. 
               Using this code assumes you want an array of single scalar values 
               as returned by the given SQL statement.
=cut

sub execute_simple {
  my ( $self, @args ) = @_;
	my ($sql, $params, $callback) = rearrange([qw(sql params callback)], @args);
	my $has_return = 1;
	my $use_hashrefs = 0;
	$callback ||= $self->_mappers()->{first_element};
	return $self->_execute($sql, $callback, $has_return, $use_hashrefs, $params);
}

=pod

=head2 execute_no_return()

  Arg [SQL]           : string sql
  Arg [CALLBACK]      : CodeRef The callback to use for mapping a row to a data point;
                        we assume you are assigning into a data structure which
                        has requirements other than simple translation into an
                        array
  Arg [USE_HASHREFS]  : boolean If set to true will cause HashRefs to be returned 
                        to the callback and not ArrayRefs
  Arg [PARAMS]        : ArrayRef The binding parameters to the SQL statement
  Returntype : None
  Exceptions : If errors occur in the execution of the SQL
  Status     : Stable
  Description: Whilst all other execute methods will return something; this assumes that the
               given mapper subroutine will be performing the business of placing values
               somewhere or doing something with them.

               There is a huge temptation to nest queries using this method; do not! Execute
               the values into an array using one of the other methods then run your subqueries
               on them; or make a better first query. SQL is flexible; so use it.

=cut

sub execute_no_return {
	my ( $self, @args ) = @_;
	my ($sql, $callback, $use_hashrefs, $params) = rearrange([qw(sql callback use_hashrefs params)], @args);
	throw('No callback defined but this is a required parameter for execute_no_return()') if ! $callback;
	my $has_return = 0;
	my $prepare_params = [];
	$self->_execute( $sql, $callback, $has_return, $use_hashrefs, $params);
	return;
}

=pod

=head2 execute_into_hash()

  Arg [SQL]           : string $sql
  Arg [CALLBACK]      : CodeRef The callback to use for mapping to a value in a hash
                        keyed by the first element in your result set; 
                        leave blank for a default mapping to a scalar value
                        of the second element
  Arg [PARAMS]        : The binding parameters to the SQL statement
  Returntype : HashRef keyed by column 1 & value is the return of callback
  Exceptions : If errors occur in the execution of the SQL
  Status     : Stable
  Description: A variant of the execute methods but rather than returning a list of
               mapped results, this will assume the first column of a returning map and
               the calling subroutine will map the remainder of your return as the
               hash's key.

               This code can handle simple queries to hashes, complex value mappings
               and repeated mappings for the same key.



  Example:

  my $sql    = 'select key, one, two from table where something =?';
  my $mapper = sub {
    # Argument one is a row from the result, argument two is any previously seen value for the first column of the row
    my ( $row, $value ) = @_;
    #Ignore field 0 as that is being used for the key
    my $obj = Some::Obj->new( one => $row->[1], two => $row->[2] );
    return $obj;
  };

  my $hash = $helper->execute_into_hash( -SQL      => $sql,
                                         -CALLBACK => $mapper,
                                         -PARAMS   => ['val'] );

  #Or the default simplistic invocation
  my $sql = 'select biotype, count(gene_id) from gene group by biotype';
  my $biotype_hash = $conn->execute_into_hash( -SQL => $sql );
  print $biotype_hash->{protein_coding} || 0, "\n";

  # More complicated mapping, result hash will be keyed on "meta_key"
  # Hash will contain lists of values linked to their meta_key
  my %args = ( -SQL => 'select meta_key, meta_value from meta '
                 . 'where meta_key =? order by meta_id',
               -PARAMS => ['species.classification'] );

  my $hash = $helper->execute_into_hash(
    %args,
    -CALLBACK => sub {
      my ( $row, $value ) = @_;
      $value = [] if !defined $value;
      push( @{$value}, $row->[1] );
      return $value;
    } );

  # Add new values to an already seen existing meta_key
  $hash = $helper->execute_into_hash(
    %args,
    -CALLBACK => sub {
      my ( $row, $value ) = @_;
      if ( defined $value ) {
        # Calling code is dealing with $row->[0], so we only have to handle the remaining columns
        push( @{$value}, $row->[1] );
        return;
      }
      my $new_value = [ $row->[1] ];
      return $new_value;
    } );

=cut

sub execute_into_hash {
	my ( $self, @args ) = @_;
	my ($sql, $callback, $params) = rearrange([qw(sql callback params)], @args);
	my $hash = {};
	
	#If no callback then we execute using a default one which sets value to 2nd element
	if(!defined $callback) {
	 $callback = $self->_mappers()->{second_element};
	}
	
	#Default mapper uses the 1st key + something else from the mapper
	my $mapper = sub {
		my $row = shift @_;
		my $key = $row->[0];
		my $value = $hash->{$key};
		my $new_value = $callback->($row, $value);
		if(defined $new_value) {
		 $hash->{ $key } = $new_value;
		}
		return;
	};
	
	$self->execute_no_return(
	  -SQL => $sql, 
	  -CALLBACK => $mapper,
	  -PARAMS => $params
	);
	
	return $hash;
}

=pod

=head2 execute_single_result()

  Arg [SQL]           : string $sql
  Arg [CALLBACK]      : CodeRef The callback to use for mapping a row to a data point; 
                        leave blank for a default scalar mapping
  Arg [USE_HASHREFS]  : boolean If set to true will cause HashRefs to be returned 
                        to the callback & not ArrayRefs
  Arg [PARAMS]        : ArrayRef The binding parameters to the SQL statement
  Arg [NO_ERROR]      : Boolean Flag to indicate that the code should not 
                        throw an error when row counts are not equal to 1
  Returntype : Scalar
  Exceptions : If errors occur in the execution of the SQL, if the query 
               returned more than 1 row and if we found no rows.
  Status     : Stable
  Example    : 
  my $meta_count =
    $helper->execute_single_result(
                -SQL => 'select count(*) from meta where species_id =?',
                -PARAMS => [1] );

  Description : Very similar to execute() except it will raise an exception if we have more 
                or less than one row returned
=cut

sub execute_single_result {
	my ( $self, @args ) = @_;
	my ($sql, $callback, $use_hashrefs, $params, $no_error) = rearrange(
	 [qw(sql callback use_hashrefs params no_error)], @args);
	
	my $results = $self->execute_simple( 
	  -SQL => $sql, 
	  -CALLBACK => $callback, 
	  -USE_HASHREFS => $use_hashrefs, 
	  -PARAMS => $params
	);
	
	my $result_count = scalar(@{$results});
	if(! $no_error && $result_count != 1) {
	  $params = [] if ! $params;
	  my $type = ($result_count == 0) ? 'No' : 'Too many';
		my $msg = "${type} results returned. Expected 1 but got $result_count for query '${sql}' with params [";
		$msg .= join( ',', map {(defined $_) ? $_ : '-undef-';} @{$params} );
		$msg .= ']';
		throw($msg);
	}
	return $results->[0] if defined $results->[0];
	return;
}

=pod

=head2 transaction()

  Arg [CALLBACK]      : CodeRef The callback used for transaction isolation; once 
                        the subroutine exists the code will decide on rollback
                        or commit. Required
  Arg [RETRY]         : integer the number of retries to attempt with this 
                        transactional block. Defaults to 0. 
  Arg [PAUSE]         : integer the time in seconds to pause in-between retries.
                        Defaults to 1. Fractions are allowed as use delegate to
                        Time::HiRes' sleep function
  Arg [CONDITION]     : CodeRef allows you to inspect the exception raised
                        and should your callback return true then the 
                        retry will be attempted. If not given then all 
                        exceptions mean attempt a retry (if specified)
  Returntype : Return of the callback
  Exceptions : If errors occur in the execution of the SQL
  Status     : Stable
  Example :
  my $val = $helper->transaction(
    -CALLBACK => sub {
      my ($dbc) = @_;
      #Do something
      return 1;
    } );

  # Or without named arguments
  my $val = $helper->transaction(
    sub {
      my ($dbc) = @_;
      #Do something
      return 1;
    } );

  # To allow retries (use sparingly, and see description)
  my $val = $helper->transaction(
    -RETRY => 3, 
    -PAUSE => 2,
    -CALLBACK => sub {
      my ($dbc) = @_;
      #Do something
      return 1;
    } );
 
  # Only retry when we have an error containing the phrase "deadlock"
  my $val = $helper->transaction(
    -RETRY => 3, -PAUSE => 2,
    -CALLBACK => sub {
      my ($dbc) = @_;
      #Do something
      return 1;
    },
    -CONDITION => sub {
      my ($error) = @_;
      return ( $error =~ /deadlock/ ) ? 1 : 0;
    }
  );

  Description :
Creates a transactional block which will ensure that the connection is
committed when your submmited subroutine has finished or will rollback
in the event of an error occuring in your block.

The code will always force AutoCommit off but will restore it to its
previous setting. If your DBI/DBD driver does not support manual
commits then this code will break. The code will turn off the
disconnect_when_idle() method to allow transactions to work as
expected.

An effect of using REPEATABLE READ transaction isolation (InnoDB's
default) is that your data is as fresh as when you started your current
transaction. To ensure the freshest data use SELECT ... FROM ... LOCK
IN SHARE MODE> or SELECT ... FROM ... LOCK FOR UPDATE if you are
going to issue updates.

Creating a transaction within a transaction results in the commit
rollback statements occuring in the top level transaction. That way any
block of code which is meant to to be transaction can be wrapped in
this block ( assuming the same instance of SQLHelper is passed around and
used).

You can also request the -RETRY of a transactional block of code which is
causing problems. This can indicate your programming model is broken, 
so use with care.
    
The -RETRY argument indicates the number of times we attempt the transaction 
and -PAUSE indicates the time in-between attempts. These retries will
only occur in the root transaction block i.e. you cannot influence the 
retry system in a sub-transaction. You can influence if the retry is done with
the -CONDITION argument which accepts a Code reference (same as the
-CALLBACK parameter). This allows you to inspect the error thrown to
retry only in some situations e.g.

=cut

sub transaction {
  my ($self, @args) = @_;
   
  my ($callback, $retry, $pause, $condition) = rearrange([qw(callback retry pause condition)], @args);
  
  throw('-CALLBACK was not a CodeRef. Got a reference of type ['.ref($callback).']. Check your parameters') 
    unless check_ref($callback, 'CODE');
  
  #Setup defaults
  $retry = 0 unless defined $retry;
  $pause = 1 unless defined $pause;
  if(! defined $condition) {
    $condition = sub {
      return 1;
    };
  }
  
  assert_ref($condition, 'CODE', '-CONDITION');
 
  my $dbc = $self->db_connection();
  my $original_dwi;
  my $ac;
  
  my $error;
  my $result;
  
  #If we were already in a transaction then we do not do any management of the
  #session & wait for the parent transaction(s) to finish
  my $perform_transaction = $self->_perform_transaction_code();
  if($perform_transaction) {
    ($original_dwi, $ac) = $self->_enable_transaction();
  }
  else {
    #If we were in a transaction then ignore any attempts at retry here
    $retry = 0;
  }
    
  for(my $iteration = 0; $iteration <= $retry; $iteration++) {
    eval {
      $result = $callback->($dbc);
      $dbc->db_handle()->commit() if $perform_transaction;
    };
    $error = $@;
    #If we were allowed to deal with the error then we apply rollbacks & then
    #retry or leave to the remainder of the code to throw
    if($perform_transaction && $error) {
      eval { $dbc->db_handle()->rollback(); };
      #If we were not on our last iteration then warn & allow the retry
      if($iteration != $retry) {
        if($condition->($error)) {
          warn("Encountered error on attempt ${iteration} of ${retry} and have issued a rollback. Will retry after sleeping for $pause second(s): $error");
          Time::HiRes::sleep $pause;
        }
        else {
          last; #break early if condition of error was not matched
        }
      }
    }
    
    #Always break the loop if we had a successful attempt
    last if ! $error;
  }
  
  if($perform_transaction) {
    $self->_disable_transaction($original_dwi, $ac);
  }
  
  throw("ABORT: Transaction aborted because of error: ${error}") if $error;
  
  return $result;
}

=pod

=head2 execute_update()

  Arg [SQL]           : string $sql
  Arg [CALLBACK]      : CodeRef The callback to use for calling methods on the 
                        DBI statement handle or DBConnection object after an 
                        update command
  Arg [PARAMS]        : ArrayRef The binding parameters to the SQL statement
  Arg [PREPARE_PARAMS] : ArrayRef Parameters to bind to the prepare() StatementHandle call
  Returntype : integer - Number of rows affected
  Exceptions : If errors occur in the execution of the SQL
  Status     : Stable
  Description: Used for performing updates but conforms to the normal execute statement subroutines.
  Example    : 
  use DBI qw(:sql_types);
  $helper->execute_update(-SQL => 'update tab set name = ? where id =?',
                          -PARAMS => [ 'andy', [ 1, SQL_INTEGER ] ] );

#If you need to do something a bit more advanced with your database then you can
#give the method a closure and this will be called after the execute has been
#issued i.e.

  my $obj;
  $helper->execute_update(
    -SQL      => 'insert into tab (name) values(?)',
    -CALLBACK => sub {
      my ( $sth, $dbh, $rv ) = @_;
      $obj->{id} = $dbh->{mysql_insertid};
    },
    -PARAMS => [ $obj->name() ] );

#This lets us access the statement handle, the database handle and
#the return value from $sth->execute, to access other properties such as
#the last identifier inserted.

=cut

sub execute_update {
  my ($self, @args) = @_;
  my ($sql, $callback, $params, $prepare_params) = rearrange([qw(sql callback params prepare_params)], @args);
  my $rv = 0;
  my $sth;
  eval {
    my @prepare_params;
    @prepare_params = @{$prepare_params} if check_ref($prepare_params, 'ARRAY');
    $sth = $self->db_connection()->prepare($sql, @prepare_params);
    $self->_bind_params($sth, $params);
    $rv = $sth->execute();
    $callback->($sth, $self->db_connection()->db_handle(), $rv) if $callback;
  };
  my $error = $@;
  $self->_finish_sth($sth);
  if($error) {
    my $params = join ' ', map { (defined $_) ? $_ : q{undef} } @{$params};
    throw("Cannot apply sql '${sql}' with params '${params}': ${error}");
  }
  return $rv;
}

=head2 execute_with_sth()

  Arg [SQL]             : string $sql
  Arg [CALLBACK]        : CodeRef The callback to use for working with the statement
                          handle once returned. This is not a mapper.
  Arg [PARAMS]          : ArrayRef The binding parameters to the SQL statement
  Arg [PREPARE_PARAMS]  : ArrayRef Used to pass parameters to the statement handle 
                          prepare method
  Description : Run query without worrying statement handles and such.
Very similar to execute() except this gives you full control over the
lifecycle of the statement handle and how you wish to proceed with working
with a statement handle. This is for situations where you believe going through
the mappers causes too much of a slow-down (since we have to execute a
subroutine for every row in order to map it correctly).

However please benchmark before adopting this method as it increases the 
complexity of your code and the mapper slowness only becomes apparent when
working with very large numbers of rows.
  Returntype  : Anything you wish to return from the callback
  Exceptions  : If errors occur in the execution of the SQL
  Status      : Stable
  Example     : 
  my $meta_count = $helper->execute_with_sth(
    -SQL      => 'select count(*) from meta where species_id =?',
    -PARAMS   => [1],
    -CALLBACK => sub {
      my ($sth) = @_;
      my $count;
      $sth->bind_columns( \$count );
      while ( $sth->fetch ) {
        print $count, "\n";
      }
      return $count;
    } );



=cut

sub execute_with_sth {
  my ($self, @args) = @_;
  my ($sql, $callback, $params, $prepare_params) = rearrange([qw(sql callback params prepare_params)], @args);
  my $sth = $self->_base_execute( $sql, $params, $prepare_params, $callback );
  my $result = eval {$callback->($sth)};
  my $error = $@;
  $self->_finish_sth($sth);
  die $error if $error;
  return $result;
}

=pod

=head2 batch()

  Arg [SQL]           : string $sql
  Arg [CALLBACK]      : CodeRef The callback to use for working with the statement
                        handle once returned; specify this or -DATA
  Arg [DATA]          : ArrayRef The data to insert; specify this or -CALLBACK
  Arg [COMMIT_EVERY]  : Integer defines the rate at which to issue commits to
                        the DB handle. This is important when working with 
                        InnoDB databases since it affects the speed of rollback
                        (larger gaps inbetween commits means more to rollback).
                        
                        Ignored if using the callback version.
  Arg [PREPARE_PARAMS]  : ArrayRef Used to pass parameters to the statement handle 
                          prepare method
  Returntype : integer rows updated
  Exceptions : If errors occur in the execution of the SQL
  Status     : Stable
  Example    :
  my $alotofdata = getitfromsomewhere();
  $helper->batch(
    -SQL      => 'insert into table (one,two) values(?,?)',
    -CALLBACk => sub {
      my ( $sth, $dbc ) = @_;
      foreach my $data (@alotofdata) {
        $sth->execute( @{$data} );
      }
    } );

  #Or for a 2D array data driven approach
  $helper->batch( -SQL  => 'insert into table (one,two) values(?,?)',
                  -DATA => $alotofdata );


  Description: Takes in a sql statement and a code reference. Your SQL is converted into a 
prepared statement and then given as the first parameter to the closure. The
second parameter is the DBH which created the statement. This is intended
to let you do mass insertion into a database without the need to
re-preparing the same statement.

This can be combined with the transaction() code to provide a construct
which does batch insertion and is transactionally aware.

We can also use data based batch insertions i.e.

  #Needs to be like:
  #   [ [1,2], [3,4] ]
  #Or if using the DBI types:
  #  [ [ [ 1, SQL_INTEGER ], [ 2, SQL_INTEGER ] ],
  #    [ [ 3, SQL_INTEGER ], [ 4, SQL_INTEGER ] ] ];

  my $alotofdata = getitfromsomewhere();
  $helper->batch( -SQL  => 'insert into table (one,two) values(?,?)',
                  -DATA => $alotofdata );

This does exactly the same as the previous example.

All batch statements will return the value the callback computes. If you are 
using the previous example with a data array then the code will return the
number affected rows by the query.

=cut

sub batch {
  my ($self, @args) = @_;
  my ($sql, $callback, $data, $commit_every, $prepare_params) = 
    rearrange([qw(sql callback data commit_every prepare_params)], @args);
  
  if(! defined $callback && ! defined $data) {
    throw('You need to define a callback for insertion work or the 2D data array');
  }
  
  my $result;
  if(defined $callback) {
    $result = $self->_callback_batch($sql, $callback, $prepare_params);
  }
  else {
    $result = $self->_data_batch($sql, $data, $commit_every, $prepare_params);
  }
  return $result if defined $result;
  return;
}

#------- Internal methods

my $default_mappers = {
  first_element => sub {
    my ($row) = @_;
    return $row->[0];
  },
  second_element => sub {
    my ($row) = @_;
    return $row->[1];
  },
  array_ref => sub {
    my $row = shift @_;
    return [@{$row}]; #copy of array done because DBI's array is read only
  },
  hash_ref => sub {
    my $row = shift @_;
    return {%{$row}}; #applying same logic as above
  }
}; 

sub _mappers {
  my ($self) = @_;
  return $default_mappers;
}

sub _perform_transaction_code {
  my ($self) = @_;
  return $self->{_transaction_active}->{$PROCESS_ID} ? 0 : 1;
}

sub _enable_transaction {
  my ($self) = @_;
  my $dbc = $self->db_connection();
  my $original_dwi = $dbc->disconnect_when_inactive();
  $dbc->disconnect_when_inactive(0);
  my $ac = $dbc->db_handle()->{'AutoCommit'};
  $dbc->db_handle()->{'AutoCommit'} = 0;
  $self->{_transaction_active}->{$PROCESS_ID} = 1;
  return ($original_dwi, $ac);
}

sub _disable_transaction {
  my ($self, $original_dwi, $ac) = @_;
  my $dbc = $self->db_connection();
  $dbc->db_handle()->{'AutoCommit'} = $ac;
  $dbc->disconnect_when_inactive($original_dwi);
  delete $self->{_transaction_active}->{$PROCESS_ID};
  return;
}

sub _bind_params {
	my ( $self, $sth, $params ) = @_;
	
	return if ! defined $params; #Return quickly if we had no data
	
	if(! check_ref($params, 'ARRAY')) {
	  throw(qq{The given parameters reference '${params}' is not an ARRAY; wrap in an ArrayRef});
	}
	
	my $count = 1;
	foreach my $param (@{$params}) {
		if ( check_ref($param, 'ARRAY') ) {
			$sth->bind_param( $count, @{$param} );
		}
		else {
			$sth->bind_param( $count, $param );
		}
		$count++;
	}
	return;
}

sub _execute {
	my ( $self, $sql, $callback, $has_return, $use_hashrefs, $params, $prepare_params, $iterator ) = @_;

	throw('Not given a mapper. _execute() must always been given a CodeRef') unless check_ref($callback, 'CODE');
	
  my $sth = $self->_base_execute($sql, $params, $prepare_params);
  
  my $sth_processor;
  if($use_hashrefs) {
    $sth_processor = sub {
      while( my $row = $sth->fetchrow_hashref() ) {
        my $v = $callback->($row, $sth);
        return $v if $has_return;
      }
      $self->_finish_sth($sth);
      return;
    };
  }
  else {
    $sth_processor = sub {
      while( my $row = $sth->fetchrow_arrayref() ) {
        my $v = $callback->($row, $sth);
        return $v if $has_return;
      }
      $self->_finish_sth($sth);
      return;
    };
  }
	
  my $iter = Bio::EnsEMBL::Utils::Iterator->new($sth_processor);
  if($has_return) {
    return $iter if $iterator;
    return $iter->to_arrayref();
  }
  else {
    #Force iteration if we had no return since the caller is expecting this
    $iter->each(sub {});
  }
  return;
}

sub _base_execute {
  my ( $self, $sql, $params, $prepare_params) = @_;
	
	$params = [] unless $params;
	
	my $conn = $self->db_connection;
	
	my $sth;
	eval {
	  my @prepare_params;
	  @prepare_params = @{$prepare_params} if check_ref($prepare_params, 'ARRAY');
		$sth = $conn->prepare($sql, @prepare_params);
		throw("Cannot continue as prepare() did not return a handle with prepare params '@prepare_params'") 
		  unless $sth;
		$self->_bind_params( $sth, $params );
		$sth->execute();
	};
	
	my $error = $@;
	if($error) {
  	throw("Cannot run '${sql}' with params '@{$params}' due to error: $error") if $error;
	}
	
	return $sth;
}

sub _finish_sth {
  my ($self, $sth) = @_;
  eval { $sth->finish() if defined $sth; };
  warn('Cannot finish() the statement handle: $@') if $@;
  return;
}

sub _callback_batch {
  my ($self, $sql, $callback, $prepare_params) = @_;
  my $error;
  my $sth;
  my $closure_return;
  eval {
    my @prepare_params;
    @prepare_params = @{$prepare_params} if check_ref($prepare_params, 'ARRAY');
    $sth = $self->db_connection()->prepare($sql, @prepare_params); 
    $closure_return = $callback->($sth, $self->db_connection());
  };
  $error = $@;
  $self->_finish_sth($sth);
	throw("Problem detected during batch work: $error") if $error;
  
  return $closure_return if defined $closure_return;
  return;
}

sub _data_batch {
  my ($self, $sql, $data, $commit_every, $prepare_params) = @_;
  
  #Input checks
  assert_ref($data, 'ARRAY', '-DATA');
  my $data_length = scalar(@{$data});
  return 0 unless $data_length > 0;
  my $first_row = $data->[0];
  throw('I expect to work with a 2D ArrayRef but this is not one') unless check_ref($first_row, 'ARRAY');

  my $callback = sub {
    my ($sth, $dbc) = @_;
    my $total_affected = 0;
    #Iterate over each data point
    for(my $data_index = 0; $data_index < $data_length; $data_index++) {
      my $row = $data->[$data_index];
      $self->_bind_params($sth, $row);
      my $affected = eval {$sth->execute()};
      if($@) {
        throw("Problem working with $sql with params @{$row}: $@");
      }
      my $num_affected = ($affected) ? $affected :  0; #Get around DBI's 0E0
      $total_affected += $num_affected;
      
      #Lets us do a commit once every x rows apart from 0. We also finish
      #off with a commit if the code told us we were doing it
      if($commit_every) {
        if( ($data_index % $commit_every == 0) && $data_index != 0) {
          $dbc->db_handle()->commit();
        }
      }
    }
    
    #finish off with a commit if the code told us we were doing it
    if($commit_every) {
      $dbc->db_handle()->commit();
    }
    
    return $total_affected || 0;
  };
  
  return $self->_callback_batch($sql, $callback, $prepare_params)
}

1;
