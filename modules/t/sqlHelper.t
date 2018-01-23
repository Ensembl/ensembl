# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#A set tests used to prod the SqlHelper class

use strict;
use warnings;

use Test::More;
use Test::Warnings;
use Test::Exception;
use Scalar::Util qw(isweak);

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Utils::SqlHelper;

#Redefine the WARN sig to note the errors (most are just from transaction retry)
$SIG{__WARN__} = sub {
  note @_;
  return 1; 
};

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $dba = $multi->get_DBAdaptor( 'core' );
ok( $dba, 'Test database instatiated' );

#Now start testing the Helper
dies_ok { Bio::EnsEMBL::Utils::SqlHelper->new() } 'Expect to die when no DBConnection was given'; 
dies_ok { Bio::EnsEMBL::Utils::SqlHelper->new(-DB_CONNECTION => $dba) } 
  'Expect to die when we do not give SqlHelper a DBConncetion'; #was given a DBAdaptor
ok ( 
  isweak(Bio::EnsEMBL::Utils::SqlHelper->new(-DB_CONNECTION => $dba->dbc())->{db_connection}),
  'Checking DBConnection reference is weak when we ask for it' 
);

my $helper = Bio::EnsEMBL::Utils::SqlHelper->new(-DB_CONNECTION => $dba->dbc());
ok ( $helper, 'SqlHelper instance was created' );

my $meta_key = 'species.common_name';
note("Meta key queries working with ${meta_key}. If the tests fail then check for it in the DB dumps");

is( 
  $helper->execute_single_result(-SQL => qq{select count(*) from meta where meta_key = '$meta_key'}),
  1,
  'Checking count of meta key is right with no params'
);

is( 
  $helper->execute_single_result(-SQL => 'select count(*) from meta where meta_key =?', -PARAMS => [$meta_key]),
  1,
  'Checking count of meta key is right with params'
);


throws_ok { $helper->execute_single_result(-SQL => 'select * from meta') } qr/Too many results/, 'More than 1 row causes an error';
throws_ok { $helper->execute_single_result(-SQL => 'select * from meta where species_id =?', -PARAMS => [-1]) } qr/No results/, 'Less than 1 row causes an error';

is(
  $helper->execute_single_result(-SQL => 'select meta_id from meta order by meta_id', -NO_ERROR => 1),
  1,
  'Checking if we have more than row we will not error if we ask to ignore it'
);

is(
  $helper->execute_single_result(-SQL => 'select meta_id from meta where species_id =?', -PARAMS => [-1], -NO_ERROR => 1),
  undef,
  'Checking if we less than row we will not error if we ask to ignore it'
);

is_deeply( 
  $helper->execute(-SQL => 'select count(*), 3 from meta where meta_key =?', -PARAMS => [$meta_key])->[0],
  [1,3],
  'Checking 2D mapping of meta key count works'
);

#EXECUTE_INTO_HASH() CHECKS

my $meta_count_hash = $helper->execute_into_hash(
  -SQL => 'select meta_key, count(*) from meta group by meta_key'
);

is($meta_count_hash->{$meta_key}, 1, 'Checking hash comes back correctly');

{
  my $count = 0;
  my %args = (
    -SQL => 'select meta_key, meta_value from meta where meta_key =? order by meta_id',
    -PARAMS => ['species.classification']
  );
  my $expected_hash = {
    'species.classification' => [
      'Homo sapiens', qw(Hominidae Catarrhini Primates Eutheria Mammalia Vertebrata Chordata Metazoa Eukaryota)
    ]
  };
  
  #Checking explicit returning of the hash
  my $explicit_hash = $helper->execute_into_hash(
    %args,
    -CALLBACK => sub {
      my ($row, $value) = @_;
      if(!$count) {
        ok(! defined $value, 'Checking value is undefined for the first call');
      }
      $value = [] if ! defined $value ;
      push(@{$value}, $row->[1]);
      $count++;
      return $value;
    }
  );
  
  is_deeply($explicit_hash, $expected_hash, 'Checking HASH building allows for callbacks with same data structure');
  
  #Checking when we do an empty return undef
  my $undef_hash = $helper->execute_into_hash(
    %args,
    -CALLBACK => sub {
      my ($row, $value) = @_;
      if(defined $value) {
        push(@{$value}, $row->[1]);
        return;
      }
      my $new_value = [$row->[1]];
      return $new_value;
    }
  );
  
  is_deeply($explicit_hash, $expected_hash, 'Checking HASH building allows for callbacks with same data structure with undef returns');
}

my $zero_count_hash = $helper->execute_into_hash(
  -SQL => 'select 1,0'
);

is($zero_count_hash->{1}, 0, 'Checking hash contains key for zero value');

my $null_count_hash = $helper->execute_into_hash(
  -SQL => 'select 1,NULL'
);

ok(!exists $null_count_hash->{1}, 'Checking hash doesnt contain key for NULL value');

##### CHECKING WORKING WITH A STH DIRECTLY
{
  my $v = $helper->execute_with_sth(-SQL => 'select count(*) from meta', -CALLBACK => sub {
    my ($sth) = @_;
    my $count;
    $sth->bind_col(1, \$count);
    $sth->fetch();
    return $count;
  });
  cmp_ok($v, '>', 0, 'Asserting we found data from meta using execute_with_sth()');
}

#TRANSACTION() CHECKS
my $meta_table_count = $helper->execute_single_result(-SQL => 'select count(*) from meta');
my $meta_memoize = $helper->execute(-SQL => 'select * from meta');

is(scalar(@{$meta_memoize}), $meta_table_count, 'All meta items are returned');

$dba->dbc()->do('alter table meta engine=InnoDB') if $dba->dbc->driver() eq 'mysql';

ok($helper->_perform_transaction_code(), 'This level should do all transaction work');

my $get_value = sub {
  return $helper->execute_single_result(-SQL => 'select meta_value from meta where meta_key =?', -PARAMS => [$meta_key]);
};

{
  #transaction isolation checks
  throws_ok {
    $helper->transaction(-CALLBACK => sub {
      my $sql = 'insert into meta (species_id, meta_key, meta_value) values (?,?,?)';
      $helper->execute_update(
        -SQL => $sql,
        -PARAMS => [2, 'm', '1']
      );
      $helper->execute_update(
        -SQL => $sql,
        -PARAMS => [2, 'm', '2']
      );
      
      my $count = $helper->execute_single_result(-SQL => 'select count(*) from meta where species_id =?', -PARAMS => [2]);
      is($count, 2, 'Count should be 2');
      die 'Dead now';
    });
  }qr/Dead now/, 'Died as expected';
  my $count = $helper->execute_single_result(-SQL => 'select count(*) from meta where species_id =?', -PARAMS => [2]);
  is($count, 0, 'Count should be 0 as we reset the transaction');
}

#Testing multiple level isolation (or more that the framework ignores them)
{
  my $new_meta_value = 'test';
  throws_ok {
  $helper->transaction( -CALLBACK => sub {
    $helper->execute_update(-SQL => 'update meta set meta_value =? where meta_key =?', -PARAMS => [$new_meta_value, $meta_key]);
    eval {
      $helper->transaction(-CALLBACK => sub {
        ok(!$helper->_perform_transaction_code(), 'This level should not be doing any transaction work');
        die 'This will not cause the transaction to be aborted';
      });
    };
    is($get_value->(), $new_meta_value, 'The die from the prior transaction should not have triggered a rollback');
    die('Dead now');
  });
  } qr/Dead now/, 'Expected die found';
  
  isnt($get_value->(), $new_meta_value, 'Meta value is reset as transaction was aborted');
  
  $helper->transaction( -CALLBACK => sub {
    $helper->execute_update(-SQL => 'delete from meta');
  });
  
  $helper->transaction( -CALLBACK => sub {
    $helper->batch(-SQL => 'insert into meta values (?,?,?,?)', -DATA => $meta_memoize);
  });
  
  my $new_count_hash = $helper->execute_into_hash(
    -SQL => 'select meta_key, count(*) from meta group by meta_key'
  );
  is_deeply($new_count_hash, $meta_count_hash, 'Counts of meta keys should be the same');
}


#Testing transactional retry
{
  my $new_meta_value = 'test';
  
  # First try retries until the very last attempt
  {
    my $counter = 0;
    $helper->transaction( -RETRY => 3, -PAUSE => 0.1, -CALLBACK => sub {
      #Die for the first 3 times (so we will succeed on the final attempt)
      $counter++;
      if($counter != 4) {
        die 'Throwing an error to be ignored';
      }
      
      $helper->execute_update(-SQL => 'update meta set meta_value =? where meta_key =?', -PARAMS => [$new_meta_value, $meta_key]);
    });
    is($counter, 4, 'Counter should be set to 4 as we tried 4 attempts at writing (one go & 3 retries)');
    is($get_value->(), $new_meta_value, 'Commit should have gone through after retries');
  }
  
  #Second try will fail as we exhaust our retries
  {
    my $counter = 0;
    throws_ok {
      $helper->transaction( -RETRY => 2, -PAUSE => 0.1, -CALLBACK => sub {
        $counter++;
        die 'Throwing an error 2';
      })
    } qr /Throwing an error 2/, 'Correct error thrown';
    
    is($counter, 3, 'Counter should be set to 3 as we had 3 attempts at writing (one go & 2 retries)');
    is($get_value->(), $new_meta_value, 'Commit should have done nothing');
  }
  
  #Third one says we cannot influence the retry count from a sub-transaction
  {
    my $counter = 0;
    throws_ok {
      $helper->transaction( -RETRY => 1, -PAUSE => 0.1, -CALLBACK => sub {
        $helper->transaction( -RETRY => 10, -CALLBACK => sub {
          $counter++;
          die 'Throwing an error 3';
        });
      })
    } qr /Throwing an error 3/, 'Correct error thrown';
    is($counter, 2, 'Counter should be set to 2 as we had 2 attempts at writing (one go & 1 retry)');
    is($get_value->(), $new_meta_value, 'Commit should have done nothing');
  }

  #Fourth says we only retry when we find a specific issue
  {
    my $counter = 0;
    throws_ok {
      $helper->transaction( 
        -RETRY => 4,
        -PAUSE => 0.1,
        -CALLBACK => sub {
          $counter++;
          die 'fake deadlock' if $counter <= 2;
          die 'Throwing an error 4';
        },
        -CONDITION => sub {
          my ($error) = @_;
          return ( $error =~ /deadlock/ ) ? 1 : 0;
        }
      )
    } qr /Throwing an error 4/, 'Correct error thrown';
    is($counter, 3, 'Counter should be set to 3 as we had 2 fake deadlocks & 1 real error even though we allowed more retries');
  }
  
  #Fith says we sleep for at least the amount we say
  {
    my $counter = 0;
    my $time = time();
    $helper->transaction( -RETRY => 1, -PAUSE => 2, -CALLBACK => sub {
      $counter++;
      if($counter != 2) {
        die 'Throwing an error 5';
      }
      
      $helper->execute_update(-SQL => 'delete from meta where meta_value =? and meta_key =?', -PARAMS => [$new_meta_value, $meta_key]);
    });
    my $elapsed = time() - $time;
    cmp_ok($elapsed, '>=', 2, 'Checking more than 2 seconds elapsed between retries');
    is(
      $helper->execute_single_result(
        -SQL => 'select count(*) from meta where meta_key =? and meta_value=?', 
        -PARAMS => [$meta_key, $new_meta_value]
      ), 0, 
    'Commit will have deleted the meta_key row '.$meta_key);
  }
  
  #Sixth says you cannot repeat the transaction if it worked
  {
    my $counter = 0;
    $helper->transaction( -RETRY => 5, -PAUSE => 0.1, -CALLBACK => sub {
      $helper->execute_single_result('select 1');
      $counter++;
    });
    is($counter, 1, 'Counter should be set to 1 as our transaction was good');
  }
  
  #Reset
  $helper->transaction( -CALLBACK => sub {
    $helper->execute_update(-SQL => 'delete from meta');
    $helper->batch(-SQL => 'insert into meta values (?,?,?,?)', -DATA => $meta_memoize);
  });
}

#Doing hashref checks
{
  my $sql = 'select meta_key, meta_value from meta where meta_key =?';
  my $params = ['species.common_name'];
  {
    my $array_of_hashes = $helper->execute(
      -SQL => $sql,
      -CALLBACK => sub {
        my ($row) = @_;
        return { name => $row->{meta_value} };
      },
      -USE_HASHREFS => 1,
      -PARAMS => $params
    );
    is_deeply($array_of_hashes, [ { name => 'Human' } ], 'HashRefs in a callback works');
  }
  {
    my $array_of_hashes = $helper->execute(
      -SQL => $sql,
      -USE_HASHREFS => 1,
      -PARAMS => $params
    );
    is_deeply($array_of_hashes, [ { meta_key => $params->[0], meta_value => 'Human' } ], 'HashRefs using a default callback works'); 
  }
}

$dba->dbc()->do('alter table meta engine=MyISAM') if $dba->dbc->driver() eq 'mysql';
done_testing();
