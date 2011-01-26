#A set tests used to prod the SqlHelper class

use strict;
use warnings;

use Test::More;
use Test::Exception;
use Scalar::Util qw(isweak);

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Utils::SqlHelper;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $dba = $multi->get_DBAdaptor( 'core' );
ok( $dba, 'Test database instatiated' );

#Now start testing the Helper
dies_ok { Bio::EnsEMBL::Utils::SqlHelper->new(-DB_CONNECTION => $dba) } 
  'Expect to die when we do not give SqlHelper a DBConncetion'; #was given a DBAdaptor
ok ( 
  isweak(Bio::EnsEMBL::Utils::SqlHelper->new(-DB_CONNECTION => $dba->dbc())->{db_connection}),
  'Checking DBConnection reference is weak when we ask for it' 
);

my $helper = Bio::EnsEMBL::Utils::SqlHelper->new(-DB_CONNECTION => $dba->dbc());
ok ( $helper, 'SqlHelper instance was created' );


my $meta_key = 'species.common_name';
diag("Meta key queries working with ${meta_key}. If the tests fail then check for it in the DB dumps");

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
      qw(sapiens Homo Hominidae Catarrhini Primates Eutheria Mammalia Vertebrata Chordata Metazoa Eukaryota)
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


#TRANSACTION() CHECKS
my $meta_table_count = $helper->execute_single_result(-SQL => 'select count(*) from meta');
my $meta_memoize = $helper->execute(-SQL => 'select * from meta');

is(scalar(@{$meta_memoize}), $meta_table_count, 'All meta items are returned');

$dba->dbc()->do('alter table meta engine=InnoDB');

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

#Doing hashref checks
{
  my $sql = 'select meta_key, meta_value from meta where meta_key =?';
  my $callback = sub {
    my ($row) = @_;
    return { name => $row->{meta_value} };
  };
  my $array_of_hashes = $helper->execute(
    -SQL => $sql,
    -CALLBACK => $callback,
    -USE_HASHREFS => 1,
    -PARAMS => ['species.common_name']
  );
  is_deeply($array_of_hashes, [ { name => 'Human' } ], 'HashRefs in a callback works');
}

$dba->dbc()->do('alter table meta engine=MyISAM');
done_testing();
