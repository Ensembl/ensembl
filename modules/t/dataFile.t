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

use strict;
use warnings;

use Test::More;
use Test::Warnings;
use Test::Exception;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::ApiVersion;

use Bio::EnsEMBL::DataFile;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db = $multi->get_DBAdaptor( 'core' );

$multi->hide('core', 'data_file');

my $dfa = $db->get_DataFileAdaptor();
isa_ok($dfa, 'Bio::EnsEMBL::DBSQL::DataFileAdaptor', 'Checking DataFileAdaptor instance');

my $csa = $db->get_CoordSystemAdaptor();
my $aa = $db->get_AnalysisAdaptor();

my $base = '/tmp';
my $expected_base = $base.'/homo_sapiens/NCBI33';
my $cs = $csa->fetch_by_dbID(1);
my $a = $aa->fetch_by_dbID(1);

my %base_args = (
  -COORD_SYSTEM => $cs,
  -ANALYSIS => $a,
  -NAME => 'wibble',
  -VERSION_LOCK => 0,
  -ABSOLUTE => 0,
  -URL => undef,
  -FILE_TYPE => 'BAM',
  -ADAPTOR => $dfa
);

{
  my $df = new_ok('Bio::EnsEMBL::DataFile' => [
    %base_args
  ], 'data file');
  is($df->path($base), $expected_base.'/core/wibble.bam', 'Checking non-absolute path');
  is_deeply($df->get_all_paths($base), [ $expected_base.'/core/wibble.bam', $expected_base.'/core/wibble.bam.bai' ], 'Checking all non-abs paths (BAM)');
}

# now test non-abs paths with bigwig file type
$base_args{-FILE_TYPE} = 'BIGWIG';
{
  my $df = new_ok('Bio::EnsEMBL::DataFile' => [
    %base_args
  ], 'data file');
  is($df->path($base), $expected_base.'/core/wibble.bw', 'Checking non-absolute path');
  is_deeply($df->get_all_paths($base), [ $expected_base.'/core/wibble.bw' ], 'Checking all non-abs paths (BIGWIG)');
}

# now test non-abs paths with bamcov file type
$base_args{-FILE_TYPE} = 'BAMCOV';
{
  my $df = new_ok('Bio::EnsEMBL::DataFile' => [
    %base_args
  ], 'data file');
  is($df->path($base), $expected_base.'/core/wibble.bam', 'Checking non-absolute path');
  is_deeply($df->get_all_paths($base), [ $expected_base.'/core/wibble.bam', $expected_base.'/core/wibble.bam.bai', $expected_base.'/core/wibble.bam.bw'], 'Checking all non-abs paths (BAMCOV)');
}

$base_args{-FILE_TYPE} = 'BAM';

{
  my $df = new_ok('Bio::EnsEMBL::DataFile' => [
    %base_args,
    -VERSION_LOCK => 1,
  ], 'data file');
  
  is($df->path($base), $expected_base.'/'.software_version().'/core/wibble.bam', 'Checking non-absolute version locked path');
}

{
  my $url = 'http://www.google.co.uk/wibble.bam';
  my $df = new_ok('Bio::EnsEMBL::DataFile' => [
    %base_args,
    -ABSOLUTE => 1,
    -URL => $url
  ], 'data file');
  
  is($df->path(), $url, 'Checking absolute path');
}

{
  my $supercontig_cs = $csa->fetch_by_dbID(2);
  my $df = new_ok('Bio::EnsEMBL::DataFile' => [
    %base_args,
    -COORD_SYSTEM => $supercontig_cs
  ], 'data file');
  is($df->path($base), $expected_base.'/core/wibble.bam', 'Checking non-absolute unversioned cs path');
}

{
  my %exts = (BAM => ['bam', 'bam.bai'], BIGWIG => ['bw'], VCF => ['vcf.gz', 'vcf.gz.tbi'], 'BAMCOV' => ['bam', 'bam.bai', 'bam.bw']);
  while( my ($type, $ext) = each %exts ) {
    is_deeply($dfa->DataFile_to_extensions(new_ok('Bio::EnsEMBL::DataFile'=>[%base_args, -FILE_TYPE => $type])), $ext, 'Checking '.$type.' extension');
  }
}

{
  my %local_args = %base_args;
  delete $local_args{-ADAPTOR};
  my $df = new_ok('Bio::EnsEMBL::DataFile' => [ %local_args ], 'data file');
  $dfa->store($df);
  cmp_ok($df->dbID(), '>=', 1, 'Checking it was assigned an ID higher than 1');
  is_deeply($dfa->fetch_by_dbID($df->dbID()), $df, 'Checking retrieved data is the same as what we currently hold');
  
  $df->absolute(1);
  $dfa->update($df);
  is_deeply($dfa->fetch_by_dbID($df->dbID()), $df, 'Checking retrieved data is the same as what we currently hold');
  
  is_deeply($dfa->fetch_all_by_Analysis($a), [$df], 'Checking retrieved data is the same as what we currently hold');
  is_deeply($dfa->fetch_all_by_CoordSystem($cs), [$df], 'Checking retrieved data is the same as what we currently hold');
  is_deeply($dfa->fetch_by_name_and_type('wibble', 'BAM'), $df, 'Checking retrieved data is the same as what we currently hold');    
}

{
  my %type2adaptor = 
    ( 'BAM'    => 'Bio::EnsEMBL::IO::Adaptor::BAMAdaptor',
      'BIGBED' => 'Bio::EnsEMBL::IO::Adaptor::BigBedAdaptor',
      'BIGWIG' => 'Bio::EnsEMBL::IO::Adaptor::BigWigAdaptor',
      'VCF'    => 'Bio::EnsEMBL::IO::Adaptor::VCFAdaptor',
      'BAMCOV' => 1
    );
   
  for my $type (keys %type2adaptor) {
    my $df = new_ok('Bio::EnsEMBL::DataFile'=>[%base_args, -FILE_TYPE => $type]);

    if ($type eq 'BAMCOV') {
      throws_ok { $dfa->DataFile_to_adaptor($df, $base, 'BIGBED') }
	qr/handler found/, 'Request for type incompatible with BAMCOV';
      
      #
      # Cannot actually do the following since this would require
      # a dependency of this test case and hence of the core API
      # on the various external data file adaptors, e.g. BAM/BIGWIG etc.
      #
      # my $eda = $dfa->DataFile_to_adaptor($df, $base, 'BAM');
      # isa_ok($eda, 'Bio::EnsEMBL::IO::Adaptor::BAMAdaptor', 'Checking BAM data adaptor');
      # $eda = $dfa->DataFile_to_adaptor($df, $base, 'BAMCOV');
      # isa_ok($eda, 'Bio::EnsEMBL::IO::Adaptor::BAMAdaptor', 'Checking BAMCOV data adaptor');
      # $eda = $dfa->DataFile_to_adaptor($df, $base, 'BIGWIG');
      # isa_ok($eda, 'Bio::EnsEMBL::IO::Adaptor::BigWigAdaptor', 'Checking BIGWIG data adaptor');

    } else {
      for my $requested_type (keys %type2adaptor) {
	if ($requested_type ne $type) {
	  throws_ok { $dfa->DataFile_to_adaptor($df, $base, $requested_type) }
	    qr/but file is of type/, "Request for $requested_type adaptor, file type $type"
	  } else {
	    #
	    # See the previous comment for the BAMCOV case
	    # 
	    # my $eda = $dfa->DataFile_to_adaptor($df, $base, $requested_type);
	    # isa_ok($eda, $type2adaptor{$type}, "Checking $requested_type data adaptor");
	  }
       
      }
      
    }
  }
  
}

$multi->restore();

done_testing();
