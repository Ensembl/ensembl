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
use Test::MockObject;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils qw/mock_object/;

# "Globals"
my $CHR           = '20';
my $START         = 30_270_000;
my $END           = 30_270_200;
my $STRAND        = 1;
my $LENGTH				= 201;	
my $REGISTRY			= 'Bio::EnsEMBL::Registry';

# Get a human core DB so we have something in the registy to mess about with
my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $dba = $multi->get_DBAdaptor('core');
my $empty_dba = $multi->get_DBAdaptor('empty');

sub mock_adaptor {
	my ($residue) = @_;
	my $mock_sequence_adaptor = Test::MockObject->new();
	$mock_sequence_adaptor->mock('fetch_by_Slice_start_end_strand', sub {
		my ( $self, $slice, $start, $end, $strand ) = @_;
		my $seq = $residue x $slice->length(); # return a bogus piece of sequence
		return \$seq;
	});
}

# Create the mock sequence adaptor which responds to one method
my $m_mock_adaptor = mock_adaptor('M');
my $x_mock_adaptor = mock_adaptor('X');

#Get the slice & setup assertion methods
my $slice = $dba->get_SliceAdaptor()->fetch_by_region('toplevel', $CHR, $START, $END);
sub normal_slice {
	my ($local_slice) = @_;
	note 'Asserting normal slice content';
	$local_slice ||= $slice;
	my $seq = $local_slice->seq();
	is(length($seq), $LENGTH, 'Expected retrieved sequence length');
	like($seq, qr/^[ACTG]+$/, 'Checking we just have ACTG');
}

sub m_switchable_slice {
	my ($local_slice) = @_;
	note 'Asserting mocked slice content of just M';
	$local_slice ||= $slice;
	my $seq = $local_slice->seq();
	is(length($seq), $LENGTH, 'Expected retrieved sequence length');
	like($seq, qr/^M+$/, 'Checking that we use the mock adaptor and it only reports the residue M');
}

sub x_switchable_slice {
	my ($local_slice) = @_;
	note 'Asserting mocked slice content of just X';
	$local_slice ||= $slice;
	my $seq = $local_slice->seq();
	is(length($seq), $LENGTH, 'Expected retrieved sequence length');
	like($seq, qr/^X+$/, 'Checking that we use the mock adaptor and it only reports the residue X');
}

#Test normality
normal_slice();

#Replace the adaptor & see what occurs
$REGISTRY->add_switchable_adaptor($dba->species(), $dba->group(), 'sequence', $m_mock_adaptor);
ok($REGISTRY->has_switchable_adaptor($dba->species(), $dba->group(), 'sequence'), 'Has a sequence switchable adaptor should respond with true');
ok(!$REGISTRY->has_switchable_adaptor($dba->species(), $dba->group(), 'someotheradaptor'), 'Has switchable adaptor should respond with false for a bogus adaptor');
m_switchable_slice();

#Attempt again and catch exceptions
dies_ok { $REGISTRY->add_switchable_adaptor($dba->species, $dba->group, 'sequence', $x_mock_adaptor) } 'Dies if we try to set it without resetting';
dies_ok { $REGISTRY->add_switchable_adaptor($dba->species, $dba->group, 'sequence', $m_mock_adaptor) } 'Dies if we try to set it to the same adaptor without resetting';

#Now switch with a force and assert
$REGISTRY->add_switchable_adaptor($dba->species, $dba->group, 'sequence', $x_mock_adaptor, 1);
x_switchable_slice();

#Remove and ensure it goes back to normal
$REGISTRY->remove_switchable_adaptor($dba->species, $dba->group, 'sequence');
normal_slice();
ok(!$REGISTRY->has_switchable_adaptor($dba->species(), $dba->group(), 'sequence'), 'Has switchable adaptor should respond with true now we have removed the adaptor');

###############################################################
##### Now test the expected functionality of the delegating 
##### DNADB interface
###############################################################
is($empty_dba->group(), 'empty', 'Make sure empty DB group is empty');
is($empty_dba->dnadb(), $empty_dba, 'Empty DNADB is itself since it is an unknown group to ConfigRegistry');

# If we call normal_slice on the empty DBA then we should get the original sequence back since it's a duplication
my $empty_dba_slice = $empty_dba->get_SliceAdaptor()->fetch_by_region('toplevel', $CHR, $START, $END);
normal_slice($empty_dba_slice);

# Set the dnadb to core
$empty_dba->dnadb($dba);

#Set the switchable adaptor on both empty dba and dba. empty dba should win
$REGISTRY->add_switchable_adaptor($dba->species, $dba->group, 'sequence', $m_mock_adaptor);
$REGISTRY->add_switchable_adaptor($empty_dba->species, $empty_dba->group, 'sequence', $x_mock_adaptor);
# my $empty_dba_slice = $empty_dba->get_SliceAdaptor()->fetch_by_region('toplevel', $CHR, $START, $END);
x_switchable_slice($empty_dba_slice);
# M should now win
$REGISTRY->remove_switchable_adaptor($empty_dba->species, $empty_dba->group, 'sequence');
m_switchable_slice($empty_dba_slice);
$REGISTRY->remove_switchable_adaptor($dba->species, $dba->group, 'sequence');

###############################################################
##### Test the methods on DBAdaptor for switching including 
##### automatic re-switching and remembering last switched
##### adaptors
###############################################################
note 'Running DBAdaptor level tests';
normal_slice();
ok(!$dba->has_switched_adaptor('sequence'), 'Verifying the instance is not being switched out');
$dba->switch_adaptor('sequence', $m_mock_adaptor);
ok($dba->has_switched_adaptor('sequence'), 'Verifying the instance is being switched out');
note 'Try switching to X. Will fail';
m_switchable_slice();
dies_ok { $dba->switch_adaptor('sequence', $x_mock_adaptor); } 'Forcing a switch to an already switched adaptor causes an error';
note 'Still set to M';
m_switchable_slice();
$dba->revert_adaptor();
normal_slice();

note 'Running callback methods';
$dba->switch_adaptor('sequence', $m_mock_adaptor, sub {
	m_switchable_slice();
});
normal_slice();

note 'Forcing a switch using callbacks';
$dba->switch_adaptor('sequence', $m_mock_adaptor);
$dba->switch_adaptor('sequence', $x_mock_adaptor, sub {
	x_switchable_slice();
}, 1);
normal_slice();

note 'Trying a switch callback and die. The revert method should still fire';

dies_ok {$dba->switch_adaptor('sequence', $m_mock_adaptor, sub {
	die "Raising an exception";
});
} 'Dies are propagated through correctly';
ok(! $dba->has_switched_adaptor('sequence'), 'We have reverted the switch');
normal_slice();



done_testing();
