requires 'DBI';
requires 'DBD::mysql';
requires 'HTTP::Tiny';
requires 'IO::Compress::Gzip';
requires 'URI::Escape';

test_requires 'Test::Warnings';
test_requires 'Test::Differences';
test_requires 'Test::Exception';
test_requires 'Test::MockObject';
test_requires 'Test::More';
test_requires 'Devel::Peek';
test_requires 'Devel::Cycle';
test_requires 'Error';
test_requires 'PadWalker';
test_requires 'Test::Builder::Module';
test_requires 'IO::String';
test_requires 'Test::Perl::Critic';
test_requires 'Perl::Critic::Utils';

=cut
feature 'assembly_mapping', 'Assembly mapper support' => sub {
  requires 'Algorithm::Diff';
  requires 'Tie::IxHash';
};

feature 'xref_mapping', 'Xref mapping pipeline' => sub {
  requires 'Config::IniFiles';
  requires 'Digest::MD5';
  requires 'Text::Glob';
  requires 'XML::LibXML';
  requires 'XML::Simple';
};
=cut
