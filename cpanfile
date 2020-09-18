requires 'DBI';
requires 'DBD::mysql';
requires 'HTTP::Tiny';
requires 'IO::Compress::Gzip';
requires 'URI::Escape';
requires 'Config::IniFiles';

test_requires 'Test::Warnings';
test_requires 'Test::Differences';
test_requires 'Test::Exception', '>0.42';
test_requires 'Test::MockObject';
test_requires 'Test::Deep';
test_requires 'Test::More';
test_requires 'Devel::Peek';
test_requires 'Devel::Cycle';
test_requires 'Error';
test_requires 'PadWalker';
test_requires 'Test::Builder::Module';
test_requires 'IO::String';
test_requires 'Test::Perl::Critic';
test_requires 'Test::Pod::Coverage';
test_requires 'Perl::Critic::Utils';

feature 'assembly_mapping', 'Assembly mapper support' => sub {
  requires 'Algorithm::Diff';
  requires 'Tie::IxHash';
};

feature 'xref_mapping', 'Xref mapping pipeline' => sub {
  requires 'DBIx::Class';
  requires 'Digest::MD5';
  requires 'LWP::UserAgent';
  requires 'List::Util', '>= 1.45';
  requires 'Moose';
  requires 'SQL::Translator', '>= 0.11018';
  requires 'Text::CSV';
  recommends 'Text::CSV_XS';
  requires 'Text::Glob';
  requires 'URI';
  requires 'XML::LibXML';

  test_requires 'Config::General';
  test_requires 'Perl::Critic::Moose';
};
