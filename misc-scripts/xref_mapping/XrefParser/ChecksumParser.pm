=head1 LICENSE

See the NOTICE file distributed with this work for additional information
regarding copyright ownership.

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

package XrefParser::ChecksumParser;

# Input format looks like:
#
# UPI0001B45C00 71B80D7A684B1F2DEDDA7B5AEE1D029E
# UPI0002473BEA 4542D97F3AB3F7B656ABB941AED3F2BB
# UPI00024743AF A69E7EEE820CA54100AD43E86BE823E4

use strict;
use warnings;
use Carp;
use English qw( -no_match_vars );
use IO::File;
use base qw( XrefParser::BaseParser );

my $TABLE_NAME = 'checksum_xref';

sub run {
  my ($self, $ref_arg) = @_;
  
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $files        = $ref_arg->{files};
  my $verbose      = $ref_arg->{verbose};
  
  if((!defined $source_id) or (!defined $species_id) or (!defined $files) ){
    croak "Need to pass source_id, species_id and files as pairs";
  }
  $verbose ||=0;

  # FIXME: this will fail if the input file is in a read-only directory (ENSCORESW-3197)
  my $target_file = $files->[0].'.mysqlinput';
  my $input_fh = $self->get_filehandle($files->[0]);
  if(-f $target_file) {
    print "Target file '${target_file}' already exists; removing\n" if $verbose;
    unlink $target_file;
  }
  my $output_fh = IO::File->new($target_file, 'w')
    || croak "Failed to open ${target_file} for writing: ${OS_ERROR}";
  
  $self->_transfer_contents($input_fh, $output_fh, $source_id);
  
  close($input_fh);
  close($output_fh);
  
  $self->_load_table($target_file, $verbose, $source_id);
  
  return;
}

sub _transfer_contents {
  my ($self, $input_fh, $output_fh, $source_id) = @_;
  my $dbh = $self->dbi();
  my ($counter) = $dbh->selectrow_array('select max(checksum_xref_id) from '.$TABLE_NAME );
  while(my $line = <$input_fh>) {
    chomp $line;
    my ($upi, $checksum) = split(/\s+/, $line);

    # Use an ID one higher than the last. Obvious? Perhaps - except before
    # the commit adding this comment the code only incremented $counter
    # AFTER using it.
    $counter += 1;

    my @output = ($counter, $source_id, $upi, $checksum);
    print $output_fh join("\t", @output);
    print $output_fh "\n";
  }
  return;
}

sub _load_table {
  my ($self, $file, $verbose, $source_id) = @_;
  my $dbh = $self->dbi();
  my ($count) = $dbh->selectrow_array('select count(*) from '.$TABLE_NAME . ' WHERE source_id = ' . $source_id);
  if($count) {
    print "'$TABLE_NAME' has rows for $source_id; deleting\n" if $verbose;
    $dbh->do('delete from ' . $TABLE_NAME . ' WHERE source_id = ' . $source_id);
  }
  print "Loading data into '$TABLE_NAME' from '$file'\n" if $verbose;
  my $load = sprintf(q{LOAD DATA LOCAL INFILE '%s' INTO TABLE %s}, $file, $TABLE_NAME);
  $dbh->do($load);
  print "Finished loading data into '$TABLE_NAME'\n" if $verbose;
  return;
}

1;
