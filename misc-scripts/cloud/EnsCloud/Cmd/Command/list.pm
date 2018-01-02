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

package EnsCloud::Cmd::Command::list;
use Moose;
use MooseX::StrictConstructor;

extends qw(MooseX::App::Cmd::Command);
use CHI;

# ABSTRACT: list the application's commands


has 'region_alias' => (
    is     => 'rw',
    isa    => 'Str',
    traits => ['Getopt'],
    #     cmd_aliases   => "h",
    documentation => "asia, useast, uswest or eu",
    required      => 1,
    default => 'useast',
);

has 'refresh' => (
    is => 'ro',
    isa => 'Str',
);


with 'EnsCloud::Describer';
sub execute {
    my ( $self, $opt, $arg ) = @_;
    
#     my $edescriber = EnsCloud::Describer->new( region_alias => $self->region_alias );
    my $cache = CHI->new(
        driver   => 'File',
        root_dir => $ENV{HOME} . '/' . '/ec2cache/' . $self->region_alias,
    );
    my $name           = 'useast_instances';
    my $instance_table = $cache->get($name);
    $DB::single = 1;
    if ( !defined $instance_table || $self->refresh ) {
#}         $instance_table = $edescriber->list_instances;
         $instance_table = $self->list_instances;

        $cache->set( $name, $instance_table, "10 minutes" );
    }
    print $instance_table;
}

sub abstract {

    return 'list cloud things';

}
__PACKAGE__->meta->make_immutable;

1;
__END__

