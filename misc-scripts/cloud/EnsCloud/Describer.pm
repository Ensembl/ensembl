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

package EnsCloud::Describer;

use Moose::Role;
requires 'region_alias';

use Net::Amazon::EC2;
use Text::SimpleTable;
our $col_widths_headings = {
    tags          => [ 22, 'Tag' ],
    instance_id   => [ 10, 'Instance' ],
    ip_address    => [ 14, 'Ip Address' ],
    dns_name      => [ 57, 'DNS Name' ],
    instance_type => [ 9,  'Type' ],
    launch_time   => [ 16, 'Launch' ],
    image_id      => [ 12, 'Ami id' ],
    name          => [ 45, 'Name' ],
    description   => [ 45, 'Description' ],
    image_state   => [ 7,  'State' ],

    # snap_creation =>    [ 23, 'Snap Creation' ],
    #     [ 13, 'Vol Creation' ]
    size        => [ 5,  'Size' ],
    create_time => [ 20, 'Create Time' ],

    #     [ 20, 'Attach Time' ],
    #     [ 13, 'Volume Id' ],
    #     [ 10, 'Device' ],
    #     [ 10, 'Vol. status' ],
    #     [ 8,  'Attach status' ],
    # zone =>    [ 10, 'Zone' ]
    #     [ 13, 'Snapshot Id' ],
    #     [ 19, 'Start Time' ],
    #     [ 8,  'Progress' ],
    # status =>    [ 10, 'Status' ],
    #     [ 12, 'owner' ],
    #     [ 11, 'owner alias' ],
};

has ec2_wanted_version => (
    is      => 'rw',
    default => '2011-02-28'

);

#has 'ec2' => (
#    is         => 'ro',
#    isa        => 'Net::Amazon::EC2',
#lazy_
#);

# has 'region_alias' => (
#     is  => 'rw',
#     isa => 'Str',
#     required =>1
# );

sub ec2 {
    my $self          = shift;
    my $region_lookup = {
        'eu'     => 'eu-west-1',
        'useast' => 'us-east-1',
        'uswest' => 'us-west-1',
        'asia'   => 'ap-southeast-1',
    };
    return Net::Amazon::EC2->new(
        AWSAccessKeyId  => $ENV{AWS_ACCESS_KEY_ID},
        SecretAccessKey => $ENV{AWS_ACCESS_KEY_SECRET},
        debug           => 0,
        region          => $region_lookup->{ $self->region_alias },
        version         => $self->ec2_wanted_version
    );

}

sub list_volumes {

    my $self           = shift;
    my $volumes        = $self->ec2->describe_volumes;
    my @sorted_volumes = sort { $a->{size} <=> $b->{size} } @$volumes;

    my $instances = $self->ec2->describe_instances;

    my $instance2ip;

    foreach my $instance (@$instances) {

        #         print dump $instance->instances_set;

        my $instance_set = $instance->instances_set->[0];
        $instance2ip->{ $instance_set->{instance_id} } =
          $instance_set->{ip_address};

    }

    #     die dump $volumes;
    my $table = Text::SimpleTable->new(
        [ 5,  'Size' ],
        [ 20, 'Create Time' ],
        [ 20, 'Attach Time' ],
        [ 13, 'Volume Id' ],
        [ 13, 'Instance Id' ],
        [ 15, 'Ip Address' ],
        [ 10, 'Device' ],
        [ 10, 'Vol. status' ],
        [ 8,  'Attach status' ],
        [ 13, 'Snapshot Id' ],
        [ 10, 'Zone' ]
    );

    foreach my $volume (@sorted_volumes) {

        #         print dump $instance->instances_set;

        my $attachments = $volume->attachments;

        print
"Volume [$volume->{volume_id}] is attached to multiple instances: this not handled by this script\n"
          and exit
          if @$attachments > 1;
        $volume->{create_time} =~ s/.000Z$//;
        $attachments->[0]->{attach_time} =~ s/.000Z$//;

        #         die dump $instance_set;
        $table->row(
            $volume->{size},
            $volume->{create_time},
            $attachments->[0]->{attach_time} || '',
            $volume->{volume_id},
            $attachments->[0]->{instance_id}                   || '',
            $instance2ip->{ $attachments->[0]->{instance_id} } || '',
            $attachments->[0]->{device}                        || '',
            $volume->{status},
            $attachments->[0]->{status} || '',

            $volume->{snapshot_id} || 'no snapshot',
            $volume->{zone}
        );
    }

    print $table->draw;

}

sub list_instances {
    my ($self) = @_;

    my $instances = $self->ec2->describe_instances;
    my $table     = Text::SimpleTable->new(
        $col_widths_headings->{tags},
        $col_widths_headings->{instance_id},
        $col_widths_headings->{ip_address},
        $col_widths_headings->{dns_name},
        $col_widths_headings->{instance_type},
        $col_widths_headings->{launch_time},
        $col_widths_headings->{image_state},
        $col_widths_headings->{image_id}
    );
    foreach my $instance (@$instances) {

        #    print dump $instance->instances_set;

        my $instance_set = $instance->instances_set->[0];
        $instance_set->{launch_time} =~ s/:\d+\.000Z$//;

        #die dump $instance_set;
        $table->row(
            $instance_set->{tags} || '',
            $instance_set->{instance_id},
            $instance_set->{ip_address} || '',
            $instance_set->{dns_name}   || '',
            $instance_set->{instance_type},
            $instance_set->{launch_time},
            $instance_set->{instance_state}->name,
            $instance_set->{image_id},
        );
    }
    return $table->draw;
}

# __PACKAGE__->meta->make_immutable;

1;

__END__
