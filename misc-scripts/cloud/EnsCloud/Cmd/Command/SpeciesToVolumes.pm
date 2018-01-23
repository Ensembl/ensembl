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

package EnsCloud::Cmd::Command::SpeciesToVolumes;
use Moose;
use autodie qw(:all);
use Data::Dump qw(dump);
use Data::Dumper;
use File::Spec;
use Moose::Util::TypeConstraints;
use EnsCloud::Image;
use EnsCloud::Image::VolumeBundle;
use EnsCloud::Image::VolumeBundle::Volume;
use EnsCloud::Image::VolumeBundle::Volume::DatabaseDetails;
use namespace::autoclean;
use IPC::Cmd qw[can_run run];

# use List::Utils qw(first);
use Log::Log4perl qw(:easy);
with 'MooseX::Log::Log4perl';
extends qw(MooseX::App::Cmd::Command);

BEGIN {
    Log::Log4perl->easy_init();
}

# ABSTRACT: list the application's commands
sub abstract {

	return "build ensembl MySQL instances by species";

    }


# Instructions
#
# * This expects a volume containing the Ensembl MYDs attached to the instance (Public Snapshot of 65 dbs = snap-56c9ab32)
# - default path is /vols/ensembl_mysql_data
# Run like this: ( N.B --base_snapshot is the base AMI with the OS)
#
# ecloud speciestovolumes   --base_snapshot snap-e36fde86  --species saccharomyces_cerevisiae (omit --species to do all)
# or by dbtype
# ecloud speciestovolumes   --base_snapshot snap-e36fde86  --species homo_sapiens --dbtype variation
# 
# * Let it run until completion.  then wait until all PENDING snapshots are complete.
# Once all PENDING snapshots are complete you will have one snapshot per species
# (Do not start process again while there are are still PENDING snapshots) 
# Now run it again with the same parameters - this will loop the completed snapshots to combine each with the base_image OS 

# has 'type' => (
#     is     => 'rw',
#     isa    => 'Str',
#     traits => ['Getopt'],
#
#     #     cmd_aliases   => "h",
#     documentation => "instances volumes or images",
#     required      => 1,
#
#     #     default       => sub { die "bucket name required" },
# );

has 'region_alias' => (
    is     => 'rw',
    isa    => 'Str',
    traits => ['Getopt'],

    #     cmd_aliases   => "h",
    documentation => "asia, useast, uswest or eu",
    required      => 1,
    default       => 'useast',
);

has 'volume_path' => (

    is     => 'ro',
    isa    => 'Str',
    traits => ['Getopt'],

    #     cmd_aliases   => "h",
    documentation => "path to the species DBs",
    required      => 1,
    default       => '/vols/ensembl_mysql_data',

);

enum 'CompositeGroup' => qw(core coreplusvariation corepluscompara variation compara all);
has 'dbtypes'         => (

    is     => 'ro',
    isa    => 'CompositeGroup',
    traits => ['Getopt'],

    #     cmd_aliases   => "h",
    documentation => "database types to copy: default=core",
    required      => 1,
    default       => 'all',

);

has db_type_lookup => (
    is       => 'ro',
    isa      => 'HashRef',
    required => 1,
    default  => sub {
        {
            core => [qw(core coreexpressionatlas coreexpressionest funcgen otherfeatures rnaseq vega )],
            coreplusvariation =>
              [qw(core coreexpressionatlas coreexpressionest funcgen otherfeatures rnaseq vega variation)],
            corepluscompara =>
              [qw(core coreexpressionatlas coreexpressionest funcgen otherfeatures rnaseq vega compara )],
            all =>
              [qw(core coreexpressionatlas coreexpressionest funcgen otherfeatures rnaseq vega compara variation )],
            variation => [qw(variation)],
            compara   => [qw(compara)],
        };
    },

);

has 'species' => (
    is     => 'ro',
    isa    => 'Str',
    traits => ['Getopt'],

    #     cmd_aliases   => "h",
    documentation => "the species to copy. default=all",
    required      => 1,
    default       => 'all',

);

has 'start_from' => (
    is     => 'ro',
    isa    => 'Str',
    traits => ['Getopt'],

    #     cmd_aliases   => "h",
    documentation => "the species to start dumping from",
    default       => '.*'
);

has 'end_at' => (
    is     => 'ro',
    isa    => 'Str',
    traits => ['Getopt'],

    #     cmd_aliases   => "h",
    documentation => "the species to stop dumping at",
    default       => '^$'
);
has 'volume_dump_queue' => (
    traits  => ['Array'],
    handles => {
        all_queued_volumes      => 'elements',
        volume_add_to_queue     => 'push',
        next_volume_from_queue  => 'shift',
        sort_queue_by_size_desc => [ sort_in_place => ( sub { $_[1]->total_size <=> $_[0]->total_size } ) ],
        find_tag                => 'first',
        queue_length            => 'count',
        filter_queue            => 'grep'
    },
    isa => 'ArrayRef[EnsCloud::Image::VolumeBundle::Volume]',
);

has db_species_details => (
    is => 'rw',

    #isa => 'HashRef',
    isa     => 'HashRef[EnsCloud:::DatabaseDetails]',
    default => sub { {} }
);

has this_instance_id => (
    is       => 'ro',
    isa      => 'Str',
    default  => sub { return `curl -s 169.254.169.254/latest/meta-data/instance-id` },
    required => 1

);

has myd_destination_folder => (
    is       => 'ro',
    isa      => 'Str',
    default  => '/vols/MYDCOPY',
    required => 1

);
has device => (
    is       => 'ro',
    isa      => 'Str',
    default  => '/dev/sdi',
    required => 1

);
has base_snapshot => (
    is       => 'ro',
    isa      => 'Str',
    required => 1

);
has this_zone => (
    is       => 'ro',
    isa      => 'Str',
    default  => sub { return `curl -s 169.254.169.254/latest/meta-data/placement/availability-zone/` },
    required => 1

);
with 'EnsCloud::Describer';

sub execute {
    my ($self) = @_;
    $self->build_queue();

    $self->log->fatal("No databases in the queue");
    my @image_list = ();

    $self->sort_queue_by_size_desc;

    foreach my $volume ( $self->all_queued_volumes ) {
        $self->log->info( "Doing: " . $volume->species );
        $self->update_queue;
        if ( $volume->snapshot_id && $volume->status eq 'completed' ) {

            #  $self->log->info($volume->snapshot_id): $self->log->warn("No snapshot for: " . $volume->tag);
            my $image = $self->make_image($volume);
            push @image_list, $image;
        }
        else {
            $self->make_snapshot($volume);
        }

    }
    $self->update_queue;
    my $image_count = @image_list > 0 ? scalar @image_list : 0;
    $self->log->info( "Finished, made ", $self->queue_length, " Volumes" );
    $self->log->info("Finished, made $image_count Images");
    my @no_snapshots = $self->filter_queue( sub { $_->status =~ /no snapshot/ } );
    $self->log->info( scalar @no_snapshots . " Snapshots still waiting" );
}

sub update_queue {
    my ($self) = @_;

    my $ec2_snapshots = $self->ec2->describe_snapshots( Owner => 'self' );
    foreach my $ec2_snap (@$ec2_snapshots) {
        my $description = $ec2_snap->{description} || next;
        $description =~ s/\s+$//;
        if ( my $has_snapshot = $self->find_tag( sub { $_->tag eq $description } ) ) {
            $has_snapshot->snapshot_id( $ec2_snap->{snapshot_id} );
            $has_snapshot->status( $ec2_snap->{status} );
        }

    }

}

sub make_image {
    my ( $self, $snapshot ) = @_;

    $self->log->debug("MAKE IMAGE");

    my $tag =
        "Ensembl"
      . $snapshot->ensembl_release . " "
      . $snapshot->species . " ["
      . ( join " ", map { $_->type } $snapshot->all_dbs ) . "]";
    my $name = "Ensembl" . $snapshot->ensembl_release . " " . $snapshot->species . " MySQL AMI";

    # has it already been created
    my $existing_images = $self->ec2->describe_images( Owner => 'self' );
    foreach my $i (@$existing_images) {
        if ( $i->{description} eq $tag ) {
            $self->log->warn( "Image " . $i->{image_id} . " already exists for " . $tag );
            return undef;
        }
    }
    my $create_image = $self->ec2->register_image(
        Name               => $name,
        Description        => $tag,
        Architecture       => 'x86_64',
        KernelId           => 'aki-427d952b',
        RootDeviceName     => '/dev/sda1',
        BlockDeviceMapping => [
            {
                deviceName => '/dev/sda1',
                ebs        => {
                    snapshotId          => $self->base_snapshot,
                    deleteOnTermination => 'true'
                }
            },
            {
                deviceName => '/dev/sdh',
                ebs        => { snapshotId => $snapshot->{snapshot_id}, deleteOnTermination => 'true' }
            },
        ]
    );
    if ( $create_image->can('errors') ) {
        $self->log->error( "[Error Creating Image] " . $self->pp_ec2_errors( $create_image->errors ) );
        return;
    }
    $self->log->info("Created Image $create_image for $tag");
    return $create_image;
}

sub make_snapshot {
    my ( $self, $bag_of_dbs ) = @_;
    my $size_as_float = $bag_of_dbs->total_size;
    $self->log->debug( "Making SNAPSHOT " . $bag_of_dbs->tag . " Current status=" . $bag_of_dbs->status );

    # round up volumesize to nearest Gb +1
    # rounding to nearest Gig doesn't seem to enough
    my $round_up_size = int( $size_as_float + 2 );

    #$self->clear_ec2;
    my $volume = $self->ec2->create_volume( Size => $round_up_size, AvailabilityZone => $self->this_zone );

    if ( $volume->can('errors') ) {
        $self->log->error( "[Error Creating Volume] " . $self->pp_ec2_errors( $volume->errors ) );
        return;
    }
    else {
        $self->log->info( "Created ", $volume->volume_id, " ", $volume->size, " Gb" );
    }
    my $do_attach = $self->ec2->attach_volume(
        VolumeId   => $volume->volume_id,
        InstanceId => $self->this_instance_id,
        Device     => $self->device
    );
    $self->log->info( "Attaching " . $volume->volume_id . " Device " . $self->device );
    if ( $do_attach->can('errors') ) {
        $self->log->error( "[Error attaching Volume] " . $self->pp_ec2_errors( $do_attach->errors ) );
        $self->log->info( "Deleting " . $volume->volume_id );
        $self->ec2->delete_volume( VolumeId => $volume->volume_id );
        return;
    }
    else {
        my $wait_time     = 0;
        my $attach_status = '';
        while ( $attach_status ne 'attached' && $wait_time < 30 ) {

            $self->log->info("Waiting 10 seconds for volume to become available");
            sleep 10;
            $wait_time += 10;
            $attach_status =
              $self->ec2->describe_volumes( VolumeId => $volume->volume_id )->[0]->attachments->[0]->{status};
            $self->log->info("Volume $attach_status");
        }
        unless ( $attach_status eq 'attached' ) {
            $self->log->error("ERROR ATTACHING VOLUME AFTER 30 seconds... DETACHING AND DELETING");
            $self->clean_up_by_volume($volume);

            # todo, check for errors in the delete call
            $self->log->fatal("COULD NOT ATTACH VOLUME. THIS IS BAD. EXITING") && die;
        }

    }

    # Make filesystem on the new volume
    $self->log->info( "Making filesystem on " . $self->device );
    my $mkfs_path = can_run('mkfs.xfs') or $self->log->warn('mkfs.xfs not installed!');
    my $mkfs_cmd = [ 'sudo', $mkfs_path, $self->device ];
    unless ( $self->run_command($mkfs_cmd) ) {
        $self->log->info("Cleaning up after failed mkfs");
        $self->clean_up_by_volume($volume);
        return;

    }

    # Make the directory
    my $mkdir_cmd = [ 'sudo', 'mkdir', '-p', $self->myd_destination_folder ];
    unless ( $self->run_command($mkdir_cmd) ) {
        $self->log->info("Cleaning up after failed mkdir");
        $self->clean_up_by_volume($volume);
        return;

    }

    # Mount new volume to it
    $self->log->info( "Mounting " . $self->device );
    my $mount_path = can_run('mount') or $self->log->warn('mount not installed!');
    my $mount_cmd = [ 'sudo', $mount_path, $self->device, $self->myd_destination_folder ];
    unless ( $self->run_command($mount_cmd) ) {
        $self->log->info("Cleaning up after failed mount");
        $self->clean_up_by_volume($volume);
        return;
    }

    # Copy the each MYD dir
    foreach my $myd_dir ( $bag_of_dbs->all_dbs ) {
        my $copy_cmd = [ 'sudo', 'cp', '-r', $myd_dir->myd_path, $self->myd_destination_folder ];
        $self->log->debug( join " ", @$copy_cmd );
        unless ( $self->run_command($copy_cmd) ) {

            $self->log->info("Cleaning up after failed copy");
            my $umount_path = can_run('umount') or $self->log->warn('umount not installed!');
            my $umount_cmd = [ 'sudo', $umount_path, $self->device ];
            unless ( $self->run_command($umount_cmd) ) {
                $self->log->fatal( "CANNOT UMOUNT " . $self->device . "EXITING" ) && die;
            }
            $self->clean_up_by_volume($volume);
            return;

        }
        $myd_dir->is_copied(1);

        #         $snapshot_description .= $myd_dir->name . " ";
    }

    $self->log->info( "umounting " . $self->device );
    my $umount_path = can_run('umount') or $self->log->warn('umount not installed!');
    my $umount_cmd = [ 'sudo', $umount_path, $self->device ];
    return unless $self->run_command($umount_cmd);

    my $wait_time = 0;
    $self->ec2->detach_volume( VolumeId => $volume->volume_id );
    while ( defined eval { $self->ec2->describe_volumes( VolumeId => $volume->volume_id )->[0]->attachments }
        && $wait_time < 60 )
    {
        $self->log->info("Waiting 10 seconds for volume to detach");
        sleep 10;
        $wait_time += 10;
    }

    my $snapshot_description = $bag_of_dbs->{tag};
    $self->log->info("Creating Snapshot");
    my $snapshot = $self->ec2->create_snapshot( VolumeId => $volume->volume_id, Description => $snapshot_description );

    if ( $snapshot->can('errors') ) {
        $self->log->error( "[Snapshot Creation Error] " . $self->pp_ec2_errors( $snapshot->errors ) );
        return;
    }
    else {
        $self->log->info( "Created Snapshot ", $snapshot->snapshot_id, " ", $volume->size, " Gb" );
        $self->log->info( "Deleting " . $volume->volume_id );
        $self->ec2->delete_volume( VolumeId => $volume->volume_id );
        $self->log->info("Tagging the Snapshot");
        my $tag_path = can_run('ec2-create-tags') or $self->log->warn('ec2-describe-tags not found');
        my $tag_cmd = [ $tag_path, $snapshot->snapshot_id, '-t', "Name=$snapshot_description" ];
        return unless $self->run_command($tag_cmd);

        # $self->log->info("Detaching Temp Volume");
        # $self->ec2->detach_volume(VolumeId => $volume->volume_id);
        # $self->log->info("Deleting Temp Volume");
        # $self->ec2->delete_volume(VolumeId => $volume->volume_id);

    }

}

sub clean_up_by_volume {
    my ( $self, $volume ) = @_;

    #detach volume
    $self->log->info( "CLEANUP: Detaching " . $volume->volume_id );
    my $wait_time = 0;
    my $do_detach = $self->ec2->detach_volume( VolumeId => $volume->volume_id );

    if ( $do_detach->can('errors') ) {
        $self->log->error( "[Error detaching Volume] " . $self->pp_ec2_errors( $do_detach->errors ) );

        $self->log->fatal( "CANNOT DETACH VOLUME DURING CLEANUP. THIS IS BAD. EXITING :" . $volume->volume_id )
          && die;
    }

    while (
        defined eval { $self->ec2->describe_volumes( VolumeId => $volume->volume_id )->[0]->attachments }

        && $wait_time < 60
      )
    {

        #                 $self->log->info("Volume $attach_status");

        $self->log->info("Waiting 10 seconds for volume to detach");
        sleep 10;
        $wait_time += 10;
    }

    # delete it
    $self->log->info( "CLEANUP: Deleting " . $volume->volume_id );
    $self->ec2->delete_volume( VolumeId => $volume->volume_id );
    return;

}

sub pp_ec2_errors {
    my ( $self, $error_obj ) = @_;
    my $full_message;
    foreach my $error (@$error_obj) {
        $full_message .= $error->message . "\n";

    }
    return $full_message;
}

sub run_command {
    my ( $self, $cmd ) = @_;

    ### in list context ###
    my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) = run( command => $cmd, verbose => 0 );

    if ($success) {

        $self->log->info( join " ", @$cmd, " Successful" );

        #        return $stdout_buf;
        #         print " is what the command printed:\n";
        my $stdout_str = join "", @$stdout_buf;
        return length $stdout_str > 0 ? $stdout_str : $success;
    }
    else {

        $self->log->error( "Failed!!: $error_message\n" . ( join "", @$stderr_buf ) );
        return 0;
    }

}

sub build_queue {
    my ($self) = @_;

    my $db_list;
    my @du = ( "du", $self->volume_path, );

    # Get a list of files
    $self->log->info("Getting Directory Sizes");
    my $du_path = can_run('du') or $self->log->warn('du not installed!');
    my $du_command = [ 'du', $self->volume_path, '|', 'sort', '-n', '-k', '1' ];
    my $du_output = $self->run_command($du_command);

    #     my @dir_size = qx{@du};
    my $compara_db;
    my $species_filter = $self->species eq 'all' ? '.*' : $self->species;
    my $start_from_filter = $self->start_from ? $self->start_from : '.*';
    my $dbgroup_to_dbnames = {
        core => [qw(core coreexpressionatlas coreexpressionest funcgen otherfeatures rnaseq vega )],
        coreplusvariation =>
          [qw(core coreexpressionatlas coreexpressionest funcgen otherfeatures rnaseq vega variation)],
        corepluscompara => [qw(core coreexpressionatlas coreexpressionest funcgen otherfeatures rnaseq vega compara )],
        all => [qw(core coreexpressionatlas coreexpressionest funcgen otherfeatures rnaseq vega compara variation )],
        variation => [qw(variation)],
        compara   => [qw(compara)],

    };
    my $wanted_dbs = $dbgroup_to_dbnames->{ $self->dbtypes };

    my @lines;
    my ( $start, $end ) = ( $self->start_from, $self->end_at );
    map { push @lines, $_ if ( /$start/ .. /$end/ ) } split "\n", $du_output;
    $self->log->fatal( "Nothing found at " . $self->volume_path ) && die if @lines == 0;

    my $db_hash;

    foreach my $du_line (@lines) {
        my ( $db_size, $db_path ) = split /\s+/, $du_line;
        my ( $volume, $directories, $db_name ) = File::Spec->splitpath($db_path);
        my ( $db_species, $db_release, $db_type );
        if ( ( $db_species, $db_type, $db_release ) = $db_name =~ /^([a-z]+_[a-z]+)_([a-z]+)_(\d+)_\w+$/ ) {
            next unless $db_species =~ /$species_filter/;
            next unless $db_type ~~ @$wanted_dbs;
            push @{ $db_hash->{$db_species} }, [ $db_size, $db_path, $db_type, $db_release, $db_name ];
        }
        elsif ( $db_name =~ /ensembl_compara_(\d+)/ && 'compara' ~~ @$wanted_dbs ) {
            push @{ $db_hash->{$db_name} }, [ $db_size, $db_path, 'compara', $1, $db_name ];
        }
    }
    foreach my $species ( keys %$db_hash ) {

        my $dbs = $db_hash->{$species};
        my $volume = EnsCloud::Image::VolumeBundle::Volume->new( species => $species );
        foreach my $db (@$dbs) {
            my $db_detail = EnsCloud::Image::VolumeBundle::Volume::DatabaseDetails->new(
                myd_path => $db->[1],
                name     => $db->[4],
                type     => $db->[2],
                size     => $db->[0] / 1024 / 1024

            );
            $volume->add_db($db_detail);
            $volume->add_size( $db_detail->size );
            my $tag = "e$db->[3] MYD $species [" . ( join " ", map { $_->type } $volume->all_dbs ) . "]";

            #             my $tag = join " ", map { $_->name } $volume->all_dbs;
            $volume->tag($tag);
            $volume->ensembl_release( $db->[3] );
        }
        $volume->sort_in_place_curried;
        $self->volume_add_to_queue($volume);
    }
    return;
}

__PACKAGE__->meta->make_immutable;

1;
__END__



