#!/usr/bin/env perl
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


# Designed to work on data retrieved from 
# ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20110106_recombination_hotspots/
#
# Imports the recombination hotspots from 1000 Genomes into Ensembl

use strict;
use warnings;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::SimpleFeature;
use Bio::EnsEMBL::DBSQL::SimpleFeatureAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis;

use Getopt::Long;
use File::Fetch;
use IO::Dir;
use IO::File;

my ($url,$db_name,$db_host,$db_user,$db_pass,$db_port,$db_version,$help,$species);
$species = "human";

GetOptions ("url=s" => \$url,
            "db_name=s" => \$db_name,
            "db_host=s" => \$db_host,
            "db_user=s" => \$db_user,
            "db_pass=s" => \$db_pass,
            "db_port=s" => \$db_port,
            "db_version=s" => \$db_version,
            "species=s" => \$species,
            "h!"        => \$help,
            "help!"     => \$help,
);

if ($help) {&usage; exit;}
unless ($url and $db_name and $db_host) {print "Insufficient arguments\n"; &usage; exit;}


my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -species => $species,
    -group => 'core',
    -dbname => $db_name,
    -host => $db_host,
    -user => $db_user,
    -pass => $db_pass,
    -port => $db_port,
    -db_version => $db_version,
);

Bio::EnsEMBL::Registry->add_DBAdaptor($species,'core',$dba);

#Bio::EnsEMBL::Registry->load_registry_from_db(
#    -host => $db_host,
#    -user => $db_user,
#    -pass => $db_pass,
#    -port => $db_port,
#    -db_version => $db_version,
#);

my $file_fetch = File::Fetch->new(uri=>$url);
my $archive_path = $file_fetch->fetch() or die "Unable to get data from given URL. ".$file_fetch->error;

system('tar','-xzf',$archive_path);

my %directory;
tie %directory, 'IO::Dir',".";
foreach my $file (keys %directory) {
    if ($file =~ /\.txt$/) {
        my $fh = IO::File->new($file);
        &process_file($fh);
        $fh->close;
    }
}
print "Finished. Feel free to delete downloaded data in this directory.\n";

sub process_file {
    my $file_handle = shift;
    <$file_handle>; #strip header
    
    my $simple_feature_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species,'core', 'SimpleFeature');
    my $slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species,'core', 'Slice');
    
    my $analysis = new Bio::EnsEMBL::Analysis (
        -description => "The map was generated using the HapMap Phase II data and human genome assembly NCBI35 using LDhat as described in the 2007 HapMap paper (Nature, 18th Sept 2007).

The map was then converted from NCBI35 to GRCh37 coordinates and inspected for regions in which
the genome assembly had be rearranged.

See ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20110106_recombination_hotspots/",
        -display_label => 'HapMap Phase II genetic recombination map',
        -displayable => 1,
        -logic_name => 'human_1kg_hapmap_phase_2',
    );

    my $previous_chromosome;
    my $slice;
    my @features;
    while (my $line = <$file_handle>) {
        my ($chromosome,$position,$score) = split(/\t+/,$line);
        $chromosome =~ s/^chr//;  # chr prefix is for UCSC
        $chromosome =~ s/\_.+$//; # for removing PAR info
        if (!$slice || $previous_chromosome ne $chromosome) {
            $slice = $slice_adaptor->fetch_by_region('toplevel', $chromosome);
        }
        my $simple_feature = Bio::EnsEMBL::SimpleFeature->new(
            -start => $position,
            -end => $position,
            -score => $score,
            -analysis => $analysis,
            -slice => $slice,
            -strand => 1,
            -display_label => 'Recombination hotspot',
        );
        
        push @features, $simple_feature;
        $previous_chromosome = $chromosome;
    }
    $simple_feature_adaptor->store(@features);
}

sub usage {
    print "Launching instructions:
    Run from a folder you are happy to have filled with files.
    
Options:
    
    -url      Supply the URL to download from
    -db_name  The DB to add these features to
    -db_host  Hostname for the DB
    -db_user
    -db_pass
    -db_port
    -db_version
    -species
    
    -help
";    
}
