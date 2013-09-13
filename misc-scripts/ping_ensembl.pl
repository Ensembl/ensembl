#!/usr/bin/env perl

=head1 LICENSE

   Copyright (c) 1999-2013 The European Bioinformatics Institute and
   Genome Research Limited, and others.  All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

   1. Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer
      in the documentation and/or other materials provided with the
      distribution.

   3. The name "Ensembl" must not be used to endorse or promote products
      derived from this software without prior written permission.  For
      written permission, please contact helpdesk@ensembl.org

   4. Products derived from this software may not be called "Ensembl"
      nor may "Ensembl" appear in their names without prior written
      permission of the Ensembl developers.

   5. Redistributions in any form whatsoever must retain the following
      acknowledgement:

        "This product includes software developed by Ensembl
        (http://www.ensembl.org/)"

   THIS SOFTWARE IS PROVIDED BY THE ENSEMBL GROUP ``AS IS'' AND ANY
   EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
   PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE ENSEMBL GROUP OR ITS
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
   OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=head1 CONTACT

  Please email comments or questions to the Ensembl help desk
  <helpdesk@ensembl.org>.

=head1 NAME

  ping_ensembl.pl

=head1 SYNOPSIS

  # print usage 
  $ ping_ensembl.pl -h 

  # ping Ensembl with default species (Human)
  $ ping_ensembl.pl

  # ping Ensembl with user provided species
  $ ping_ensembl.pl -s "dog"

  # ping Ensembl with a different version (Human)
  $ ping_ensembl.pl -db_version 70

  # ping the US Ensembl mirror
  $ ping_ensembl.pl -ue

  # ping Ensembl Genomes with default species (arabidopsis thaliana)
  $ ping_ensembl.pl -eg

  # ping Ensembl Genomes with user provided species 
  $ ping_ensembl.pl -eg -s "oryza sativa japonica"

=head1 DESCRIPTION

  This script is used to detect if you can contact the Ensembl database 
  server with your current setup. The program will attempt to print out
  helpful hints about how to resolve your problems. If they still persist
  then please contact helpdesk@ensembl.org.

=head1 SUBROUTINES

=cut

use strict;
use warnings;

use File::Temp qw/tempfile/;
use Net::FTP;
use Getopt::Long;

#
# Default option values
#
my $help = 0;
my $host = 'ensembldb.ensembl.org';
my $user = 'anonymous';
my $port = 5306;
my $db_version = -1;

my $useast = 0;
my $ensembl_genomes = 0;
my $species = undef;
my $api_version = -1;

#
# Parse command-line arguments
#
my $options_ok = 
  GetOptions(
    "ue"            => \$useast,
    "eg"            => \$ensembl_genomes,
    "species=s"     => \$species,
    "db_version=i"  => \$db_version,
    "help"          => \$help);
($help or !$options_ok) && usage();

$useast and $ensembl_genomes and
 die "Cannot test Ensembl Genomes on the US mirror.\n" .
  "Options \"ue\" and \"eg\" are mutually exclusive\n";

$useast and $host = "useastdb.ensembl.org";

if ($ensembl_genomes) {
  $host = "mysql.ebi.ac.uk";
  $port = 4157;
  $species = "arabidopsis thaliana"
    unless defined $species;
}

eval {
  require DBI;
  require DBD::mysql;
  require Bio::Perl;
  require Bio::EnsEMBL::Registry;
  require Bio::EnsEMBL::ApiVersion;
  require Bio::EnsEMBL::LookUp if $ensembl_genomes;
  $api_version = Bio::EnsEMBL::ApiVersion::software_version();
  $db_version = $api_version if $db_version == -1; #if it was still -1 then it wasn't set. Default is current API version
  Bio::EnsEMBL::Registry->load_registry_from_db(
    -host       => $host,
    -port       => $port,
    -user       => $user,
    -db_version => $db_version,
  );
  $species = "human" unless defined $species;
  my $species_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor("$species", 'core');  
  print "Installation is good. Connection to Ensembl works and you can query the $species core database\n";
};
my $error = $@;

# If no error found then see if we've got all of our external modules available
if(!$error) {
  $error = '';
  eval {
    foreach my $module (qw/Compara Variation Funcgen/) {
      my $full_module = "Bio::EnsEMBL::${module}::DBSQL::DBAdaptor";
      eval "require $full_module;";
      if($@) {
        $error .= "\tMissing the checkout $module\n"; 
      }
    }
  };
}

# Check the current release of datafiles from the FTP site
my $ftp_version = -1;
eval {
  my $ftp = Net::FTP->new('ftp.ensembl.org', Debug => 0);
  $ftp->login("anonymous",'-anonymous@');
  $ftp->cwd('/pub');
  my ($fh, $filename) = tempfile();
  close($fh);
  $ftp->get('current_README', $filename);
  $ftp->quit();
  open($fh, '<', $filename);
  local $/ = undef;
  my $ftp_readme = <$fh>;
  close($fh);
  ($ftp_version) = $ftp_readme =~ /Ensembl Release (\d+) Databases/;
};

# Print all the errors which could have occured 
if($error) {
  print "ERROR: Error detected when connecting to Ensembl!\n";
  if($error =~ /DBI/) {
    print "\tCannot find the DBI perl module. Please install this using your package management system, cpan or cpanm. Please consult http://www.ensembl.org/info/docs/api/api_installation.html\n";
  }
  if($error =~ /mysql/) {
    print "\tCannot find the DBD::mysql perl module. Please install this using your package management system, cpan or cpanm. Also install the mysql libs. Please consult http://www.ensembl.org/info/docs/api/api_installation.html\n";
  }
  if($error =~ /Can't locate Bio\/E/) {
    print "\tLooks like you need to setup your PERL5LIB with the Ensembl API. Please consult http://www.ensembl.org/info/docs/api/api_installation.html\n";
  }
  if($error =~ /Can't locate Bio\/Perl/) {
    print "\tLooks like you need to setup your PERL5LIB with BioPerl. Please consult http://www.ensembl.org/info/docs/api/api_installation.html\n";
  }
  if($error =~ /Cannot connect to/) {
    print "\tCannot seem to contact EnsemblDB at '$host' with the username '$user'. Try running 'ping $host' or asking your systems about firewalls against port $port\n";
  }
  if($error =~ /internal name/) {
    print "\tSpecies was not found. You may have accidentally download the HEAD API version (found API release $db_version & public FTP release is $ftp_version). Please consult http://www.ensembl.org/info/docs/api/api_installation.html\n";
  }
  if($error =~ /Species not defined/) {
    print "\tSpecies was not found. You may have accidentally download the HEAD API version (told to load release $db_version, API version is $api_version & public FTP release is $ftp_version). Please consult http://www.ensembl.org/info/docs/api/api_installation.html\n";
  }
  if($error =~ /Missing the checkout/) {
    print "\tYour core installation was good but supplementary modules cannot be found. If you wish to access these other Ensembl resources add the libraries to your PERL5LIB:\n";
    print $error;
    #bail early
    exit 0;
  }
  print '='x80, "\n";
  print "If the problem persists please send the following error message to helpdesk\@ensembl.org\n";
  print $error;
  print '='x80, "\n";
  exit 1;
}


=head2 usage

  Arg []         : None
  Returntype     : None
  Example        :
  Description    : Print script usage string
  Exceptions     : None
  Caller         : General

=cut

sub usage {
  my $prog = `basename $0`; chomp($prog);
    
  print "Usage: $prog [OPTIONS]\n\n";
  print "Options:\n";
  print "  -ue                    Ping Ensembl US mirror\n";
  print "  -eg                    Ping Ensembl Genomes (can't be used together with \"ue\")\n";
  print "  -species <species>     Use species <species> (use double quotes if species name contains spaces)\n";
  print "  -db_version <version>  Use the specified version of Ensembl not the API version\n";
  print "  -help                  Print this message\n";
  print "\n\n";

  exit 1;
}
