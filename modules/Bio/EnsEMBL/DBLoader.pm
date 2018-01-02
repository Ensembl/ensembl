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


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::DBLoader - Run time database loader

=head1 SYNOPSIS

    $db =
      Bio::EnsEMBL::DBLoader->new( "Bio::EnsEMBL::DBSQL::DBAdaptor/"
        . "host=localhost;"
        . "dbname=homo_sapiens_core_19_34a;"
        . "user=ensro;" );

    # $db is a database object
    $db = Bio::EnsEMBL::DBLoader->standard();

    # equivalent to
    # Bio::EnsEMBL::DBLoader->new( $ENV{'ENSEMBL_DATABASE'} );

=head1 DESCRIPTION

This system provides a run-time loading of the database for ensembl,
allowing two things

    a) Only "using" the database module which is required for a
       particular implementation

    b) Providing a simple string method to indicate where the database
       is, allowing per sites defaults and other things as such


The string is parsed as follows:

Before the / is the Perl database object to load, after are the
parameters to pass to that database.  The parameters are series of
key=values separated by semi-colons.  These are passed as a hash to the
new method of the database object

=head1 METHODS

=cut

package Bio::EnsEMBL::DBLoader;

use strict;


=head2 new

  Arg [1]    : string $string
               An Ensembl database locator string.
  Example    : Bio::EnsEMBL::DBSQL::DBLoader->new("Bio::EnsEMBL::DBSQL::DBAdaptor/host=localhost;dbname=homo_sapiens_core_19_34a;user=ensro;"
  Description: Connects to an Ensembl database using the module specified in
               the locator string.
  Returntype : The module specified in the load string is returned.
  Exceptions : thrown if the specified module cannot be instantiated or the
               locator string cannot be parsed
  Caller     : ?
  Status     : Stable

=cut

sub new{
   my ($class,$string) = @_;
   my ($module,%hash);

   $string =~ /(\S+?)\/([\S+\s*]+)/ || die "Could not parse [$string] as a ensembl database locator. Needs database_module/params";
   $module = $1;
   my $param = $2;

   &_load_module($module);
   my @param = split(/;/,$param);
   foreach my $keyvalue ( @param ) {
       $keyvalue =~ /(\S+?)=([\S*\s*]*)/ || do { warn("In loading $keyvalue, could not split into keyvalue for loading $module. Ignoring"); next; };

       my $key = $1;
       my $value = $2;

       $hash{"-$key"} = $value;
   }
   
   my @kv = %hash;

   return "$module"->new(%hash);
}


sub _load_module{
  my ($modulein) = @_;
  my ($module,$load,$m);

  $module = "_<$modulein.pm";
  $load = "$modulein.pm";
  $load =~ s/::/\//g;
  
  return 1 if $main::{$module};
  eval {
    require $load;
  };
  if( $@ ) {
    print STDERR <<END;
$load: cannot be found
Exception $@

END
  ;
    return;
  }
  return 1;
}

1;
