
#
# BioPerl module for Bio::EnsEMBL::DBLoader
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBLoader - Run time database loader

=head1 SYNOPSIS

    $db = Bio::EnsEMBL::DBLoader->new("Bio::EnsEMBL::DBSQL::Obj/host=localhost;dbname=ensdev;user=ensro;");

    # $db is a database object
    $db = Bio::EnsEMBL::DBLoader->standard();
    # equivalent to Bio::EnsEMBL::DBLoader->new($ENV{'ENSEMBL_DATABASE'});


=head1 DESCRIPTION

This system provides a run-time loading of the database for ensembl, allowing two things

    a) Only "using" the database module which is required for a particular implementation

    b) Providing a simple string method to indicate where the database is, allowing per sites
defaults and other things as such


The string is parsed as follows:

Before the / is the Perl database object to load, after are the parameters to pass
to that database. The parameters are series of key=values separated by semi-colons.
These are passed as a hash to the new method of the database object

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DBLoader;

use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DBOLD::Obj;
use Bio::EnsEMBL::TimDB::Obj;
use strict;

=head2 new

 Title   : new
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub new{
   my ($class,$string) = @_;
   my ($module,%hash);

   $string =~ /(\S+?)\/(\S+)/ || die "Could not parse [$string] as a ensembl database locator. Needs database_module/params";
   $module = $1;
   my $param = $2;

   &_load_module($module);
   my @param = split(/;/,$param);
   foreach my $keyvalue ( @param ) {
       $keyvalue =~ /(\S+?)=(\S*)/ || do { warn("In loading $keyvalue, could not split into keyvalue for loading $module. Ignoring"); next; };

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
