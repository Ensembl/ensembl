use strict;
use warnings;

use Getopt::Long;

my ($file, $user, $password, $verbose, $force, $help, $schema);

GetOptions ('file=s'      => \$file,
            'schema=s'    => \$schema,
            'user=s'      => \$user,
            'password=s'  => \$password,
            'verbose'     => \$verbose,
            'force'       => \$force,
            'help'        => sub { &show_help(); exit 1;} );

my $file = shift;

open FILE, $file or die( "Please provide config file for the transfer" );

while( $line = <> ) {
  next if $line =~ /^#/;
  next if !$line;
  
  my ( $species, $host, $source_db_name, $target_db_name ) = 
    split( "\t", $line );

  eval {
    require SeqStoreConverter::$species;

    push( @all_species_converters, 
          SeqStoreConverter::$species->new( $user, password, $host, 
                                            $source_db_name, $target_db_name, 
                                            $schema, $force, $verbose ));
  };
  if( $@ ) {
    die ( "Species $species produced error: \n$@" );
  }
}

for my $converter ( @all_species_converters ) {
  
  $converter->create_coord_systems();
  $converter->create_seq_region();
  $converter->create_assembly();

  $converter->transfer_genes();
  $converter->transfer_prediction_transcripts();
  $converter->transfer_features();

  $converter->copy_tables();
}

