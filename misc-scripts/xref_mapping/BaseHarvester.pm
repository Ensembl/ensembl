# Superclass providing common functionality for harvesters
# as well as the ability to run them all

package BaseHarvester;

use strict;

# --------------------------------------------------------------------------------

if (!defined(caller())) {

  find_and_run_harvesters();

}

# --------------------------------------------------------------------------------
# Find subclasses and run them

sub find_and_run_harvesters {

  foreach my $file (glob("*.pm")) {

    if (eval { require $file }) {

      $file =~ s/\..*$//;
      my $obj = $file->new();
      if ($obj->isa("BaseHarvester") && $file !~ /BaseHarvester/) {
	print "Running " . $file . "\n";
	$obj->run();
      }

    } else {

      warn("Error during require of " . $file . "\n");

    }
  }

}

# --------------------------------------------------------------------------------

sub taxonomy_ids {

  my $self = shift;

  return (4530, 55529);

}

# --------------------------------------------------------------------------------

sub new {

  my $self = {};
  bless $self, "BaseHarvester";
  return $self;

}

# --------------------------------------------------------------------------------

1;
