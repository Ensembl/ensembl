# EnsEMBL Gene reading writing adaptor for mySQL
#
# Copyright EMBL-EBI 2001
#
# Author: Graham McVicker
# based partially on
# CPAN module ex::implements by PDCAWLEY 
# 
# 
# Date : 06.08.2002
#

=head1 NAME

implements - a Pragma for perl 5.6.0 that allows Perl programmers to 
define and implement interfaces.  Throws compile-time errors when 
interface is not fully implemented. 

=head1 SYNOPSIS

use implements qw(Interface1 Interface2);

=head1 CONTACT

  Arne Stabenau: stabenau@ebi.ac.uk
  Graham McVicker: mcvicker@ebi.ac.uk
  Ewan Birney  : birney@ebi.ac.uk

=head1 APPENDIX

=cut


package implements;

use strict;

use Carp;

my %__INTERFACES;

sub import {
  my ($class, @interfaces) = @_;

  my $caller = caller(0);

  #print "storing interface defs for $caller\n";

  #put interfaces this caller is implementing into a hash, keyed by caller
  $__INTERFACES{$caller} = @interfaces;

  no strict 'refs';
  #update the inheritence of this caller to include its interfaces
  push @{"$caller\::ISA"}, @interfaces;
  use strict 'refs';
}



CHECK {
  my @errors = ();

  #validate each module that imported 'implements'
  foreach my $implementor (keys %__INTERFACES) {

    #print "Validating $implementor\n";

    my %seen_packages = ();
    my @checked_packages = ();
    
    no strict 'refs';
    my @unchecked_packages = @{"$implementor\::ISA"};
    use strict 'refs';

    #traverse the isa heirarchy of each implementor to construct 
    #a complete list (including super class interfaces) of interfaces 
    #this caller must implement
    while(@unchecked_packages) {
      my $package = pop @unchecked_packages;
      
      #skip this package if we've seen it before
      next if($seen_packages{$package});
    
      #we haven't seen this package before so add it to the list and obtain 
      #its superclass interfaces
      eval "require $package";
      if($@) {
	croak "$implementor implements unknown interface $package";
      }

      push @checked_packages, $package;
      no strict 'refs';
      push @unchecked_packages, @{"$package\::ISA"};
      use strict 'refs';
    }

    #check the implementation of all the superclass interfaces
    foreach my $interface (@checked_packages) {
      
      #obtain the interface's methods from it's symbol table
      no strict 'refs';
      my %symbol_table = %{"$interface\::"};
      use strict 'refs';
      
      while(my ($method, $value) = each (%symbol_table)) {
	local (*alias);
	#use type globbing to determine discard non-method symbols
	*alias = $value;
	next unless defined &alias;
	
	#print "implements: checking $implementor for implementation of "
	#  . "$interface\::$method()\n";

	#check the implentor's symbol table to see if this method exists
	no strict 'refs';
	my %implementation_symbols = %{"$implementor\::"};
	use strict 'refs';

	local *a;
	*a = $implementation_symbols{$method};
	unless(defined &a) {
	  push @errors, "Error: $implementor implements interface " . 
                        "'$interface' but does not define method '$method'\n";
	}
	
	#print "implementor contains symbol $method\n";

      }
      
    }
  }
  
  if(@errors) {
    croak @errors;
  }

  return 1;
}

1;
  
