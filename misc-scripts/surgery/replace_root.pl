#!/usr/bin/perl

use warnings;
use strict;

foreach my $file (@ARGV) {

  my $tmpFile = $file. "_tmp";
  open (IN, $file);
  open (OUT, ">$tmpFile");

  while (<IN>) {

    next if ($_ =~ /\@ISA\s+=\s+(.*Bio::EnsEMBL::Root.*)/); # remove @ISA line completely

    # Import
    if ($_ =~ /use\s+Bio::EnsEMBL::Root/) {
      print OUT $_;    # note inheritance kept in for now. Remove this line in due course
      print OUT "use Bio::EnsEMBL::Utils::Exception qw(throw warning);\n";
      print OUT "use Bio::EnsEMBL::Utils::Argument  qw(rearrange);\n";
    }

    # Methods - throw, warn, _rearrange
    elsif ($_ =~ /(\$[a-zA-Z]+)->(throw|warn|_rearrange)\((.*)/   ) {

      my $obj = $1;
      $obj =~ s/\$self\-\>//g; # call static methods instead of via $self
      my $met = $2;
      my $arg = $3;

      $met =~ s/warn/warning/;           # avoid clash with Perl builtin
      $met =~ s/_rearrange/rearrange/;   # no need for this to be "private"

      #print "Object: $obj  Method: $met   Arg: $arg\n";
      my $replacement = $obj . "->" . "$met(" . $arg ;
      #print "Changing $_ to $replacement \n";
      print OUT $replacement;

    }

    else {

      print OUT $_;

    }

  }

  close(IN);
  close(OUT);

  rename $tmpFile, $file;
}

# use something like
# find . -type f -name "*.pm" -o -name "*.t" -o -name "*.pl" | xargs perl /scratch/replaceRoot.pl
