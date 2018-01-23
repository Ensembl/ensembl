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

Bio::EnsEMBL::Utils::ScriptUtils;

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 METHODS

=cut

package Bio::EnsEMBL::Utils::ScriptUtils;

use strict;
use warnings;
no warnings 'uninitialized';

use Exporter;
our @ISA = qw(Exporter);

our @EXPORT_OK = qw(
  user_proceed
  commify
  sort_chromosomes
  parse_bytes
  directory_hash
  path_append
  dynamic_use
  inject
);


=head2 user_proceed

  Arg[1]      : (optional) String $text - notification text to present to user
  Example     : # run a code snipped conditionally
                if ($support->user_proceed("Run the next code snipped?")) {
                    # run some code
                }

                # exit if requested by user
                exit unless ($support->user_proceed("Want to continue?"));
  Description : If running interactively, the user is asked if he wants to
                perform a script action. If he doesn't, this section is skipped
                and the script proceeds with the code. When running
                non-interactively, the section is run by default.
  Return type : TRUE to proceed, FALSE to skip.
  Exceptions  : none
  Caller      : general

=cut

sub user_proceed {
  my ($text, $interactive, $default) = @_;

  unless (defined($default)) {
    die("Need a default answer for non-interactive runs.");
  }

  my $input;

  if ($interactive) {
    print "$text\n" if $text;
    print "[y/N] ";
  
    $input = lc(<>);
    chomp $input;
  } else {
    $input = $default;
  }
  
  if ($input eq 'y') {
    return(1);
  } else {
    print "Skipping.\n" if ($interactive);
    return(0);
  }
}


=head2 sort_chromosomes

  Arg[1]      : (optional) Hashref $chr_hashref - Hashref with chr_name as keys
  Example     : my $chr = { '6-COX' => 1, '1' => 1, 'X' => 1 };
                my @sorted = $support->sort_chromosomes($chr);
  Description : Sorts chromosomes in an intuitive way (numerically, then
                alphabetically). If no chromosome hashref is passed, it's
                retrieve by calling $self->get_chrlength()
  Return type : List - sorted chromosome names
  Exceptions  : thrown if no hashref is provided
  Caller      : general

=cut

sub sort_chromosomes {
    my @chromosomes = @_;
    
    return (sort _by_chr_num @chromosomes);
}


=head2 _by_chr_num

  Example     : my @sorted = sort _by_chr_num qw(X, 6-COX, 14, 7);
  Description : Subroutine to use in sort for sorting chromosomes. Sorts
                numerically, then alphabetically
  Return type : values to be used by sort
  Exceptions  : none
  Caller      : internal ($self->sort_chromosomes)

=cut

sub _by_chr_num {
    my @awords = split /-/, $a;
    my @bwords = split /-/, $b;

    my $anum = $awords[0];
    my $bnum = $bwords[0];

    if ($anum !~ /^[0-9]*$/) {
        if ($bnum !~ /^[0-9]*$/) {
            return $anum cmp $bnum;
        } else {
            return 1;
        }
    }
    if ($bnum !~ /^[0-9]*$/) {
        return -1;
    }

    if ($anum <=> $bnum) {
        return $anum <=> $bnum;
    } else {
        if ($#awords == 0) {
            return -1;
        } elsif ($#bwords == 0) {
            return 1;
        } else {
            return $awords[1] cmp $bwords[1];
        }
    }
}


=head2 commify

  Arg[1]      : Int $num - a number to commify
  Example     : print "An easy to read number: ".$self->commify(100000000);
                # will print 100,000,000
  Description : put commas into a number to make it easier to read
  Return type : a string representing the commified number
  Exceptions  : none
  Caller      : general
  Status      : stable

=cut

sub commify {
  my $num = shift;

  $num = reverse($num);
  $num =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;

  return scalar reverse $num;
}


sub parse_bytes {
  my $bytes = shift;

  my @suffixes = qw(bytes kb Mb Gb Tb);

  my $length = length($bytes);
  my $order = int(($length-1)/3);

  my $parsed = sprintf('%.1f', $bytes/10**(3*$order));

  return "$parsed ".$suffixes[$order];
}


sub directory_hash {
  my $filename = shift;

  my (@md5) = md5_hex($filename) =~ /\G(..)/g;
  return join('/', @md5[0..2]);
}


sub path_append {
  my $path1 = shift;
  my $path2 = shift;

  # default to current directory
  $path1 = '.' unless (defined($path1));

  my $return_path = "$path1/$path2";

  unless (-d $return_path) {
    system("mkdir -p $return_path") == 0 or
      die("Unable to create directory $return_path: $!\n");
  }
  
  return $return_path;
}


=head2 inject

  Arg [1]    : String $classname - The name of the class to require/import
  Example    : $self->inject('Bio::EnsEMBL::DBSQL::DBAdaptor');
  Description: Requires and imports the methods for the classname provided,
               checks the symbol table so that it doesnot re-require modules
               that have already been required.
  Returntype : true on success
  Exceptions : Warns to standard error if module fails to compile
  Caller     : internal

=cut

sub inject {
  my $classname = shift;
  my ($parent_namespace, $module) = $classname =~/^(.*::)(.*)$/ ?
                                      ($1,$2) : ('::', $classname);
  no strict 'refs'; ## no critic

  # return if module has already been imported
  return 1 if $parent_namespace->{$module.'::'};
  
  eval "require $classname"; ## no critic
  die("Failed to require $classname: $@") if ($@);

  $classname->import();
  
  return 1;
}


sub dynamic_use {
  return inject(@_);
}

1;

