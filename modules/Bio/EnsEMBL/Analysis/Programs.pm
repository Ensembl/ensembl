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

package Bio::EnsEMBL::Analysis::Programs;
use strict;
use vars qw( %Program_Paths );
use Carp;
use Cwd;
use Sys::Hostname;
use Bio::EnsEMBL::Utils::Exception qw ( throw ) ; 

sub import {
    my $pkg = shift;
    foreach (@_) {
        #print STDERR "importing: $_\n";
        $Program_Paths{ $_ } = 0;
    }
    my( $home, @PATH, @missing );
    
    $home = cwd() or die "Can't save cwd";
    my $H = [ $home, 1 ];
    
    @PATH = split /:/, $ENV{'PATH'};
    foreach (@PATH) {
        s|/?$|/|; # Append / to each path
    }
    
    # For each program, check there is an executable
    foreach my $program (keys %Program_Paths) {
        
        # Deal with paths
        if ($program =~ m|/|) {
            _go_home( $H );
            my $path = $program;
            # Deal with tildes
            $path =~ s{^~([^/]*)}{ $1 ? (getpwnam($1))[7]
                                      : (getpwuid($>))[7] }e;
            if (my $real = _is_prog( $H, $path )) {
                $Program_Paths{ $program } = $real;
            }            
        }
        # Or search through all paths
        else {
            foreach my $path (@PATH) {
                _go_home( $H );
                if (my $real = _is_prog( $H, $path, $program )) {
                    $Program_Paths{ $program } = $real;
                    last;
                }
            }
        }
    }
    _go_home( $H ); # Return to home directory
    
    # Make a list of all missing programs
    foreach my $program (keys %Program_Paths) {
        push( @missing, $program ) unless $Program_Paths{ $program };
    }
    
    # Give informative death message if programs weren't found
    if (@missing) { 
        throw("Unable to locate the following programs as '". (getpwuid($<))[0]. "' on host '". hostname(). "' :\t".
        join ( " --> " , @missing )) ; 
    }
}

# Recursive function which follows links, or tests final destination
sub _is_prog {
    my( $h, $path, $prog ) = @_;

    # Need to split path if $prog not provided
    unless ($prog) {
        ($path, $prog) = $path =~ m|(.*?)([^/]+)$|;
    }
    
    if (-l "$path$prog") {
        # Follow link
        _follow( $h, $path ) or return;
        unless (-x readlink($prog)) {
          confess "Can't read link '$path$prog' : $!";
	}
	my $link = $prog;
	$path = cwd() or confess "Can't determine cwd";
        return "$path/$prog";
    } elsif (-f _ and -x _) {
        # Return full path
        _follow( $h, $path ) or return;
        $path = cwd() or confess "Can't determine cwd";
        return "$path/$prog";
    } else {
        # Not a link or an executable plain file
        return;
    }
}

# To avoid unnecessary chdir'ing
sub _follow {
    my( $H, $path ) = @_;
    
    # Chdir without arguments goes to home dir.
    # Can't use defined in test since $path may contain
    # a real null string.
    if ( ! $path and $path ne '0' ) {
        return 1;
    } elsif (chdir($path)) {
        $H->[1] = 0;
        return 1;
    } else {
        return;
    }
}
sub _go_home {
    my( $H ) = @_;
    
    # Go home unless we're already there
    if ($H->[1] == 0) {
        if (chdir( $H->[0] )) {
            $H->[1] = 1;
        } else {
            confess "Can't go home to [ ", $H->[0], ' ]';
        }
    }
}

1;

__END__

=head1 NAME Programs

=head1 SYSNOPSIS

    use Bio::EnsEMBL::Analysis::Programs qw( efetch getz est2genome 
                     /usr/local/bin/this_one
                     ~me/some/path/my_prog
                     ~/../jane/bin/her_prog );

    # Can also do at run time
    Bio::EnsEMBL::Analysis::Programs->import( $someProg );

    $path_to_prog = $Bio::EnsEMBL::Analysis::Programs::Program_Paths{ $prog };

=head1 DESCRIPTION

B<Programs> is used to check at compile time for the
presence of executables which will be called from your
script.  Arguments passed via the use statement
can be just the program name, or an absolute or
relative path to the program.  Tildes are expanded
correctly (I<not> using "glob").  Failure to find any
one program is fatal, and a list of all failures is
printed, along with the host''s name.

If you want to check for a program during run time,
the import funtion can be called directly, as shown above.

The paths to each program found are stored in the
B<%Program_Paths> hash, which is keyed on the original
arguments passed.

=head1 BUGS

If the executable is in the root directory, then it''s found
path will appear as "//prog" in %Program_Paths, not "/prog".

=head1 AUTHOR

B<James Gilbert> Email jgrg@sanger.ac.uk











