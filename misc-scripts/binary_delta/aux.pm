package aux;

# $Id$
#
# aux.pm
#
# Auxiliary subroutines for the apply.pl and build.pl programs.
#
# Author: Andreas Kahari, <andreas.kahari@ebi.ac.uk>
#

require Exporter;
@ISA = qw(Exporter);

@EXPORT_OK   = qw( make_human_readable make_checksum
		   usage_apply do_decompress
		   usage_build do_compress );

%EXPORT_TAGS = (
    default  => [ qw( make_human_readable make_checksum ) ],
    apply    => [ qw( usage_apply do_decompress ) ],
    build    => [ qw( usage_build do_compress ) ]
);

use Digest::MD5;
use Compress::Zlib;

#================= :default

# Compute the MD5 checksum of a file.  Returns the checksum as a
# hex string.
sub make_checksum
{
    my $file_path = shift;

    my $digest = new Digest::MD5;

    open FILE, $file_path or die $!;
    binmode FILE;

    $digest->addfile(*FILE);

    my $hex = $digest->hexdigest;
    close FILE;

    return $hex;
}

# Converts a byte count to a form more easly read by humans.
# Returns a string consisting of an float (two decimal places),
# a space, and a suffix.
sub make_human_readable
{
    my $bytes = shift;

    my @prefix = qw(b Kb Mb Gb Tb Pb);
    my $step = 0;

    while ($bytes > 10000) {
	$bytes /= 1024;
	++$step;
    }

    return sprintf("%.2f %s", $bytes, $prefix[$step]);
}

#================= :apply

# Display usage information for the apply.pl program.
sub usage_apply
{
    my $opts = shift;

    print STDERR <<EOT;
Usage:  $0 [options] [--] database old_v new_v

database    The database to work on, e.g. "homo_sapiens_core".
old_v       The older version, e.g. "11_31".
new_v       The newer version, e.g. "12_31".

The options may be any of these:

-c cmd  Path to xdelta executable.
        Default: "$opts->{'c'}".
-s path Path to the directory where the delta directory is stored.
        Default: "$opts->{'s'}"
-d path Path to the directory holding the old version of the
        database, and where the new version of the database
        should be created.
        Default: "$opts->{'d'}"

EOT
}

# Decompress a file.
sub do_decompress
{
    my $zfile_path = shift;
    my $file_path  = shift;

    open(OUT, '>' . $file_path) or die $!;
    binmode OUT;

    my $gz = gzopen($zfile_path, "r");

    if (!defined($gz)) {
	close OUT;
	die $gzerrno;
    }

    my $buffer;
    while ((my $bytesread = $gz->gzread($buffer)) != 0) {
	print OUT substr($buffer, 0, $bytesread);
    }

    $gz->gzclose();
    close OUT;
}

#================= :build

# Display usage information for the build.pl program.
sub usage_build
{
    my $opts = shift;

    print STDERR <<EOT;
Usage:  $0 [options] [--] database old_v new_v

database    The database to work on, e.g. "homo_sapiens_core".
old_v       The older version, e.g. "11_31".
new_v       The newer version, e.g. "12_31".

The options may be any of these:

-c cmd  Path to xdelta executable.
        Default: "$opts->{'c'}".
-s path Path to the directory where the databases are stored.
        Default: "$opts->{'s'}"
-d path Path to the directory within which the delta
        directory should be created.
        Default: "$opts->{'d'}"

EOT
}

# Compress a file.
sub do_compress
{
    my $file_path  = shift;
    my $zfile_path = shift;

    open(IN, $file_path) or die $!;
    binmode IN;

    my $gz = gzopen($zfile_path, "w");

    if (!defined($gz)) {
	close IN;
	die $gzerrno;
    }

    while (defined(my $input = <IN>)) {
	$gz->gzwrite($input) or die $gzerrno;
    }

    $gz->gzclose();
    close IN;
}

1;
