#!/usr/local/bin/perl 

BEGIN {
    unshift(@INC,"../ensembl/modules");
    unshift(@INC,"../bioperl-live");
}

use EnsWeb;
require "ensembl-cgi-lib.pl";
use Getopt::Long;

my $number = 20;
my $help;
my $db;

&GetOptions( 
	     'dbtype:s'  => \$dbtype,
	     'host:s'    => \$host,
	     'port:n'    => \$port,
	     'dbname:s'  => \$dbname,
	     'dbuser:s'  => \$dbuser,
	     'dbpass=s'  => \$dbpass,
	     'module:s'  => \$module,
	     'number:n'  => \$number,
	     'h|help'    => \$help,
	     );

if ($help) {
    exec('perldoc', $0);
}



eval {
     my $locator = &EnsWeb::get_locator();
     $db =  Bio::EnsEMBL::DBLoader->new($locator);
};


if( $@ ) {
    print "<p>Warning! Exception<p>\n<pre>\n$@\n</pre>\n";
    exit(0);
}
     

my $sth     = $db->prepare("select count(*),hid from supporting_feature where name = 'hmmpfam' group by hid order by 1 desc limit $number");

$sth->execute;

my @hits;

while (my $rowhash = $sth->fetchrow_hashref) {
    my $tmp;
    my $count  = $rowhash->{'count(*)'};
    my $domain = $rowhash->{'hid'};

    $tmp->{count}  = $count;
    $tmp->{domain} = $domain;

    push(@hits,$tmp);
}

hits2html(@hits);

1;

############################################################################################
sub hits2html {
    my (@hits) = @_;

    my $numhits = scalar(@hits);
    my $date    = `date`;
    chomp($date);

    #print make_cgi_header();

    print("<h1>The Top $numhits Pfam Domains in the Ensembl Database</h1>\n");
    print("This table was generated on $date<p>\n");

    print("<center><table BORDER=2 CELLPADDING=5>\n");
    print("<tr><TD>Domain name</TD><td>Number of Ensembl hits</td></tr>\n");

    foreach my $hit (@hits) {
	my $name  = $hit->{domain};
	my $count = $hit->{count};

	print("<tr><td><A href=\"/perl/pfamview?pfamentry=$name\">$name</a></td><td>$count</td></tr>\n");

    }

    print("</table></center><br><BR>\n");
    
    #print make_cgi_footer();
}
############################################################################################


