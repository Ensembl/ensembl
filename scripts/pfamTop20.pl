#!/usr/local/bin/perl 

BEGIN {
    unshift(@INC,"../modules");
    unshift(@INC,"~/bioperl-live");
}

use Bio::EnsEMBL::DBLoader;
use Getopt::Long;

my $host   = 'localhost';
my $port   = '410000';
my $dbname = 'ensembl';
my $dbuser = 'ensro';
my $module = "Bio::EnsEMBL::DBSQL::Obj";
my $dbpass = undef;
my $number = 20;
my $help;

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

my $locator = "$module/host=$host;port=$port;dbname=$dbname;user=$dbuser;pass=$dbpass";
my $db      =  Bio::EnsEMBL::DBLoader->new($locator);
    

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

sub hits2html {
    my (@hits) = @_;

    my $numhits = scalar(@hits);
    my $date    = `date`;
    chomp($date);

    print("<h1>The Top $numhits Pfam Doains in the EnsEMBL Database</h1>\n");
    print("This table was generated on $date<p>\n");

    print("<center><table BORDER=2 CELLPADDING=5>\n");
    print("<tr><TD>Domain name</TD><td>Number of EnsEMBL hits</td></tr>\n");

    foreach my $hit (@hits) {
	my $name  = $hit->{domain};
	my $count = $hit->{count};

	print("<tr><td><A href=\"/cgi-bin/pfamview.pl?pfamentry=$name\">$name</a></td><td>$count</td></tr>\n");

    }

    print("</table></center>\n");
    
}


