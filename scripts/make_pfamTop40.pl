#!/usr/local/bin/perl 



use EnsWeb;
use Bio::EnsEMBL::DBLoader;
use Getopt::Long;

my $number = 100;
my $help;
my $db;

&GetOptions( 
	     'h|help'    => \$help,
	     );

if ($help) {
    exec('perldoc', $0);
}



eval {
     my $locator = &EnsWeb::get_locator();
#    my $locator = "Bio::EnsEMBL::DBSQL::Obj/host=obi-wan;port=3306;dbname=ensembl;user=ensro;pass=";
     $db =  Bio::EnsEMBL::DBLoader->new($locator);
};


if( $@ ) {
    print "<p>Warning! Exception<p>\n<pre>\n$@\n</pre>\n";
   exit(0);
}
     

my $sth     = $db->prepare("select sf.hid,tr.gene,count(*) from "     .
			   "transcript as tr,exon_transcript as et, " .
			   "supporting_feature as sf "                .
			   "where "                                   . 
			   "sf.exon = et.exon and "                   .
			   "et.transcript = tr.id and "               .
			   "sf.name = 'hmmpfam' "  
                    .
			   "group by sf.hid,tr.gene");

$sth->execute;
my $domain;

my @hits;

while (my $rowhash = $sth->fetchrow_hashref) {
    my $count  = $rowhash->{'count(*)'};
    my $dom    = $rowhash->{'hid'};

    $domain->{$dom}{genes}++;
    $domain->{$dom}{count} += $count;
}



hits2html($domain,$number);

1;

############################################################################################
sub hits2html {
    my ($domain,$number) = @_;

    my $numhits = scalar(keys %$domain);
    my $date    = `date`;
    chomp($date);

    my @domids = sort { $domain->{$b}{genes} <=> $domain->{$a}{genes}} keys %$domain;
    print("<h3>Top $number Pfam Domains in the Ensembl Database</h3>\n");
    print("This table was generated on $date<p>\n");

    print("<center><table BORDER=1 CELLPADDING=5 width=80%>\n");
    print("<tr class=\"violet3\"><Th>No.</Th><Th>Domain name</Th><th>Number of genes</th><th>Number of Ensembl hits</th><Th>No.</Th><Th>Domain name</Th><th>Number of genes</th><th>Number of Ensembl hits</th></tr>\n");

    for (my $i = 0; $i< $number/2; $i++){
	my $tmpdom1 = $domain->{$domids[$i]};
	my $tmpdom2 = $domain->{$domids[$i + $number/2]};

	my $name1  = $domids[$i];
	my $gene1  = $tmpdom1->{genes};
	my $count1 = $tmpdom1->{count};

	my $name2  = $domids[$i + $number/2];
	my $gene2  = $tmpdom2->{genes};
	my $count2 = $tmpdom2->{count};

        my $cnt1 = $i+1;
        my $cnt2 = $i+($number/2)+1;

	my $acc1 = pfamid2acc($name1);
	my $acc2 = pfamid2acc($name2);

	print("<tr align=\"center\"><td class=\"violet1\"><B>$cnt1</B></td>" .
	      "<td><A href=\"http://www.sanger.ac.uk/cgi-bin/Pfam/getacc?$acc1\">$name1</a></td>" .
	      "<td><A href=\"/perl/pfamview?pfamentry=$name1\">$gene1</a></td>".
	      "<td>$count1</td>");
	print("<td class=\"violet1\"><B>$cnt2</B></td><td>" .
	      "<td><A href=\"http://www.sanger.ac.uk/cgi-bin/Pfam/getacc?$acc2\">$gene2</a></td>" .
	      "<A href=\"/perl/pfamview?pfamentry=$name2\">$name2</a></td>" .
	      "<td>$count2</td></tr>\n");

    }

    print("</table></center><br><BR>\n");
    
}

sub pfamid2acc {
    my ($domain) = @_;

    my $pfamid = `/usr/local/pubseq/bin/getz -f acc '[pfam-id:$domain]'`;

    $pfamid =~ /\s*AC(\s+)(\S+)(.*)/;
    $pfamid = $2 if ($2);

    return $pfamid;
}
############################################################################################


