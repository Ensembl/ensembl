#!/usr/local/bin/perl 

# makes GFF stuff for a contig.

package Apache::EnsEMBL::contig2gff;

BEGIN {
    push(@INC,"/mysql/ensembl/www/apache_1.3.6/cgi-test/ensembl/modules");
    push(@INC,"/mysql/ensembl/www/apache_1.3.6/cgi-test/bioperl-live");
}

use Apache::Constants qw(:response :methods :http);
use Apache::Log ();
use CGI;
use Bio::EnsEMBL::DB::Obj;
use Apache::EnsEMBL::Header;
use Apache::EnsEMBL::Footer;
use DBI;
use strict;


sub handler {

    my $r = shift;

    $r->err_header_out('EnsEMBL-Error'=>"Problem in module Apache::EnsEMBL::contig2gff");
    $r->custom_response(SERVER_ERROR, "/Crash");

    if($ENV{PERL_SEND_HEADER}) {
        print "Content-type: text/html\n\n";
    }
    else {
        $r->content_type('text/html');
        $r->send_http_header;
    }
    if ($r->header_only) {
		return OK;
	}

	my ($header, $footer);
    &Apache::EnsEMBL::Header::make_ensembl_header(\$r, \$header);
    &Apache::EnsEMBL::Footer::make_ensembl_footer(\$r, \$footer);
	
	my %params = $r->method eq 'POST' ? $r->content : $r->args;
	my $contigid = $params{'contig'};

	my @features;
	eval {
	    my $db||=new Bio::EnsEMBL::DB::Obj( -user => 'root', -db => 'pog' , -host =>
		'caldy.sanger.ac.uk', -debug => 5);
	    #$r->warn("DB Handle: $db");
	    my $contig = $db->get_Contig($contigid);
	    @features = $contig->get_all_SeqFeatures;
	};

    $r->print($header);

	
	if( $@ ) {
	    print "<P>	<BLOCKQUOTE CLASS=\"error\">Warning! Exception
						<pre>\n$@\n</pre>
					</BLOCKQUOTE>
				</P>";
	} else {

	    print "<BLOCKQUOTE><pre>\n";
	    foreach my $sf ( @features ) {
		# $sf is Bio::SeqFeature::Generic object.
		print $sf->gff_string, "\n";
	    }
	    print "</pre></BLOCKQUOTE>\n";
	}
	
    #$r->warn("HTTP dump:\n", $r->as_string);
	
    $r->print($footer);
	return OK;


} # end of handler
  
1;

__END__
