package Apache::EnsEMBL::NotFound;
       
use strict;
use Apache::Constants qw(:response :methods :http);
use Apache::File ();
use Apache::Log ();
use CGI qw(:html);
use Apache::EnsEMBL::Header;
use Apache::EnsEMBL::Footer;

sub handler {

	my $r = shift;
	my $error = $ENV{'REDIRECT_ERROR_NOTES'};
	unless ($error){
		 $error = 'unknown (no specific information available)';
	}
	$error =~ s/\/mysql\/ensembl\/www\/apache_1.3.6\/htdocs//i;	

	$r->content_type('text/html');
	$r->send_http_header;
	return OK if $r->header_only;
	
	my $original_request = $r->prev;
	my $original_uri = $original_request ? $original_request->uri : '';
	my $admin = $r->server->server_admin;
	
	my $header = "";
        &Apache::EnsEMBL::Header::make_ensembl_header(\$r, \$header);
	
	$r->print($header);
	
	$r->print (
	
		start_html(	-title   => 'Server Error', 
				-bgcolor => 'white',
				-style   => {-src => '/EnsEMBL.css'}
			   ),
		h1('EnsEMBL Server: Page Not Found'),
		p("Sorry, the page you requested was not found on this server.
		Please check that you have typed in the correct URL or else use the site search
		facility to try and locate information you require.
		If you think an error has occurred please send email to the server administrator 
		using the link below."),

		p(b("The error was:")),
		p(
			blockquote({-class=> "error"},
				strong($error),
			)
		),

		p(b("URL:")),
		p(
			blockquote({-class=> "error"},
				strong($ENV{'REDIRECT_URL'}),
			)
		),

		p(b("HTTP Status Code:")),
		p(
			blockquote({-class=> "error"},
				strong($ENV{'REDIRECT_STATUS'}), 
			)
		),

		p(b("Request Method:")),
		p(
			blockquote({-class=> "error"},
				strong($ENV{'REDIRECT_REQUEST_METHOD'}), 
			)
		),
		p(b("Query String (if known):")),
		p(
			blockquote({-class=> "error"},
				strong($ENV{'REDIRECT_QUERY_STRING'}),
			)
		),
		
                )	
                ;
		
	my $footer = "";
        &Apache::EnsEMBL::Footer::make_ensembl_footer(\$r, \$footer);
	$r->print($footer);

	return OK;
}

1;

__END__

		
	
	
	
	
	);


} # end of hndler
