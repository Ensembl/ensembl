package Apache::EnsEMBL::ServerError;
       
use strict;
use Apache::Constants qw(:response :methods :http);
use Apache::File ();
use Apache::Log ();
use CGI qw(:html);

sub handler {

	my $r = shift;
	my $error = $r->err_header_out('EnsEMBL-Error');
	unless ($error){
		 $error = 'unknown (no specific information available)';
	}
	
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
		h1('EnsEMBL Server Error'),
		p("Sorry, an error occured while the server was processing your request.
		Please email a report , quoting any additional information given below,
		along with the URL, to the server administrator using the link below."),

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

		p(b("Error (if known):")),
		p(
			blockquote({-class=> "error"},
				strong($ENV{'REDIRECT_ERROR_NOTES'}),
			)
		),
			
		end_html
		);
		
        my $footer = "";
        &Apache::EnsEMBL::Footer::make_ensembl_footer(\$r, \$footer);
        $r->print($footer);

	return OK;
}

1;

__END__

		
	
	
	
	
	);


} # end of hndler
