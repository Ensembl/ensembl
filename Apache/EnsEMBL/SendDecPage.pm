package Apache::EnsEMBL::SendDecPage;
       
use strict;
use Apache::Constants qw(:response :methods :http);
use Apache::File ();
use Apache::Log ();
use Apache::EnsEMBL::Header;
use Apache::EnsEMBL::Footer;
use Apache::EnsEMBL::NavBar;
use Apache::EnsEMBL::NavTrail;

#############################################################
# Mod_perl request handler all /htdocs pages
#############################################################
sub handler {
    my $r = shift;

    $r->err_header_out('EnsEMBL-Error'=>"Problem in module Apache::EnsEMBL::SendDecPage");
    $r->custom_response(SERVER_ERROR, "/Crash");
    
    
    if ($r->content_type ne 'text/html') {
        return DECLINED;
    }
    
    if ((my $rc = $r->discard_request_body) != OK) {
        return $rc;
    }
    
    if ($r->method_number == M_INVALID) {
        $r->log->error("Invalid method in request ", $r->the_request);
        return NOT_IMPLEMENTED;
    }

    if ($r->method_number == M_OPTIONS) {
        return DECLINED;
    }

    if ($r->method_number == M_PUT) {
        return HTTP_METHOD_NOT_ALLOWED;
    }

    if (-d $r->finfo) {
        return DECLINED;
    }

    unless (-e $r->finfo) {
        $r->log->error("File does not exist: ", $r->filename);
        return NOT_FOUND;
    }

    if ($r->method_number != M_GET) {
        return HTTP_METHOD_NOT_ALLOWED;
    }

    my $fh = Apache::File->new($r->filename);
    unless ($fh) {
        $r->log->error("File permissions deny server access: ", $r->filename);
        return FORBIDDEN;
    }

    $r->update_mtime(-s $r->finfo);
    $r->set_last_modified;
    $r->set_etag;
    if((my $rc = $r->meets_conditions) != OK) {
        return $rc;
    }

    #$r->set_content_length;

    if($ENV{PERL_SEND_HEADER}) {
        print "Content-type: text/html\n\n";
    }
    else {
        $r->content_type('text/html');
        $r->send_http_header;
    }
    
    my ($trail, $nav_start, $nav_end, $nav_menu);
    my ($title, $gif);  # currently unused
    my ($header, $footer);
    
    &Apache::EnsEMBL::NavTrail::make_nav_trail(\$r, \$trail);
    
    &Apache::EnsEMBL::NavBar::make_nav_menu (\$r, \$nav_menu);
    &Apache::EnsEMBL::NavBar::make_nav_start(\$r, \$nav_start, \$nav_menu);
    &Apache::EnsEMBL::NavBar::make_nav_end  (\$r, \$nav_end);

    &Apache::EnsEMBL::Header::make_ensembl_header(\$r, \$header);
    &Apache::EnsEMBL::Footer::make_ensembl_footer(\$r, \$footer);
    
    unless ($r->header_only) {
        $r->print($header);
        $r->print($trail);

        $r->print($nav_start);
        while (<$fh>){
            $r->print($_);
        }
        $r->print($nav_end);

        $r->print($footer);
    }

    
    close $fh;

    my $http_header = $r->headers_in->{'X-Forwarded-For'};
        if( my $ip = (split /,\s*/, $http_header)[-1] ) {
            $r->connection->remote_ip($ip);
        }
        
    #$r->warn("HTTP dump:\n", $r->as_string);
    
    
    return OK;
} # end of handler


#############################################################

1;

__END__

#
# EnsEMBL module for Apache::EnsEMBL::SendDecPage
#
# Cared for by Tony Cox <avc@sanger.ac.uk>
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs after the code

=head1 NAME

Apache::EnsEMBL::SendDecPage - Apache Mod_perl module to serve "decorated" HTML files

=head1 SYNOPSIS

=head2 General

This mod_perl module is used to take a skeleton HTML file and wrap it in a header,
footer and side navigation bar.

=head1 DESCRIPTION (some features not yet implemented)

This module takes standard HTTP requests received by Apache and locates the
corresponding HTML file which is then returned wrapped suitably. It is fully 
HTTP/1.1 compliant. This module also includes an optimization which checks 
whether the browser already has the page in its disk cache and is only
requesting the page if it has changed. 

The header displays a "bread crumb" trail to indicate where the current document 
lies in the server tree structure and a reference to a CSS which is 
used to apply a style. URLs can optionally be mapped to a "friendly" name by placing
and entry in the "/urlmappings.txt" file in the server "/htdocs" directory.

The navigation side bar menu is created dynamically by incorporating the contents of a 
"nav.conf" file in the same directory as the HTML file (if it exists). This contains
a list of navigation options that will be written to the side bar. 
If the config file does not exist the default navigation bar is used (/def_nav.conf).

The footer shows a last modification date of the HTML file and a contact email address.


=head1 RELATED MODULES

See also: Apache::EnsEMBL::SendPage, Apache::EnsEMBL::Header, Apache::EnsEMBL::Footer

=head1 FEED_BACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
EnsEMBL modules. Send your comments and suggestions to one of the
EnsEMBL mailing lists.  Your participation is much appreciated.

  http://ensembl.ebi.ac.uk/?     - About the mailing lists

=head2 Reporting Bugs

Report bugs to the EnsEMBL bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via
email or the web:

  ?@ensembl.ebi.ac.uk
  http://ensembl.ebi.ac.uk/?

=head1 AUTHOR - Tony Cox

Email - avc@sanger.ac.uk


=cut
