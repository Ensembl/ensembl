package Apache::EnsEMBL::SendDecPage;
       
use strict;
use Apache::Constants qw(:response :methods :http);
use Apache::File ();
use Apache::Log ();
use Apache::EnsEMBL::Header;
use Apache::EnsEMBL::Footer;

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
    
    &make_nav_trail(\$r, \$trail);
    &make_nav_menu (\$r, \$nav_menu);
    &make_nav_start(\$r, \$nav_start, \$nav_menu);
    &make_nav_end  (\$r, \$nav_end);

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
# Make menu HTML
#############################################################
sub make_menu {

    my ($req_ref, $file) = @_;
    my %menu =();
    my @menu_order;
    my $html;
    
    my $oldRS = $/;
    $/ = "";
    open (MENU, $file) or $$req_ref->log->error("Can't open nav menu conf file [$file]: $!");
    my @menus = <MENU>;
    
    
    foreach my $m (@menus){
        chomp ($m);
        next if (/^#/);
        
        my @this_menu = reverse(split (/\n/,$m));
        my $header = pop(@this_menu);
        my ($h,$u) = split (/;/,$header);
        
        $html .=<<EOH;
<TR>
 <TD>&nbsp;</TD>
 <TD CLASS="navbarhead"><A HREF="$u">$h</A></TD>
 <TD>&nbsp;</TD>
</TR>
EOH

        foreach my $item (reverse(@this_menu)){
            my ($h,$u) = split (/;/,$item);
            $html .=<<EOI;
<TR>
 <TD>&nbsp;</TD>
 <TD CLASS="navbar"><A HREF="$u">$h</A>&nbsp;&nbsp;&nbsp;</TD>
 <TD>&nbsp;</TD>
</TR>
EOI

        } # end of foreach my $item

        $html .=<<EOH;
<TR>
 <TD COLSPAN="3"><IMG HEIGHT="10" WIDTH="1" SRC="/icons/nothing.gif" ALT=""></TD>
</TR>
EOH

    } # end of foreach my $m

    $/ = $oldRS;

    close (MENU) or $$req_ref->log->error("Can't close nav menu conf file [$file]: $!");

    return ($html)
    
} # end of sub

#############################################################
# Construct the nav menu
#############################################################
sub make_nav_menu {

my ($req_ref, $navmenu_ref) = @_;

    my $default_filename = '/def_nav.conf';
    my $local_filename = 'nav.conf';
    my $menufile;
    my %nav = ();
    my $uri = $$req_ref->filename;
    
    $uri =~ s/\w+\.*\w+$//;
    if (-e "$uri$local_filename"){
        $menufile = "$uri$local_filename";
        $$navmenu_ref = &make_menu($req_ref, $menufile);
    }
    else{
        $menufile = $$req_ref->document_root.$default_filename;
        $$navmenu_ref = &make_menu($req_ref, $menufile);
    }

    
} # end of sub

#############################################################
# Construct the nav bar top
#############################################################
sub make_nav_start {

my ($req_ref, $start_ref, $menu_ref) = @_;

$$start_ref=<<EOS;

<TABLE BORDER="0" CELLPADDING="0" CELLSPACING="0">
    <TR>
        <TD VALIGN="TOP" ROWSPAN=2 BGCOLOR="#EFEFFF">
            <!-- table cell for left navbar -->
            <!-- navbar items begin here -->
            <TABLE WIDTH="90" BORDER="0" CELLPADDING="0" CELLSPACING="3">

$$menu_ref               

            </TABLE>
            
            <!-- navbar items end here -->
            <!-- end of table cell for left navbar -->
        </TD>

        <TD ROWSPAN="2" WIDTH="30">&nbsp;</TD>
        <!-- margin between navbar and main-page -->
        <TD WIDTH="80%"><IMG WIDTH="1" HEIGHT="10" SRC="/icons/nothing.gif" ALT=""></TD>
        <!-- cell just above main page to keep distance from top header -->
        <TD ROWSPAN="2" WIDTH="5%">&nbsp;</TD>
        <!-- right margin -->
    </TR>

    <TR>
        <TD ALIGN="LEFT" VALIGN="TOP">

        <!-- table cell for main page -->
<!-- ---------------page content starts here--------------- -->
        
EOS

} # end of sub"

#############################################################
# Construct the nav bar top
#############################################################
sub make_nav_end {

my ($req_ref, $end_ref) = @_;

$$end_ref=<<EOS;

<!-- ---------------page content ends here--------------- -->
        <!-- end of table cell for main page -->
        </TD>
    </TR>

</TABLE>

<!-- close table for page content -->
EOS

}# end of sub"

#############################################################
# Construct the navigation trail
#############################################################
sub make_nav_trail {

    my ($req_ref, $trail_ref) = @_;
    
    my %map = ();
    my %trail =();
    my $uri = $$req_ref->uri;
    my @dirs = split (/\//,$uri);
    my $mapfile = $$req_ref->document_root."/urlmappings.txt";
    pop(@dirs);             # remove the html filename at the end of the URL
    shift(@dirs);           # remove the root dir

    unless (open (MAP, "$mapfile")) {
        $$req_ref->log->error("Cannot open URL map file: $!");
        return;
    }
    my @lines = <MAP>;
    foreach my $l (@lines){
        chomp ($l);
        my ($url, $label) = split (/ /,$l, 2);
        $url =~ s/^\///;    # remove the leading slash from URL mapping
        $map{$url} = $label;        
    }

    my $path_tmp ="";
    $$trail_ref = "<B CLASS=\"trailbar\">You are here:</B> <A HREF=\"/\">Home</A>";
    
    for (my $i=0;$i<(scalar(@dirs));$i++){
        $path_tmp .= "$dirs[$i]/";
        if ($i < $#dirs){
            if ($map{$path_tmp}){
                $$trail_ref .= " -> <A HREF=\"/$path_tmp\">$map{$path_tmp}</A>";
            }    
            else{
                $$trail_ref .= " -> <A HREF=\"/$path_tmp\">$dirs[$i]</A>";
            }
        }
        else{
            if ($map{$path_tmp}){
                $$trail_ref .= " -> $map{$path_tmp}";
            }    
            else{
                $$trail_ref .= " -> $dirs[$i]";
            }
        }
        
    }
    
    $$trail_ref=<<EOS; 
<TABLE WIDTH="100%" BORDER="0">
 <TR>
  <TD CLASS="trailbar">
    $$trail_ref
  </TD>
 </TR>
</TABLE>
<P>
EOS


} # end of sub"

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
