package Apache::EnsEMBL::SendDecPage;
       
use strict;
use Apache::Constants qw(:response :methods :http);
use Apache::File ();
use Apache::Log ();

#############################################################
# Mod_perl request handler for Location /test2
#############################################################
sub handler {
    my $r = shift;
    my ($header, $footer);

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
    
    my ($title, $gif);  # currently unused
    my ($trail, $nav_start, $nav_end, $nav_menu);
    
    &make_nav_trail(\$r, \$trail);
    &make_nav_menu (\$r, \$nav_menu);
    #$r->log->error("Printed Nav: $nav_menu");
    &make_nav_start(\$r, \$nav_start, \$nav_menu);
    &make_nav_end  (\$r, \$nav_end);

    #$r->log->error("Top Nav: $nav_start");
    #$r->log->error("Bottom Nav: $nav_end");

    &make_ensembl_header(\$r, \$header);
    &make_ensembl_footer(\$r, \$footer);
    
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
    #$$req_ref->log->error("Opened $file!");
    my @menus = <MENU>;
    
    
    foreach my $m (@menus){
        chomp ($m);
        next if (/^#/);
        #$$req_ref->log->error("Full Menu: $m");
        
        my @this_menu = reverse(split (/\n/,$m));
        #foreach my $tm(@this_menu){
        #    $$req_ref->log->error("\tMenu items: $tm");
        #}
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
    



    #$$req_ref->log->error("NAVBAR HTML: $html");
    
    $/ = $oldRS;

    close (MENU) or $$req_ref->log->error("Can't close nav menu conf file [$file]: $!");

    return ($html)
    
} # end of sub

#############################################################
# Construct the nav menu bar
#############################################################
sub make_nav_menu {

my ($req_ref, $navmenu_ref) = @_;

    my $default_filename = '/def_nav.conf';
    my $local_filename = 'nav.conf';
    my $menufile;
    my %nav = ();
    my $uri = $$req_ref->filename;
    
    #$$req_ref->log->error("Full URL: $uri");
    $uri =~ s/\w+\.*\w+$//;
    #$$req_ref->log->error("DIR: $uri");
    #$$req_ref->log->error("Looking for: $uri$local_filename");
    if (-e "$uri$local_filename"){
        $menufile = "$uri$local_filename";
        #$$req_ref->log->error("Menufile: $menufile");
        $$navmenu_ref = &make_menu($req_ref, $menufile);
    }
    else{
        $menufile = $$req_ref->document_root.$default_filename;
        #$$req_ref->log->error("Menufile: $menufile");
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

    #$$req_ref->log->error("Mapfile: $mapfile");
    #$$req_ref->log->error("URL dirs: @dirs");

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

    #foreach my $key (keys %map){
    #    $$req_ref->log->error("HASH: $key -> $map{$key}")
    #}
    
    my $path_tmp ="";
    $$trail_ref = "<B CLASS=\"trailbar\">You are here:</B> <A HREF=\"/\">Home</A>";
    
    for (my $i=0;$i<(scalar(@dirs));$i++){
        $path_tmp .= "$dirs[$i]/";
        #$$req_ref->log->error("PATH: $path_tmp");
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
    
    #$$req_ref->log->error("Trail: $$trail_ref");
     
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
# Construct the standard EnsEMBL page header
#############################################################
sub make_ensembl_header {

    my ($req_ref, $header_ref) = @_;

# need to be able to change TITLE and header gif

    $$header_ref=<<EOS;
<!DOCTYPE HTML PUBLIC "-//IETF//DTD HTML//EN">
<HTML>
 <HEAD>
  <TITLE>EnsEMBL Genome Views Server</TITLE>    
  <LINK REL="stylesheet" HREF="/EnsEMBL.css">
 </HEAD>
<BODY TEXT="#000000" BGCOLOR="#FFFFFF">

<TABLE WIDTH="100%" BORDER="0">
 <TR>
  <TD>
   <A HREF="http://ensembl.ebi.ac.uk/">
   <IMG SRC="/icons/humembl.gif" ALIGN="CENTER" BORDER="0">
   </A>
  </TD>
 </TR>
</TABLE>
<!P>
EOS

    return;

} # end of sub"

#############################################################
# Construct the standard EnsEMBL page header
#############################################################
sub make_ensembl_footer {

    my ($req_ref, $footer_ref) = @_;
    my $modtime = localtime (time - $$req_ref->mtime);
    
    $$footer_ref=<<EOS;
<P>
<TABLE WIDTH="100%" BORDER="0">
 <TR>
  <TD CLASS="trailbar">
    &nbsp;
  </TD>
 </TR>
</TABLE>
<HR>
<!P>
<TABLE BORDER="0" WIDTH="100%">
  <TR>
    <TD CLASS="lfooter"><I>
         last modified : $modtime </I>
    </TD>
    <TD CLASS="rfooter" ALIGN="RIGHT">
        <I><A HREF=mailto:webmaster\@sanger.ac.uk>webmaster\@sanger.ac.uk</A></I>
    </TD>
  </TR>
</TABLE>
</BODY>
</HTML>
EOS

    return;

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

Apache::EnsEMBL::SendDecPage - Apache Mod_perl module to serve decorated HTML files

=head1 SYNOPSIS

=head2 General

This mod_perl module is used to take a skeleton HTML file and wrap it in a header,
footer and side navigation bar.

=head1 DESCRIPTION (some features not yet implemented)

This module takes standard HTTP GET requests received by Apache and locates the
corresponding HTML file which is then returned wrapped suitably. It is fully 
HTTP/1.1 compliant. 

The header includes a "bread crumb" trail to show (roughly) where the current page 
lies in the server document tree structure and a reference to a CSS which is 
used to apply a style. 

The navigation bar is created dynamically by incorporating the contents of a 
"nav.conf" file in the same directory as the HTML file (if it exists). This contains
a list of navigation options that will be written to the side bar. 
If the config file
does not exist the default navigation bar is used (/def_nav.conf).

The footer shows last modification date of the HTML file and a contact email address.

This module includes an optimization which checks whether the browser already
has the page in its disk cache and is only requesting the page if it has changed.

=head1 RELATED MODULES

See also: Apache::EnsEMBL::SendPage, nav.conf, def_nav.conf

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
