package Apache::EnsEMBL::SendDecPage;
       
use strict;
use Apache::Constants qw(:response :methods :http);
use Apache::File ();
use Apache::Log ();

sub handler {
    my $r = shift;
    my ($header, $footer);
    
    if ((my $rc = $r->discard_request_body) != OK) {
        return $rc;
    }
    
    if ($r->method_number == M_INVALID) {
        $r->log->error("Invalid method in request ", $r->the_request);
        return NOT_IMPLEMENTED;
    }

    if ($r->method_number == M_OPTIONS) {
        return DECLINED;                     #pass to http core default_handler()
    }

    if ($r->method_number == M_PUT) {
        return HTTP_METHOD_NOT_ALLOWED;
    }

    if (-d $r->finfo) {
        # Need to add a method here to make sure dir parsing is solid
        # ie. directory path: /test -> /test/index.html
        return DECLINED;                    #pass to http core default_handler()
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

    $r->set_content_length;

    if($ENV{PERL_SEND_HEADER}) {
        print "Content-type: text/html\n\n";
    }
    else {
        $r->content_type('text/html');
        $r->send_http_header;
    }
    
    &make_ensembl_header(\$r, \$header);
    &make_ensembl_footer(\$r, \$footer);
    my ($title, $gif);
    
    unless ($r->header_only) {
        $r->print($header);
        while (<$fh>){
            $r->print($_);
        }
        $r->print($footer);
    }

    
    close $fh;
    return OK;
}

#############################################################
# Construct the standard EnsEMBL page header
#############################################################

sub make_ensembl_header {

    my ($req_ref, $header_ref) = @_;

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
   <IMG SRC="/icons/humembl.gif" ALIGN="center" BORDER="0">
   </A>
  </TD>
 </TR>
</TABLE>
<P>
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
<HR>
<TABLE BORDER="0" WIDTH="100%">
  <TR>
    <TD ALIGN="LEFT" HALIGN="TOP"><I>
         last modified : $modtime </I>
    </TD>
    <TD ALIGN="RIGHT" HALIGN="TOP">
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

=head1 DESCRIPTION (most not yet implemented)

This module takes standard HTTP GET requests received by Apache and locates the
corresponding HTML file which is then returned wrapped suitably. It is fully 
HTTP/1.1 compliant. 

The header includes a "bread crumb" trail to show (roughly) where the page lies in
the server document tree structure and also a reference to a CSS which will be 
used later. The footer shows last modification date of the HTML file and a 
contact email address.

The navigation bar is created dynamically by incorporating the contents of a 
"nav.conf" file in the same directory as the HTML file (if it exists). This contains
a list of navigation options that will be written to the side bar. 
If the config file
does not exist the default navigation bar is used (read from /def_nav.conf) 

This module includes an optimization which checks whether the browser already
has the page in its disk cache and is only requesting the page if it has changed.

=head1 RELATED MODULES

See also: Apache::EnsEMBL::SendDecPage, nav.conf, def_nav.conf

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
