package Apache::EnsEMBL::SendPage;
       
use strict;
use Apache::Constants qw(:response :methods :http);
use Apache::File ();
use Apache::Log ();

sub handler {
    my $r = shift;
    if ((my $rc = $r->discard_request_body) != OK) {
      return $rc;
    }

    #$r->log->error("SendPage.pm handling: ", $r->filename);

    if ($r->method_number == M_INVALID) {
       $r->log->error("Invalid method in request ", $r->the_request);
       return NOT_IMPLEMENTED;
    }

    if ($r->method_number == M_OPTIONS) {
       return DECLINED;                     #http_core.c:default_handler() will pick this up
    }

    if ($r->method_number == M_PUT) {
       return HTTP_METHOD_NOT_ALLOWED;
    }

    if (-d $r->finfo) {
      #$r->log->error("Directory request: ", $r->filename);
      #printf "%s is a directory\n", $r->filename;
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
      $r->log->error("file permissions deny server access: ", 
                     $r->filename);
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

    unless ($r->header_only) {
      $r->send_fd($fh);
    }

    close $fh;
    return OK;
}

1;

__END__

#
# EnsEMBL module for Apache::EnsEMBL::SendPage
#
# Cared for by Tony Cox <avc@sanger.ac.uk>
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Apache::EnsEMBL::SendPage - Apache Mod_perl module to serve undecorated HTML files

=head1 SYNOPSIS

=head2 General

This mod_perl module is used to "echo" HTML files directly from the server document
tree. It makes no modifications to the file but simply sends it "as-is". You are
responsible for making sure the HTML is standards-compliant.

=head1 DESCRIPTION

This module takes standard HTTP GET requests received by Apache and locates the
corresponding HTML file which is then returned unmodified. It is fully HTTP/1.1
compliant. 

This module includes an optimization which checks whether the browser already
has the page in its disk cache and is only requesting the page if it has changed.

=head1 RELATED MODULES

See also: Apache::EnsEMBL::SendDecPage

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
