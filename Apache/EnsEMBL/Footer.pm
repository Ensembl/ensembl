package Apache::EnsEMBL::Footer;
       
use strict;

use vars qw/@ISA @EXPORT @EXPORT_OK/;

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

######################### This is the list of exported subroutines #######################

@EXPORT = qw(make_ensembl_footer);

##########################################################################

1;
__END__

#
# EnsEMBL module for Apache::EnsEMBL::Footer
#
# Cared for by Tony Cox <avc@sanger.ac.uk>
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs after the code

=head1 NAME

Apache::EnsEMBL::Footer - Apache Mod_perl module to generate a standard EnsEMBL
HTML page footer

=head1 SYNOPSIS

=head2 General

This mod_perl module is used to take a standard EnsEMBL
HTML page footer


=head1 RELATED MODULES

See also: Apache::EnsEMBL::SendPage, Apache::EnsEMBL::Header

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
