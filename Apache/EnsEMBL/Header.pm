package Apache::EnsEMBL::Header;
       
use strict;

use vars qw/@ISA @EXPORT @EXPORT_OK/;

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


######################### This is the list of exported subroutines #######################

@EXPORT = qw(make_ensembl_header);

##########################################################################

1;
__END__

#
# EnsEMBL module for Apache::EnsEMBL::Header
#
# Cared for by Tony Cox <avc@sanger.ac.uk>
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs after the code

=head1 NAME

Apache::EnsEMBL::Header - Apache Mod_perl module to generate a standard EnsEMBL
HTML page header

=head1 SYNOPSIS

=head2 General

This mod_perl module is used to generate a standard EnsEMBL
HTML page header


=head1 RELATED MODULES

See also: Apache::EnsEMBL::SendPage, Apache::EnsEMBL::Footer

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
