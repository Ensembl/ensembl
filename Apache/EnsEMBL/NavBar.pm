package Apache::EnsEMBL::NavBar;
       
use strict;

use vars qw/@ISA @EXPORT @EXPORT_OK/;

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

######################### This is the list of exported subroutines #######################

@EXPORT = qw(make_menu make_nav_menu make_nav_start make_nav_end);

##########################################################################

1;
__END__

#
# EnsEMBL module for Apache::EnsEMBL::NavBar
#
# Cared for by Tony Cox <avc@sanger.ac.uk>
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs after the code

=head1 NAME

Apache::EnsEMBL::NavBar - Apache Mod_perl module to generate an EnsEMBL
HTML navigation sidebar

=head1 SYNOPSIS

=head2 General

This mod_perl module is used to generate a standard EnsEMBL
HTML navigation sidebar


=head1 RELATED MODULES

See also: Apache::EnsEMBL::SendPage, Apache::EnsEMBL::Footer, Apache::EnsEMBL::Header

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
