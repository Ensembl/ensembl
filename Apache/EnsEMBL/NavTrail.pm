package Apache::EnsEMBL::NavTrail;
       
use strict;

use vars qw/@ISA @EXPORT @EXPORT_OK/;

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


######################### This is the list of exported subroutines #######################

@EXPORT = qw(make_nav_trail);

##########################################################################

1;
__END__

#
# EnsEMBL module for Apache::EnsEMBL::NavTrail
#
# Cared for by Tony Cox <avc@sanger.ac.uk>
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs after the code

=head1 NAME

Apache::EnsEMBL::NavTrail - Apache Mod_perl module to generate 
a navigation trail for an EnsEMBL HTML page

=head1 SYNOPSIS

=head2 General

This mod_perl module is used to take a standard EnsEMBL
HTML page footer


=head1 RELATED MODULES

See also: Apache::EnsEMBL::SendPage, Apache::EnsEMBL::Header, Apache::EnsEMBL::Footer,
Apache::EnsEMBL::NavBar

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
