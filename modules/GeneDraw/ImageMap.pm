
=head1 NAME

ImageMap.pm

=head1 SYNOPSIS


&ImageMap::print_map_start(\*HTML);

&ImageMap::print_map($im_fh,$x_start,$y_start,$x_end,$y_end,$name,$url,$img_map_fh);

&ImageMap::print_map_end(\*HTML);



=head1 DESCRIPTION

ImageMap contains routines that print an image map. ImageMap.pm consists of the following modules: 
 
    print_map_start   : prints the beginning of an image map
    print_map         : prints an image map   
    print_map_end     : prints the end of an image map
   



=head1 AUTHOR

Arek Kasprzyk
arek@ebi.ac.uk

=cut


package ImageMap;

use strict;


=head1 Methods 

Methods available in the ImageMap package:



=head2 print_map_start

 Title   : print map start
 Usage   : &ImageMap::print_map_start(\*HTML)
 Function: prints start of an image map
 Example : 
 Returns : 


=cut


    
    sub print_map_start{
	my ($img_map_fh)=@_;
	print $img_map_fh "<IMG SRC=\"output.gif\" BORDER=\"0\" USEMAP=\"#Ngs22\"><BR>\n";
	print $img_map_fh "<MAP Name=\"Ngs22\">\n";
    }






=head2 print_map

 Title   : print map
 Usage   : &ImageMap::print_map($im_fh,$x_start,$y_start,$x_end,$y_end,$name,$url,$img_map_fh);
 Function: prints an image map
 Example : 
 Returns : 


=cut

    
    sub print_map{
	my($im_fh,$x_start,$y_start,$x_end,$y_end,$name,$url,$img_map_fh)=@_;
	print $img_map_fh "<AREA Shape=\"Rect\" coords = \"",$x_start," ", $y_start," ",$x_end," ",$y_end," \"  HREF=\"",$url,$name,"\">","\n";
    }




=head2 print_map_end

 Title   : print map end
 Usage   : &ImageMap::print_map_end(\*HTML)
 Function: prints end of an image map
 Example : 
 Returns : 


=cut


    sub print_map_end{
	my ($img_map_fh)=@_;
	print  $img_map_fh "</MAP>\n";
    }


1;
