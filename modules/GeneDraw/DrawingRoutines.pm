

=head1 NAME

DrawingRoutines.pm

=head1 SYNOPSIS


&DrawingRoutines::draw__s2_band($im_fh,$x_start,$y_start,$x_end,$y_end,$black,$black,$name,$url,$img_map_fh);

&DrawingRoutines::draw_a_rectangle($im_fh,$x_start,$y_start,$x_end,$y_end,$lightgrey,$black,$name,$url,$img_map_fh);

&DrawingRoutines::draw_a_centromere_p($im_fh,$x_start,$y_start,$x_end,$y_end,$darkgrey,$black,$name,$url,$img_map_fh);

&DrawingRoutines::draw_a_centromere_q($im_fh,$x_start,$y_start,$x_end,$y_end,$darkgrey,$black,$name,$url,$img_map_fh);



=head1 DESCRIPTION

DrawingRoutines contains basic routines for drawing triangles and rectangles. DrawingRoutines consists of the following modules: 
 

    draw_s2_band        : draws two rectangles
    draw_a_rectangle    : draws a filled rectangle
    draw_a_centromere_p : draws a triangle pointing down
    draw_a_centromere_q : draws a triangle pointing up


=head1 AUTHOR

Arek Kasprzyk
arek@ebi.ac.uk

=cut



package DrawingRoutines;

use strict;
use GD;




=head1 Methods 

Methods available in the DrawingRoutines package:




=head2 draw_a_rectangle

 Title   : draw a rectangle
 Usage   : &DrawingRoutines::draw_a_rectangle($im_fh,$x_start,$y_start,$x_end,$y_end,$lightgrey,$black,$name,$url,$img_map_fh);
 Function: draws a rectangle
 Example : 
 Returns : 


=cut


sub draw_a_rectangle {    
    my($im_fh,$x_start,$y_start,$x_end,$y_end,$colour,$font_colour,$name,$url,$img_map_fh)=@_;
    $im_fh->filledRectangle($x_start,$y_start,$x_end,$y_end,$colour);
 
}




=head2 draw_a_s2_band

 Title   : draw a s2 band
 Usage   : &DrawingRoutines::draw__s2_band($im_fh,$x_start,$y_start,$x_end,$y_end,$black,$black,$name,$url,$img_map_fh);
 Function: draws two rectangles
 Example : 
 Returns : 


=cut



sub draw__s2_band {
    my($im_fh,$x_start,$y_start,$x_end,$y_end,$colour,$font_colour,$name,$url,$img_map_fh)=@_;    
    $im_fh->filledRectangle($x_start,$y_start+($y_end-$y_start)/5,$x_end,$y_start+(($y_end-$y_start)/5)*2,$colour);
    $im_fh->filledRectangle($x_start,$y_start+(($y_end-$y_start)/5)*3,$x_end,$y_start+(($y_end-$y_start)/5)*4,$colour);
  
}




=head2 draw_a_centromere_P

 Title   : draw a centromere_p
 Usage   :&DrawingRoutines::draw_a_centromere_p($im_fh,$x_start,$y_start,$x_end,$y_end,$darkgrey,$black,$name,$url,$img_map_fh);
 Function: draws a triangle pointing down
 Example : 
 Returns : 


=cut




sub draw_a_centromere_p { 
    my($im_fh,$x_start,$y_start,$x_end,$y_end,$colour,$font_colour,$name,$url,$img_map_fh)=@_; 
    my $poly = new GD::Polygon;
    $poly->addPt($x_start,$y_start);
    $poly->addPt($x_end,$y_start);
    $poly->addPt(($x_end-$x_start)*0.5+$x_start,$y_end);
    $im_fh->filledPolygon($poly,$colour);
}




=head2 draw_a_centromere_q

 Title   : draw a centromere_q
 Usage   :&DrawingRoutines::draw_a_centromere_q($im_fh,$x_start,$y_start,$x_end,$y_end,$darkgrey,$black,$name,$url,$img_map_fh);
 Function: draws a triangle pointing up
 Example : 
 Returns : 


=cut



sub draw_a_centromere_q { 
    my($im_fh,$x_start,$y_start,$x_end,$y_end,$colour,$font_colour,$name,$url,$img_map_fh)=@_; 			  
    my $poly = new GD::Polygon;
    $poly->addPt($x_start,$y_end);
    $poly->addPt($x_end,$y_end);
    $poly->addPt(($x_end-$x_start)*0.5+$x_start,$y_start);
    $im_fh->filledPolygon($poly,$colour);
   
}



1;





